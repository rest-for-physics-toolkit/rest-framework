/// \class TRestDetectorHitsDriftProcess
/// Converts the energy of each hit into ionization electrons (energy / wValue of the medium, with the
/// chosen fluctuations) and diffuses each electron with the TRestDetectorMedium at its drift field
/// (driftField of its charge; gaussian diffusion and attachment, from a transport table of the medium or its
/// fixed parameters; the values used are written to the metadata of the process) and drifts it to the
/// readout. Each surviving electron becomes a hit on the readout (its diffused position projected on the
/// readout plane, TRestDetectorReadout::ProjectToReadout, where the readout gives its channel) at its arrival
/// time, the time of the original hit plus its diffused distance to the readout over the drift velocity of
/// the medium, with the energy of its electrons: the electrons of a hit carry the fluctuated number of
/// electrons times the W-value (the energy of the hit without fluctuations), so the fluctuation is kept. A
/// hit whose number of electrons fluctuates to 0 or less is removed. TRestDetectorHitsToSignalProcess then
/// gives the signals; the drift coordinate is recovered from the time by
/// TRestDetectorHitsTimeProjectionProcess (e.g. to look at the drifted electrons in 3D).
///
/// With `diffusion: false` the hits are not split into electrons nor diffused: each hit is projected whole
/// on the readout at t + d / v, d its distance to the readout, with its charge fluctuated as above and
/// attenuated by the mean attachment, exp(-attachment * d) (faster, e.g. to check the geometry and the
/// readout; maxHits is not used).
///
/// The drift distance and its direction are given by TRestDetectorReadout::GetDistanceToReadout (e.g. the
/// Micromegas plane and its normal). Hits without a distance or outside the drift volume (distance < 0 or
/// beyond TRestDetectorReadout::GetMaxDistanceToReadout, the cathode) are removed. Diffused electrons are
/// kept between the readout plane and the cathode, as legacy. A readout without maximum distance (NaN, e.g.
/// a Micromegas readout without gasThickness) has no cathode: no upper cut.
///
/// In a non-uniform field (special simulations), `responseFile` is a drift response computed by restGarfield
/// in a field map (TRestDetectorDriftResponse, the `drift` key of a TRestDetectorTransportManager) and
/// replaces the drift field of the medium (its transport is not used, only wValue and fanoFactor): each
/// electron survives with the fraction of the electrons that reach the boundary plane from the position of
/// its hit (TRestDetectorDriftResponse::GetSurvival, interpolated between the nodes of the grid); a surviving
/// electron becomes a hit on the readout (the hit projected on the readout plane plus the displacement
/// sampled from the response along the axes of the readout plane, TRestDetectorDriftResponse::PlaneAxes of
/// the readout normal) at the time of the hit plus the drift time sampled from the response
/// (SampleArrival). The arrival is on the readout plane at the time of the boundary plane: the drift from
/// the boundary to the readout (through the mesh) belongs to the amplification (TRestDetectorAvalancheProcess
/// with an amplification response). Hits outside the grid of the response are removed; the cuts of the drift
/// volume are those of the response (its grid and its field map). `diffusion: false` is not available with a
/// response.
///
/// Differences with legacy: the attachment is exp(-attachment * distance) with attachment in 1/cm (as
/// Garfield), the gas parameters come from the medium instead of the process, there are no veto planes
/// (veto hits are split beforehand) and there is no microscopic Garfield drift: in a uniform field it uses
/// the same transport table as the gaussian diffusion (microscopic simulations are done with restGarfield).
///
/// Parameters: readout (empty: the data readout of the TRestDetectorReadoutManager), medium (empty:
/// the only one with `charge`), fluctuations ("fano", "poisson" or "none"), maxHits (1000, 0: no
/// limit), diffusion (true), responseFile (empty: the uniform drift field of the medium), seed (0: random
/// seed; otherwise the generator is re-seeded at each event from the seed and the event and subevent IDs,
/// TRestEvent::EventSeed, so the result does not depend on the number of threads).
///
/// Example:
///
///    hitsDrift:
///      class: "TRestDetectorHitsDriftProcess"
///      maxHits: 1000
///
/// In a field map:
///
///    hitsDrift:
///      class: "TRestDetectorHitsDriftProcess"
///      responseFile: "drift_response.root"
///      maxHits: 1000

#include "TRestDetectorHitsDriftProcess.h"

#include <algorithm>
#include <cmath>
#include <stdexcept>

#include "TFile.h"
#include "TRestDetectorHitsEvent.h"
#include "TRestDetectorMedium.h"
#include "TRestDetectorReadout.h"
#include "TRestDetectorReadoutManager.h"
#include "TRestManager.h"
#include "TRestTools.h"

using ROOT::Math::XYZVector;

static const bool TRestDetectorHitsDriftProcess_FieldsRegistered = []() {
    auto& reg = TRestMetadataFieldRegistry::Instance();
    using P = TRestDetectorHitsDriftProcess;
    reg.RegisterField<P>("readout", &P::fReadoutName);
    reg.RegisterField<P>("medium", &P::fMediumName);
    reg.RegisterField<P>("fluctuations", &P::fFluctuations);
    reg.RegisterField<P>("maxHits", &P::fMaxHits);
    reg.RegisterField<P>("responseFile", &P::fResponseFile);
    reg.RegisterField<P>("diffusion", &P::fDiffusion);
    reg.RegisterField<P>("seed", &P::fSeed);
    reg.RegisterField<P>("longDiff", &P::fLongitudinalDiffusion);
    reg.RegisterField<P>("transDiff", &P::fTransversalDiffusion);
    reg.RegisterField<P>("attachment", &P::fAttachment);
    reg.RegisterField<P>("driftVelocity", &P::fDriftVelocity);
    return true;
}();

namespace {
const bool kRegistered = []() {
    MetadataClassRegistry::Instance().Register(
        "TRestDetectorHitsDriftProcess", [](const std::string& instanceName, const YAML::Node& params) {
            return std::make_unique<TRestDetectorHitsDriftProcess>(instanceName, params);
        });
    return true;
}();
}  // namespace

TRestDetectorHitsDriftProcess::TRestDetectorHitsDriftProcess() : TRestEventProcess() {
    fName = "TRestDetectorHitsDriftProcess";
}

TRestDetectorHitsDriftProcess::TRestDetectorHitsDriftProcess(const std::string& instanceName,
                                                             const YAML::Node& node)
    : TRestEventProcess(instanceName, node) {
    LoadConfig();
}

TRestDetectorHitsDriftProcess::TRestDetectorHitsDriftProcess(const std::string& fileName,
                                                             const std::string& sectionName)
    : TRestEventProcess(fileName, sectionName) {
    LoadConfig();
}

void TRestDetectorHitsDriftProcess::LoadConfig() {
    TRestEventProcess::LoadConfig();
    if (!fNode || fNode.IsNull()) {
        RESTError << "TRestDetectorHitsDriftProcess::LoadConfig - YAML node is missing" << RESTendl;
        return;
    }
    UpdateParamsFromYAML<TRestDetectorHitsDriftProcess>(fNode);
    fFluctuations = TRestTools::NormalizeOption("fluctuations", fFluctuations, kFluctuations);
    UpdateYAMLFromParams<TRestDetectorHitsDriftProcess>(fNode);
}

void TRestDetectorHitsDriftProcess::InitProcess() {
    if (!fManager) throw std::runtime_error("TRestDetectorHitsDriftProcess: no TRestManager");

    auto* readouts =
        dynamic_cast<TRestDetectorReadoutManager*>(fManager->GetMetadataClass("TRestDetectorReadoutManager"));
    if (!readouts) {
        throw std::runtime_error("TRestDetectorHitsDriftProcess: no TRestDetectorReadoutManager found");
    }
    fReadout = readouts->SelectReadout(fReadoutName, "TRestDetectorHitsDriftProcess");

    const std::string user = "TRestDetectorHitsDriftProcess";
    const auto* medium = TRestDetectorMedium::Find(*fManager, fMediumName, "charge", user);
    if (!medium) throw std::runtime_error(user + ": needs a TRestDetectorMedium with charge");
    fCharge = &medium->GetChargeResponse(user);
    if (fCharge->GetWValue() <= 0) {
        throw std::runtime_error(user + ": the charge of the medium '" + medium->GetName() +
                                 "' needs wValue > 0");
    }
    if (fFluctuations == "fano" && fCharge->GetFanoFactor() <= 0) {
        throw std::runtime_error(user + ": fluctuations 'fano' need fanoFactor > 0 in the charge of '" +
                                 medium->GetName() + "'");
    }
    fMediumName = medium->GetName();
    // Seed 0: random seed. Otherwise re-seeded at each event in ProcessEvent
    fRandom.SetSeed(fSeed);

    if (!fResponseFile.empty()) {
        // The drift in a field map: the response replaces the transport of the medium
        if (!fDiffusion)
            throw std::runtime_error(user + ": diffusion false is not available with responseFile");
        std::unique_ptr<TFile> file(TFile::Open(fResponseFile.c_str(), "READ"));
        if (!file || file->IsZombie()) throw std::runtime_error(user + ": cannot open " + fResponseFile);
        fResponse =
            std::make_shared<const TRestDetectorDriftResponse>(TRestDetectorDriftResponse::Read(file.get()));
        if (!fResponse->fMedium.empty() && fResponse->fMedium != fMediumName) {
            RESTWarning << user << ": the drift response of " << fResponseFile << " is for the medium '"
                        << fResponse->fMedium << "', not '" << fMediumName << "'" << RESTendl;
        }
        UpdateYAMLFromParams<TRestDetectorHitsDriftProcess>(fNode);
        return;
    }
    fResponse.reset();
    fLongitudinalDiffusion = fCharge->GetLongitudinalDiffusion();
    fTransversalDiffusion = fCharge->GetTransversalDiffusion();
    fAttachment = fCharge->GetAttachment();
    fDriftVelocity = fCharge->GetDriftVelocity();
    if (fDriftVelocity <= 0) {
        throw std::runtime_error(user + ": the charge of the medium '" + medium->GetName() +
                                 "' gives no drift velocity");
    }
    // The stored metadata keeps the medium and the transport used
    UpdateYAMLFromParams<TRestDetectorHitsDriftProcess>(fNode);
}

bool TRestDetectorHitsDriftProcess::ProcessEvent(const TRestEvent& input, TRestEvent& output) {
    const auto* inputEvent = dynamic_cast<const TRestDetectorHitsEvent*>(&input);
    auto* outputEvent = dynamic_cast<TRestDetectorHitsEvent*>(&output);
    if (!inputEvent || !outputEvent) {
        throw std::runtime_error(
            "TRestDetectorHitsDriftProcess: input/output is not a TRestDetectorHitsEvent");
    }

    outputEvent->Initialize();
    outputEvent->SetEventInfo(inputEvent);
    if (fSeed != 0) fRandom.SetSeed(TRestEvent::EventSeed(fSeed, *inputEvent));

    double totalEnergy = 0;
    for (size_t n = 0; n < inputEvent->GetNumberOfHits(); n++) {
        if (inputEvent->GetEnergy(n) > 0) totalEnergy += inputEvent->GetEnergy(n);
    }
    const double wValue = fCharge->GetWValue();
    // Effective W-value: each electron carries the energy of several if there are more than fMaxHits
    const double effectiveW =
        fMaxHits > 0 && totalEnergy / wValue > fMaxHits ? totalEnergy / fMaxHits : wValue;

    for (size_t n = 0; n < inputEvent->GetNumberOfHits(); n++) {
        const double energy = inputEvent->GetEnergy(n);
        if (energy <= 0) continue;

        const XYZVector position = inputEvent->GetPosition(n);
        // With a drift response: the hits in its grid, the readout plane at the position of the hit
        if (fResponse && !fResponse->Contains(position)) continue;
        const double distance = fReadout->GetDistanceToReadout(position);
        // NaN maximum distance: no cathode, no upper cut
        const double maxDistance = fReadout->GetMaxDistanceToReadout(position);
        const bool hasCathode = !std::isnan(maxDistance);
        if (!fResponse && (std::isnan(distance) || distance < 0 || (hasCathode && distance > maxDistance))) {
            continue;
        }
        const XYZVector normal = fReadout->GetReadoutNormal(position);
        // A readout without distance (no normal) cannot place the arrivals of the response
        if (fResponse && normal.R() == 0) continue;

        double electrons = energy / wValue;
        if (fFluctuations == "fano") {
            electrons = fRandom.Gaus(electrons, std::sqrt(fCharge->GetFanoFactor() * electrons));
        } else if (fFluctuations == "poisson") {
            electrons = fRandom.Poisson(electrons);
        }
        if (!(electrons > 0)) continue;
        if (!fDiffusion) {
            // The whole hit on the readout, its charge attenuated by the mean attachment (1/cm, d in mm)
            const double charge = electrons * wValue * std::exp(-fAttachment * distance / 10.);
            outputEvent->AddHit(fReadout->ProjectToReadout(position), charge,
                                inputEvent->GetTime(n) + distance / fDriftVelocity, inputEvent->GetType(n),
                                inputEvent->GetParticleID(n));
            continue;
        }
        // The electrons of the hit carry the fluctuated charge (electrons * W)
        const long nElectrons = std::max(1L, std::lround(electrons * wValue / effectiveW));
        const double electronEnergy = electrons * wValue / nElectrons;

        if (fResponse) {
            // Survival and arrival (displacement on the readout plane and drift time) from the response
            const double survival = fResponse->GetSurvival(position);
            const XYZVector projected = fReadout->ProjectToReadout(position);
            const auto [u, v] = TRestDetectorDriftResponse::PlaneAxes(normal);
            for (long e = 0; e < nElectrons; e++) {
                if (fRandom.Uniform() >= survival) continue;
                double du = 0, dv = 0, time = 0;
                if (!fResponse->SampleArrival(position, fRandom, du, dv, time)) continue;
                outputEvent->AddHit(projected + du * u + dv * v, electronEnergy,
                                    inputEvent->GetTime(n) + time, inputEvent->GetType(n),
                                    inputEvent->GetParticleID(n));
            }
            continue;
        }
        for (long e = 0; e < nElectrons; e++) {
            XYZVector electron = position;
            if (!TRestDetectorChargeResponse::DiffuseElectron(electron, normal, distance,
                                                              fLongitudinalDiffusion, fTransversalDiffusion,
                                                              fAttachment, fRandom)) {
                continue;
            }
            // Kept between the readout plane and the cathode, then collected on the readout plane
            double diffusedDistance = std::max(0.0, fReadout->GetDistanceToReadout(electron));
            if (hasCathode) diffusedDistance = std::min(diffusedDistance, maxDistance);
            const double arrival = inputEvent->GetTime(n) + diffusedDistance / fDriftVelocity;
            outputEvent->AddHit(fReadout->ProjectToReadout(electron), electronEnergy, arrival,
                                inputEvent->GetType(n), inputEvent->GetParticleID(n));
        }
    }

    return outputEvent->GetNumberOfHits() > 0;
}
