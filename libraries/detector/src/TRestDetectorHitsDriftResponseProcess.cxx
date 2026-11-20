/// \class TRestDetectorHitsDriftResponseProcess
/// Converts the energy of each hit into ionization electrons (energy / wValue of the medium, with the chosen
/// fluctuations, as TRestDetectorHitsDriftProcess) and drifts each electron with a drift response in a
/// non-uniform field (TRestDetectorDriftResponse, computed by restGarfield from a field map, `responseFile`)
/// to the boundary plane above the readout, where the amplification starts: for special simulations with a
/// field map, instead of TRestDetectorHitsDriftProcess (uniform field, the driftField of the medium).
///
/// Each electron survives with the fraction of the electrons that reach the boundary from the position of its
/// hit (TRestDetectorDriftResponse::GetSurvival, interpolated between the nodes of the grid); a surviving
/// electron becomes a hit on the readout (the hit projected on the readout plane,
/// TRestDetectorReadout::ProjectToReadout, plus the displacement sampled from the response along the axes of
/// the readout plane, TRestDetectorDriftResponse::PlaneAxes of the readout normal) at the time of the hit
/// plus the drift time sampled from the response (SampleArrival), with the energy of its electrons as
/// TRestDetectorHitsDriftProcess. Hits outside the grid of the response are removed, and a hit whose number
/// of electrons fluctuates to 0 or less.
///
/// The arrival is on the readout plane at the time of the boundary plane: the drift from the boundary to the
/// readout (through the mesh) belongs to the amplification.
///
/// Parameters: responseFile, readout (empty: the data readout of the TRestDetectorReadoutManager), medium
/// (empty: the only one with `charge`), fluctuations ("fano", "poisson" or "none"), maxHits (1000, 0: no
/// limit), seed (0: random seed; otherwise the generator is re-seeded at each event from the seed and the
/// event and subevent IDs, TRestDetectorHitsEvent::EventSeed, so the result does not depend on the number of
/// threads).
///
/// Example:
///
///    driftResponse:
///      class: "TRestDetectorHitsDriftResponseProcess"
///      responseFile: "drift_response.root"
///      maxHits: 1000
///
/// \author JuanAn García <juanangp@unizar.es>

#include "TRestDetectorHitsDriftResponseProcess.h"

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

static const bool TRestDetectorHitsDriftResponseProcess_FieldsRegistered = []() {
    auto& reg = TRestMetadataFieldRegistry::Instance();
    using P = TRestDetectorHitsDriftResponseProcess;
    reg.RegisterField<P>("responseFile", &P::fResponseFile);
    reg.RegisterField<P>("readout", &P::fReadoutName);
    reg.RegisterField<P>("medium", &P::fMediumName);
    reg.RegisterField<P>("fluctuations", &P::fFluctuations);
    reg.RegisterField<P>("maxHits", &P::fMaxHits);
    reg.RegisterField<P>("seed", &P::fSeed);
    return true;
}();

namespace {
const bool kRegistered = []() {
    MetadataClassRegistry::Instance().Register(
        "TRestDetectorHitsDriftResponseProcess",
        [](const std::string& instanceName, const YAML::Node& params) {
            return std::make_unique<TRestDetectorHitsDriftResponseProcess>(instanceName, params);
        });
    return true;
}();
}  // namespace

TRestDetectorHitsDriftResponseProcess::TRestDetectorHitsDriftResponseProcess() : TRestEventProcess() {
    fName = "TRestDetectorHitsDriftResponseProcess";
}

TRestDetectorHitsDriftResponseProcess::TRestDetectorHitsDriftResponseProcess(const std::string& instanceName,
                                                                             const YAML::Node& node)
    : TRestEventProcess(instanceName, node) {
    LoadConfig();
}

TRestDetectorHitsDriftResponseProcess::TRestDetectorHitsDriftResponseProcess(const std::string& fileName,
                                                                             const std::string& sectionName)
    : TRestEventProcess(fileName, sectionName) {
    LoadConfig();
}

void TRestDetectorHitsDriftResponseProcess::LoadConfig() {
    TRestEventProcess::LoadConfig();
    if (!fNode || fNode.IsNull()) {
        RESTError << "TRestDetectorHitsDriftResponseProcess::LoadConfig - YAML node is missing" << RESTendl;
        return;
    }
    UpdateParamsFromYAML<TRestDetectorHitsDriftResponseProcess>(fNode);
    fFluctuations = TRestTools::NormalizeOption("fluctuations", fFluctuations, kFluctuations);
    UpdateYAMLFromParams<TRestDetectorHitsDriftResponseProcess>(fNode);
}

void TRestDetectorHitsDriftResponseProcess::InitProcess() {
    const std::string user = "TRestDetectorHitsDriftResponseProcess";
    if (!fManager) throw std::runtime_error(user + ": no TRestManager");

    auto* readouts =
        dynamic_cast<TRestDetectorReadoutManager*>(fManager->GetMetadataClass("TRestDetectorReadoutManager"));
    if (!readouts) throw std::runtime_error(user + ": no TRestDetectorReadoutManager found");
    fReadout = readouts->SelectReadout(fReadoutName, user);

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

    if (fResponseFile.empty()) throw std::runtime_error(user + ": responseFile is needed");
    std::unique_ptr<TFile> file(TFile::Open(fResponseFile.c_str(), "READ"));
    if (!file || file->IsZombie()) throw std::runtime_error(user + ": cannot open " + fResponseFile);
    fResponse =
        std::make_shared<const TRestDetectorDriftResponse>(TRestDetectorDriftResponse::Read(file.get()));
    if (!fResponse->fMedium.empty() && fResponse->fMedium != fMediumName) {
        RESTWarning << user << ": the drift response of " << fResponseFile << " is for the medium '"
                    << fResponse->fMedium << "', not '" << fMediumName << "'" << RESTendl;
    }
    UpdateYAMLFromParams<TRestDetectorHitsDriftResponseProcess>(fNode);

    // Seed 0: random seed. Otherwise re-seeded at each event in ProcessEvent
    fRandom.SetSeed(fSeed);
}

bool TRestDetectorHitsDriftResponseProcess::ProcessEvent(const TRestEvent& input, TRestEvent& output) {
    const auto* inputEvent = dynamic_cast<const TRestDetectorHitsEvent*>(&input);
    auto* outputEvent = dynamic_cast<TRestDetectorHitsEvent*>(&output);
    if (!inputEvent || !outputEvent) {
        throw std::runtime_error(
            "TRestDetectorHitsDriftResponseProcess: input/output is not a TRestDetectorHitsEvent");
    }

    outputEvent->Initialize();
    outputEvent->SetEventInfo(inputEvent);
    if (fSeed != 0) fRandom.SetSeed(TRestDetectorHitsEvent::EventSeed(fSeed, *inputEvent));

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
        if (!fResponse->Contains(position)) continue;
        // A readout without distance (no normal) cannot place the arrivals
        const XYZVector normal = fReadout->GetReadoutNormal(position);
        if (normal.R() == 0) continue;

        double electrons = energy / wValue;
        if (fFluctuations == "fano") {
            electrons = fRandom.Gaus(electrons, std::sqrt(fCharge->GetFanoFactor() * electrons));
        } else if (fFluctuations == "poisson") {
            electrons = fRandom.Poisson(electrons);
        }
        if (!(electrons > 0)) continue;
        // The electrons of the hit carry the fluctuated charge (electrons * W)
        const long nElectrons = std::max(1L, std::lround(electrons * wValue / effectiveW));
        const double electronEnergy = electrons * wValue / nElectrons;

        const double survival = fResponse->GetSurvival(position);
        const XYZVector projected = fReadout->ProjectToReadout(position);
        const auto [u, v] = TRestDetectorDriftResponse::PlaneAxes(normal);
        for (long e = 0; e < nElectrons; e++) {
            if (fRandom.Uniform() >= survival) continue;
            double du = 0, dv = 0, time = 0;
            if (!fResponse->SampleArrival(position, fRandom, du, dv, time)) continue;
            outputEvent->AddHit(projected + du * u + dv * v, electronEnergy, inputEvent->GetTime(n) + time,
                                inputEvent->GetType(n), inputEvent->GetParticleID(n));
        }
    }

    return outputEvent->GetNumberOfHits() > 0;
}
