/// \class TRestDetectorAvalancheProcess
/// Applies the avalanche gain of the amplification gap to the electrons drifted by
/// TRestDetectorHitsDriftProcess. The amplification field is E = amplificationVoltage / gap (uniform,
/// parallel plates); the mean gain G and the Polya parameter theta at E come from the gain table of
/// the TRestDetectorMedium with that gap at the pressure and temperature of the medium (`table`, its key,
/// only if there are several; a field out of the range of the table is an error). The gain of one electron
/// follows a Polya distribution; a hit with energy E_hit carries n = E_hit / W electrons (more than one with
/// maxHits in hitsDrift), whose total gain, a sum of n Polya, is drawn from a Gamma(n (theta + 1), G / (theta
/// + 1)). The energy of the hit is multiplied by its gain per electron (total gain / n), as legacy multiplied
/// it by the gain. G and theta used are written to the metadata of the process.
///
/// In a field map (special simulations, after TRestDetectorHitsDriftProcess with a drift response),
/// `responseFile` is an amplification response computed by restGarfield in the field map of a unit cell
/// (TRestDetectorAmplificationResponse, the `amplification` key of a TRestDetectorTransportManager) instead
/// of amplificationVoltage and gap: of the n = E_hit / W electrons of a hit (rounded, at least 1), each one
/// is amplified with the transparency of the response (the others are lost on the mesh or a dielectric) and
/// its gain is sampled from the histogram of the avalanche size; the energy of the hit is multiplied by the
/// total gain over n, and its time is delayed by a time sampled from the response (the time from the start
/// plane of the response, the end of the drift, to the anode). A hit with no electron amplified is
/// removed. The mean gain, theta and transparency of the response are written to the metadata.
///
/// Differences with legacy: the gain comes from the avalanche simulations of restGarfield at the voltage,
/// with Polya fluctuations, instead of a fixed detectorGain with a gaussian resolution
/// (energyReference, resolutionReference: that is TRestDetectorHitsSmearingProcess, not to be used with
/// this process).
///
/// Parameters: amplificationVoltage and gap, or responseFile; medium (empty: the only one with `charge`),
/// table (""), seed (0:
/// random seed; otherwise the generator is re-seeded at each event from the seed and the event and subevent
/// IDs, TRestEvent::EventSeed, so the result does not depend on the number of threads).
///
/// Example:
///
///    avalanche:
///      class: "TRestDetectorAvalancheProcess"
///      amplificationVoltage: 330 V
///      gap: 128 um
///
/// In a field map:
///
///    avalanche:
///      class: "TRestDetectorAvalancheProcess"
///      responseFile: "amplification_320V.root"

#include "TRestDetectorAvalancheProcess.h"

#include <algorithm>
#include <cmath>
#include <stdexcept>

#include "TFile.h"
#include "TRestDetectorHitsEvent.h"
#include "TRestDetectorMedium.h"
#include "TRestDetectorTransportManager.h"
#include "TRestManager.h"

namespace {
/// Gamma(shape, scale) drawn from a TRandom (Marsaglia and Tsang). Not std::gamma_distribution, whose
/// algorithm depends on the standard library (the same seed would give other numbers on another platform)
double Gamma(TRandom& random, double shape, double scale) {
    if (!(shape > 0) || !(scale > 0)) return 0;
    if (shape < 1) {
        // Gamma(shape + 1) U^(1 / shape)
        const double gamma = Gamma(random, shape + 1, scale);
        return gamma * std::pow(random.Rndm(), 1 / shape);
    }
    const double d = shape - 1. / 3;
    const double c = 1 / std::sqrt(9 * d);
    while (true) {
        const double x = random.Gaus(0, 1);
        double v = 1 + c * x;
        if (v <= 0) continue;
        v = v * v * v;
        const double u = random.Rndm();
        if (u < 1 - 0.0331 * x * x * x * x || std::log(u) < 0.5 * x * x + d * (1 - v + std::log(v))) {
            return d * v * scale;
        }
    }
}
}  // namespace

static const bool TRestDetectorAvalancheProcess_FieldsRegistered = []() {
    auto& reg = TRestMetadataFieldRegistry::Instance();
    using P = TRestDetectorAvalancheProcess;
    reg.RegisterField<P>("amplificationVoltage", &P::fAmplificationVoltage);
    reg.RegisterField<P>("gap", &P::fGap);
    reg.RegisterField<P>("medium", &P::fMediumName);
    reg.RegisterField<P>("table", &P::fTable);
    reg.RegisterField<P>("responseFile", &P::fResponseFile);
    reg.RegisterField<P>("seed", &P::fSeed);
    reg.RegisterField<P>("gain", &P::fGain);
    reg.RegisterField<P>("theta", &P::fTheta);
    reg.RegisterField<P>("transparency", &P::fTransparency);
    return true;
}();

namespace {
const bool kRegistered = []() {
    MetadataClassRegistry::Instance().Register(
        "TRestDetectorAvalancheProcess", [](const std::string& instanceName, const YAML::Node& params) {
            return std::make_unique<TRestDetectorAvalancheProcess>(instanceName, params);
        });
    return true;
}();
}  // namespace

TRestDetectorAvalancheProcess::TRestDetectorAvalancheProcess() : TRestEventProcess() {
    fName = "TRestDetectorAvalancheProcess";
}

TRestDetectorAvalancheProcess::TRestDetectorAvalancheProcess(const std::string& instanceName,
                                                             const YAML::Node& node)
    : TRestEventProcess(instanceName, node) {
    LoadConfig();
}

TRestDetectorAvalancheProcess::TRestDetectorAvalancheProcess(const std::string& fileName,
                                                             const std::string& sectionName)
    : TRestEventProcess(fileName, sectionName) {
    LoadConfig();
}

void TRestDetectorAvalancheProcess::LoadConfig() {
    TRestEventProcess::LoadConfig();
    if (!fNode || fNode.IsNull()) {
        RESTError << "TRestDetectorAvalancheProcess::LoadConfig - YAML node is missing" << RESTendl;
        return;
    }
    UpdateParamsFromYAML<TRestDetectorAvalancheProcess>(fNode);
    UpdateYAMLFromParams<TRestDetectorAvalancheProcess>(fNode);
}

void TRestDetectorAvalancheProcess::InitProcess() {
    const std::string user = "TRestDetectorAvalancheProcess";
    if (!fResponseFile.empty() && (fAmplificationVoltage != 0 || fGap != 0)) {
        throw std::runtime_error(user + ": responseFile or amplificationVoltage and gap, not both");
    }
    if (fResponseFile.empty() && (fAmplificationVoltage <= 0 || fGap <= 0)) {
        throw std::runtime_error(user + ": amplificationVoltage and gap must be > 0 (or responseFile)");
    }
    if (!fManager) throw std::runtime_error(user + ": no TRestManager");
    const auto* medium = TRestDetectorMedium::Find(*fManager, fMediumName, "charge", user);
    if (!medium) throw std::runtime_error(user + ": needs a TRestDetectorMedium with charge");
    fCharge = &medium->GetChargeResponse(user);
    if (fCharge->GetWValue() <= 0) {
        throw std::runtime_error(user + ": the charge of the medium '" + medium->GetName() +
                                 "' needs wValue > 0");
    }

    // Seed 0: random seed. Otherwise re-seeded at each event in ProcessEvent
    fRandom.SetSeed(fSeed);
    fMediumName = medium->GetName();

    if (!fResponseFile.empty()) {
        // The amplification in a field map: the response instead of the gain table
        std::unique_ptr<TFile> file(TFile::Open(fResponseFile.c_str(), "READ"));
        if (!file || file->IsZombie()) throw std::runtime_error(user + ": cannot open " + fResponseFile);
        fResponse = std::make_shared<const TRestDetectorAmplificationResponse>(
            TRestDetectorAmplificationResponse::Read(file.get()));
        if (fResponse->fCollected <= 0) {
            throw std::runtime_error(user + ": the amplification response of " + fResponseFile +
                                     " has no avalanche at the anode");
        }
        if (!fResponse->fMedium.empty() && fResponse->fMedium != fMediumName) {
            RESTWarning << user << ": the amplification response of " << fResponseFile
                        << " is for the medium '" << fResponse->fMedium << "', not '" << fMediumName << "'"
                        << RESTendl;
        }
        fGain = fResponse->fMeanSize;
        fTheta = fResponse->GetTheta();
        fTransparency = fResponse->GetTransparency();
        UpdateYAMLFromParams<TRestDetectorAvalancheProcess>(fNode);
        return;
    }
    fResponse.reset();
    fTransparency = 1;

    const auto& table = fCharge->GetGainTable(fGap, fTable);
    const double field = fAmplificationVoltage / fGap;  // V/mm
    fGain = table.GetGain(field);
    fTheta = table.GetTheta(field);
    if (fGain <= 0 || fTheta <= -1) {
        throw std::runtime_error(user + ": table " + table.GetKey() + " gives gain " + std::to_string(fGain) +
                                 " and theta " + std::to_string(fTheta) +
                                 " (gain > 0 and theta > -1 needed)");
    }
    fTable = table.GetKey();
    // The stored metadata keeps the medium, the table and the gain used
    UpdateYAMLFromParams<TRestDetectorAvalancheProcess>(fNode);
}

bool TRestDetectorAvalancheProcess::ProcessEvent(const TRestEvent& input, TRestEvent& output) {
    const auto* inputEvent = dynamic_cast<const TRestDetectorHitsEvent*>(&input);
    auto* outputEvent = dynamic_cast<TRestDetectorHitsEvent*>(&output);
    if (!inputEvent || !outputEvent) {
        throw std::runtime_error(
            "TRestDetectorAvalancheProcess: input/output is not a TRestDetectorHitsEvent");
    }

    outputEvent->Initialize();
    outputEvent->SetEventInfo(inputEvent);
    if (fSeed != 0) fRandom.SetSeed(TRestEvent::EventSeed(fSeed, *inputEvent));

    const double wValue = fCharge->GetWValue();
    for (size_t n = 0; n < inputEvent->GetNumberOfHits(); n++) {
        const double energy = inputEvent->GetEnergy(n);
        if (energy <= 0) continue;
        const double electrons = energy / wValue;
        if (fResponse) {
            // Each electron amplified with the transparency, its gain and the time of the hit from the
            // response
            const long count = std::max(1L, std::lround(electrons));
            double gain = 0;
            for (long e = 0; e < count; e++) {
                if (fRandom.Uniform() < fTransparency) gain += fResponse->SampleSize(fRandom);
            }
            if (!(gain > 0)) continue;
            const double time = inputEvent->GetTime(n) + fResponse->SampleTime(fRandom);
            outputEvent->AddHit(inputEvent->GetPosition(n), energy * gain / count, time,
                                inputEvent->GetType(n), inputEvent->GetParticleID(n));
            continue;
        }
        // Sum of the Polya gains of the electrons of the hit
        const double gain = Gamma(fRandom, electrons * (fTheta + 1), fGain / (fTheta + 1));
        outputEvent->AddHit(inputEvent->GetPosition(n), energy * gain / electrons, inputEvent->GetTime(n),
                            inputEvent->GetType(n), inputEvent->GetParticleID(n));
    }

    return outputEvent->GetNumberOfHits() > 0;
}
