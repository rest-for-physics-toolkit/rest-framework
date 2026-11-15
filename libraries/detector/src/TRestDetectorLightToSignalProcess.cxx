/// \class TRestDetectorLightToSignalProcess
/// Converts the photoelectrons of a TRestDetectorHitsEvent (energies in photoelectrons, one hit per input
/// hit and channel, from TRestDetectorLightCollectionProcess) into a TRestDetectorSignalEvent: each
/// photoelectron adds its charge to the signal of its channel (TRestDetectorReadout::GetChannelFromPosition
/// at the position of the hit): 1 in units of the mean single photoelectron charge, or Gauss(1, gainSpread)
/// (truncated at 0) with the relative spread of the gain of the photosensor (e.g. 0.3-0.4 for a PMT, 0.1 for
/// a SiPM). It is added at the time
///
///     t = t_hit + t_scintillation + t_wls + t_jitter
///
/// - t_hit: time of the hit (deposit time plus the propagation to the photosensor of the light collection);
/// - t_scintillation: emission time of the photon, sampled from the scintillation profile of the particle
///   of the hit (TRestDetectorLightResponse::GetScintillationTime with TRestDetectorHitsEvent::GetParticleID;
///   the default profile for unknown particles, 0 without profiles);
/// - t_wls: re-emission delay of a wavelength shifter between the scintillator and the photosensor (e.g. WLS
///   fibres), Exp(decayTime) of the `wls` response of the medium wlsMedium (none if it is not given);
/// - t_jitter: Gaussian time jitter of the photosensor (timeJitter, transit time spread).
///
/// The pulse of a photoelectron and the mean gain of the photosensor go into
/// TRestDetectorSignalShapingProcess (e.g. `exponential` or `gaus` with a shapingTime of a few ns for a PMT):
/// the pulse of an event is then the scintillation profile convolved with the single photoelectron response.
/// The electronics go into TRestDetectorSignalToRawSignalProcess. With sampling > 0 the times are rounded
/// down to multiples of it (fewer points per signal); otherwise the times are exact. Hits without channel
/// are counted in NHitsWithoutChannel, NPhotoelectrons is the number of photoelectrons of the signals.
/// Events without signals are rejected.
///
/// Parameters: readout (empty: the data readout of the TRestDetectorReadoutManager), medium (name of the
/// TRestDetectorMedium with `light`; empty: the only one), wlsMedium (name of the TRestDetectorMedium with
/// `wls`; empty: no wavelength shifter), timeJitter (0), gainSpread (0), sampling (0), seed (0: random seed;
/// otherwise the generator is re-seeded at each event from the seed and the event and subevent IDs).
///
/// Example (a scintillator with pulse shape discrimination and SiPMs):
///
///    medium:
///      class: "TRestDetectorMedium"
///      name: "LiquidScintillator"
///      state: "liquid"
///      light:
///        lightYield: 10
///        scintillation:
///          - components:
///              - {fraction: 0.8, riseTime: 0.9 ns, decayTime: 2.1 ns}
///              - {fraction: 0.2, decayTime: 20 ns}
///          - particles: ["alpha", "proton", "ion"]
///            components:
///              - {fraction: 0.5, decayTime: 2.1 ns}
///              - {fraction: 0.5, decayTime: 60 ns}
///    ...
///    lightToSignal:
///      class: "TRestDetectorLightToSignalProcess"
///      medium: "LiquidScintillator"
///      timeJitter: 0.2 ns
///      gainSpread: 0.1
///      observables: ["NPhotoelectrons", "NHitsWithoutChannel"]
///
/// New in restv2 (legacy had no scintillation time profiles).

#include "TRestDetectorLightToSignalProcess.h"

#include <algorithm>
#include <cmath>
#include <map>
#include <stdexcept>
#include <vector>

#include "TRestDetectorHitsEvent.h"
#include "TRestDetectorMedium.h"
#include "TRestDetectorReadout.h"
#include "TRestDetectorReadoutManager.h"
#include "TRestDetectorSignalEvent.h"
#include "TRestManager.h"

static const bool TRestDetectorLightToSignalProcess_FieldsRegistered = []() {
    auto& reg = TRestMetadataFieldRegistry::Instance();
    using P = TRestDetectorLightToSignalProcess;
    reg.RegisterField<P>("readout", &P::fReadoutName);
    reg.RegisterField<P>("medium", &P::fMediumName);
    reg.RegisterField<P>("wlsMedium", &P::fWLSMediumName);
    reg.RegisterField<P>("timeJitter", &P::fTimeJitter);
    reg.RegisterField<P>("gainSpread", &P::fGainSpread);
    reg.RegisterField<P>("sampling", &P::fSampling);
    reg.RegisterField<P>("seed", &P::fSeed);
    return true;
}();

namespace {
const bool kRegistered = []() {
    MetadataClassRegistry::Instance().Register(
        "TRestDetectorLightToSignalProcess", [](const std::string& instanceName, const YAML::Node& params) {
            return std::make_unique<TRestDetectorLightToSignalProcess>(instanceName, params);
        });
    return true;
}();
}  // namespace

TRestDetectorLightToSignalProcess::TRestDetectorLightToSignalProcess() : TRestEventProcess() {
    fName = "TRestDetectorLightToSignalProcess";
}

TRestDetectorLightToSignalProcess::TRestDetectorLightToSignalProcess(const std::string& instanceName,
                                                                     const YAML::Node& node)
    : TRestEventProcess(instanceName, node) {
    LoadConfig();
}

TRestDetectorLightToSignalProcess::TRestDetectorLightToSignalProcess(const std::string& fileName,
                                                                     const std::string& sectionName)
    : TRestEventProcess(fileName, sectionName) {
    LoadConfig();
}

void TRestDetectorLightToSignalProcess::LoadConfig() {
    TRestEventProcess::LoadConfig();
    if (!fNode || fNode.IsNull()) {
        RESTError << "TRestDetectorLightToSignalProcess::LoadConfig - YAML node is missing" << RESTendl;
        return;
    }
    UpdateParamsFromYAML<TRestDetectorLightToSignalProcess>(fNode);
    UpdateYAMLFromParams<TRestDetectorLightToSignalProcess>(fNode);
}

void TRestDetectorLightToSignalProcess::InitProcess() {
    const std::string user = "TRestDetectorLightToSignalProcess";
    if (fTimeJitter < 0 || fGainSpread < 0 || fSampling < 0) {
        throw std::runtime_error(user + ": timeJitter, gainSpread and sampling must be >= 0");
    }
    if (!fManager) throw std::runtime_error(user + ": no TRestManager");

    auto* readouts =
        dynamic_cast<TRestDetectorReadoutManager*>(fManager->GetMetadataClass("TRestDetectorReadoutManager"));
    if (!readouts) throw std::runtime_error(user + ": no TRestDetectorReadoutManager found");
    fReadout = readouts->SelectReadout(fReadoutName, user);

    const auto* medium = fManager->FindMetadata<TRestDetectorMedium>(fMediumName, user, "medium");
    if (!medium) throw std::runtime_error(user + ": needs a TRestDetectorMedium");
    fLight = &medium->GetLightResponse(user);
    fMediumName = medium->GetName();
    fWLS = nullptr;
    if (!fWLSMediumName.empty()) {
        const auto* wlsMedium =
            fManager->FindMetadata<TRestDetectorMedium>(fWLSMediumName, user, "wlsMedium");
        fWLS = &wlsMedium->GetWLSResponse(user);
    }
    // The stored metadata keeps the readout and medium used
    UpdateYAMLFromParams<TRestDetectorLightToSignalProcess>(fNode);

    // Seed 0: random seed. Otherwise re-seeded at each event in ProcessEvent
    fRandom.SetSeed(fSeed);

    RegisterObservable("NPhotoelectrons", fNPhotoelectrons);
    RegisterObservable("NHitsWithoutChannel", fNHitsWithoutChannel);
}

bool TRestDetectorLightToSignalProcess::ProcessEvent(const TRestEvent& input, TRestEvent& output) {
    const auto* hitsEvent = dynamic_cast<const TRestDetectorHitsEvent*>(&input);
    auto* signalEvent = dynamic_cast<TRestDetectorSignalEvent*>(&output);
    if (!hitsEvent || !signalEvent) {
        throw std::runtime_error(
            "TRestDetectorLightToSignalProcess: input/output is not a TRestDetectorHitsEvent/"
            "TRestDetectorSignalEvent");
    }

    signalEvent->Initialize();
    signalEvent->SetEventInfo(hitsEvent);
    if (fSeed != 0) fRandom.SetSeed(TRestDetectorHitsEvent::EventSeed(fSeed, *hitsEvent));
    fNPhotoelectrons = 0;
    fNHitsWithoutChannel = 0;

    // Photoelectrons of each channel versus time, sorted
    std::map<int, std::map<double, double>> charges;
    for (size_t n = 0; n < hitsEvent->GetNumberOfHits(); n++) {
        const int photoelectrons = static_cast<int>(std::lround(hitsEvent->GetEnergy(n)));
        if (photoelectrons <= 0) continue;
        const auto position = hitsEvent->GetPosition(n);
        const int daqID = fReadout->GetChannelFromPosition(position.X(), position.Y(), position.Z());
        if (daqID < 0) {
            fNHitsWithoutChannel++;
            continue;
        }
        const int particleID = hitsEvent->GetParticleID(n);
        auto& signal = charges[daqID];
        for (int p = 0; p < photoelectrons; p++) {
            double time = hitsEvent->GetTime(n) + fLight->GetScintillationTime(particleID, fRandom);
            if (fWLS) time += fWLS->GetReemissionTime(fRandom);
            if (fTimeJitter > 0) time += fRandom.Gaus(0, fTimeJitter);
            if (fSampling > 0) time = std::floor(time / fSampling) * fSampling;
            signal[time] += fGainSpread > 0 ? std::max(0.0, fRandom.Gaus(1, fGainSpread)) : 1.0;
        }
        fNPhotoelectrons += photoelectrons;
    }

    for (const auto& [daqID, points] : charges) {
        signalEvent->AddSignal(daqID, std::vector<std::pair<double, double>>(points.begin(), points.end()));
    }

    RESTDebug << "TRestDetectorLightToSignalProcess: " << signalEvent->GetNumberOfSignals() << " signals, "
              << fNPhotoelectrons << " photoelectrons" << RESTendl;
    return signalEvent->GetNumberOfSignals() > 0;
}
