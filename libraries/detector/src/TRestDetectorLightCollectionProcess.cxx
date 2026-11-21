/// \class TRestDetectorLightCollectionProcess
/// Collection of the scintillation photons of the hits of a TRestDetectorHitsEvent (energies in photons,
/// from TRestDetectorScintillationProcess) by the readout channels and their detection by the photosensor of
/// each channel (SiPM, PMT): the photons reaching it give Binomial(N, quantumEfficiency) photoelectrons.
/// Each photoelectron is an output hit, at a position that the readout assigns to its channel (so
/// TRestDetectorHitsToSignalProcess puts it in that channel), with the particle of the input hit, the
/// charge 1 (one photoelectron) and the time its photon reaches the photosensor
///
///     t = t_hit + t_transport + t_scintillation + t_wls
///
/// - t_transport: the delay to the photosensor of the model (below);
/// - t_scintillation: emission time of the photon, sampled from the scintillation profile of the particle
///   of the hit (TRestDetectorLightResponse::GetScintillationTime of the medium with `light`; the default
///   profile for unknown particles, 0 without profiles or without such a medium);
/// - t_wls: re-emission delay of the fibres, Exp(decayTime) of the `wls` response of wlsMedium.
///
/// TRestDetectorPhotosensorProcess then gives the response of the photosensor (time jitter, gain spread,
/// dark counts, crosstalk, afterpulses, saturation) and TRestDetectorHitsToSignalProcess the signals in
/// photoelectrons (its sampling bins the times);
/// the pulse of a photoelectron and the mean gain of the photosensor go into
/// TRestDetectorSignalShapingProcess (e.g. `exponential` or `gaus` with a shapingTime of a few ns for a PMT),
/// so the pulse of an event is the scintillation profile convolved with the single photoelectron response,
/// and the electronics into TRestDetectorSignalToRawSignalProcess. The way the light reaches the channels
/// depends on the scintillator (`model`):
///
/// - `localCapture`, opaque scintillator (e.g. LiquidO: the light stays near where it is produced and the
///   fibres around it capture it): the channels closer than maxDistance to the hit
///   (TRestDetectorReadout::GetChannelsInRange, e.g. the fibres of a TRestFibreReadout, in the plane
///   transverse to them) capture each photon with probability captureProbability(r), r the distance to the
///   channel; the photons of the hit are distributed among them with a multinomial (a photon is captured by
///   one channel at most; the probabilities of the channels of a position must add up to 1 at most).
///   captureProbability includes the transport in the medium and the trapping of the channel towards its
///   photosensor (for a fibre read at both ends, each end is a channel with its own probability): a ROOT
///   formula of x = r in mm (e.g. a Gaussian light ball) or a graph or histogram of a ROOT file
///   ("file.root:name", e.g. from an optical simulation of Geant4). With wlsMedium (the medium of the fibres
///   with a `wls` response), captureProbability is the capture by the fibre and the photons are guided half
///   in each direction along it (TRestDetectorReadout::GetGuidedLightPaths of the channel): a fibre read at
///   both ends gives half to each channel; read at one end, half go directly to it, a distance d, and half
///   to the other end, where a fraction endReflectivity of them is reflected back, a distance 2 L - d (L the
///   length of the fibre). Of the photons on a path of length d, Binomial(N, R exp(-d / attenuationLength))
///   arrive (R = 1 for the direct path), at t + d / lightSpeed. Without wlsMedium the captured photons
///   arrive at once (the transport in the fibre is then part of captureProbability, per channel). The output
///   hits are at the position given by the readout (for a fibre, its read end); their time is that of the
///   hit plus the delay in the fibre.
/// - `attenuation`, transparent scintillator read at one end (e.g. a veto with a light guide and a PMT):
///   all the light goes to the channel of the hit (TRestDetectorReadout::GetChannelFromPosition), the
///   photons reaching its readout are Binomial(N, exp(-d / attenuationLength)) and they arrive at
///   t + d / effectiveLightSpeed, with d the distance to the readout of the channel
///   (TRestDetectorReadout::GetDistanceToReadout, e.g. to the readoutVolume of a TRestGDMLReadout) and the
///   `light` response of the medium (attenuationLength 0: no attenuation; effectiveLightSpeed 0: no delay). A
///   channel without distance to its readout is an error. The output hits are at the position of the hit.
///
/// Hits without any channel collecting their light are counted in the observable NHitsWithoutChannel
/// (their photons are lost); NPhotoelectrons is the total of the event. Events without photoelectrons are
/// rejected.
///
/// Parameters: model (localCapture), readout (empty: the data readout of the TRestDetectorReadoutManager),
/// medium (name of the TRestDetectorMedium with `light`; empty: the only one with `light`, needed for
/// attenuation), captureProbability, maxDistance and wlsMedium (localCapture), quantumEfficiency (1), seed
/// (0: random seed; otherwise the generator is re-seeded at each event from the seed and the event and
/// subevent IDs).
///
/// Differences with the previous version of restv2: the times and charges of the photoelectrons were given
/// by TRestDetectorLightToSignalProcess (removed), now TRestDetectorHitsToSignalProcess as for the charge;
/// with wlsMedium captureProbability was per channel and all the captured light went to its read end, now it
/// is per fibre, half in each direction (the reflected half needs endReflectivity); timeJitter and
/// gainSpread were parameters of this process, now of TRestDetectorPhotosensorProcess.
///
/// Examples:
///
///    # LiquidO (TRestFibreReadout): Gaussian light ball of 10 mm, SiPMs
///    lightCollection:
///      class: "TRestDetectorLightCollectionProcess"
///      model: "localCapture"
///      captureProbability: "0.03*exp(-x*x/(2*10*10))"
///      maxDistance: 30 mm
///      wlsMedium: "Y11"         # medium of the fibres with a `wls` response
///      quantumEfficiency: 0.4   # SiPMs
///      observables: ["NPhotoelectrons", "NHitsWithoutChannel"]
///    photosensor:
///      class: "TRestDetectorPhotosensorProcess"
///      timeJitter: 0.2 ns
///      gainSpread: 0.1
///    hitsToSignal:
///      class: "TRestDetectorHitsToSignalProcess"
///      sampling: 0.04 us
///
///    # Vetoes (TRestGDMLReadout with readoutVolume), medium with attenuationLength and effectiveLightSpeed
///    lightCollection:
///      class: "TRestDetectorLightCollectionProcess"
///      model: "attenuation"
///      readout: "Veto_Readout"
///      medium: "BC408"
///      quantumEfficiency: 0.25  # PMT
///
/// The attenuation model replaces the legacy TRestDetectorLightAttenuationProcess (energy E exp(-d/λ) of the
/// veto hits), now with photons; its parameters (vetoLightAttenuationLength, vetoEffectiveLightSpeed) are
/// the light response of the medium and the distance is to the readout volume of the channel.

#include "TRestDetectorLightCollectionProcess.h"

#include <algorithm>
#include <cmath>
#include <memory>
#include <stdexcept>

#include "TF1.h"
#include "TFile.h"
#include "TGraph.h"
#include "TH1.h"
#include "TRestDetectorHitsEvent.h"
#include "TRestDetectorMedium.h"
#include "TRestDetectorReadout.h"
#include "TRestDetectorReadoutManager.h"
#include "TRestManager.h"
#include "TRestTools.h"

using ROOT::Math::XYZVector;

namespace {
/// Binomial(n, p) drawn from a TRandom in about constant time: TRandom::Binomial draws n uniform numbers,
/// too slow for the photons of a hit (thousands, for each channel around it). Not
/// std::binomial_distribution, whose algorithm depends on the standard library (the same seed would give
/// other numbers on another platform). For a mean below 10, inversion of the distribution function; above,
/// rejection from a Lorentzian envelope, as TRandom::Poisson.
int Binomial(TRandom& random, int n, double p) {
    if (n <= 0 || !(p > 0)) return 0;
    if (p >= 1) return n;
    if (p > 0.5) return n - Binomial(random, n, 1 - p);
    const double q = 1 - p;
    const double mean = n * p;
    if (mean < 10) {
        // The first k whose cumulated probability exceeds the uniform (P(0) = q^n >= exp(-20))
        double u = random.Rndm();
        double probability = std::pow(q, n);
        int k = 0;
        while (u > probability && k < n) {
            u -= probability;
            probability *= (n - k) * p / ((k + 1) * q);
            k++;
        }
        return k;
    }
    const double pi = std::acos(-1.0);
    const double logFactorialN = std::lgamma(n + 1.0);
    const double logP = std::log(p);
    const double logQ = std::log(q);
    const double width = std::sqrt(2 * mean * q);
    while (true) {
        double y, k;
        do {
            y = std::tan(pi * random.Rndm());
            k = width * y + mean;
        } while (k < 0 || k >= n + 1);
        k = std::floor(k);
        // Probability of k over the envelope (at most 1)
        const double ratio =
            1.2 * width * (1 + y * y) *
            std::exp(logFactorialN - std::lgamma(k + 1) - std::lgamma(n - k + 1) + k * logP + (n - k) * logQ);
        if (random.Rndm() <= ratio) return static_cast<int>(k);
    }
}
}  // namespace

static const bool TRestDetectorLightCollectionProcess_FieldsRegistered = []() {
    auto& reg = TRestMetadataFieldRegistry::Instance();
    using P = TRestDetectorLightCollectionProcess;
    reg.RegisterField<P>("model", &P::fModel);
    reg.RegisterField<P>("readout", &P::fReadoutName);
    reg.RegisterField<P>("medium", &P::fMediumName);
    reg.RegisterField<P>("wlsMedium", &P::fWLSMediumName);
    reg.RegisterField<P>("captureProbability", &P::fCaptureProbability);
    reg.RegisterField<P>("maxDistance", &P::fMaxDistance);
    reg.RegisterField<P>("quantumEfficiency", &P::fQuantumEfficiency);
    reg.RegisterField<P>("seed", &P::fSeed);
    return true;
}();

namespace {
const bool kRegistered = []() {
    MetadataClassRegistry::Instance().Register(
        "TRestDetectorLightCollectionProcess", [](const std::string& instanceName, const YAML::Node& params) {
            return std::make_unique<TRestDetectorLightCollectionProcess>(instanceName, params);
        });
    return true;
}();

/// Function of the distance (mm) from "file.root:name" (graph or histogram) or a formula of x
std::function<double(double)> MakeCaptureFunction(const std::string& text, double maxDistance,
                                                  const std::string& prefix) {
    const auto separator = text.rfind(':');
    const auto extension = text.find(".root");
    if (separator != std::string::npos && extension != std::string::npos && separator > extension) {
        const std::string fileName = text.substr(0, separator);
        const std::string objectName = text.substr(separator + 1);
        std::unique_ptr<TFile> file(TFile::Open(fileName.c_str()));
        if (!file || file->IsZombie()) throw std::runtime_error(prefix + "cannot open " + fileName);
        TObject* object = file->Get(objectName.c_str());
        if (auto* graph = dynamic_cast<TGraph*>(object)) {
            std::shared_ptr<TGraph> copy(static_cast<TGraph*>(graph->Clone()));
            return [copy](double r) { return copy->Eval(r); };
        }
        if (auto* histogram = dynamic_cast<TH1*>(object)) {
            std::shared_ptr<TH1> copy(static_cast<TH1*>(histogram->Clone()));
            copy->SetDirectory(nullptr);
            return [copy](double r) { return copy->Interpolate(r); };
        }
        throw std::runtime_error(prefix + "no graph or histogram '" + objectName + "' in " + fileName);
    }
    // Not added to the list of functions of ROOT (several processes may have the same formula)
    auto function =
        std::make_shared<TF1>("lightCaptureProbability", text.c_str(), 0, maxDistance, TF1::EAddToList::kNo);
    if (!function->IsValid()) {
        throw std::runtime_error(prefix + "captureProbability '" + text + "' is not a valid formula of x");
    }
    return [function](double r) { return function->Eval(r); };
}
}  // namespace

TRestDetectorLightCollectionProcess::TRestDetectorLightCollectionProcess() : TRestEventProcess() {
    fName = "TRestDetectorLightCollectionProcess";
}

TRestDetectorLightCollectionProcess::TRestDetectorLightCollectionProcess(const std::string& instanceName,
                                                                         const YAML::Node& node)
    : TRestEventProcess(instanceName, node) {
    LoadConfig();
}

TRestDetectorLightCollectionProcess::TRestDetectorLightCollectionProcess(const std::string& fileName,
                                                                         const std::string& sectionName)
    : TRestEventProcess(fileName, sectionName) {
    LoadConfig();
}

void TRestDetectorLightCollectionProcess::LoadConfig() {
    TRestEventProcess::LoadConfig();
    if (!fNode || fNode.IsNull()) {
        RESTError << "TRestDetectorLightCollectionProcess::LoadConfig - YAML node is missing" << RESTendl;
        return;
    }
    UpdateParamsFromYAML<TRestDetectorLightCollectionProcess>(fNode);
    const std::string prefix = "TRestDetectorLightCollectionProcess: ";
    fModel = TRestTools::NormalizeOption("model", fModel, {"localCapture", "attenuation"});
    if (fQuantumEfficiency < 0 || fQuantumEfficiency > 1) {
        throw std::runtime_error(prefix + "quantumEfficiency must be in [0, 1]");
    }
    fCapture = nullptr;
    if (fModel == "localCapture") {
        if (fMaxDistance <= 0) throw std::runtime_error(prefix + "localCapture needs maxDistance > 0");
        if (fCaptureProbability.empty()) {
            throw std::runtime_error(prefix + "localCapture needs captureProbability");
        }
        fCapture = MakeCaptureFunction(fCaptureProbability, fMaxDistance, prefix);
    }
    UpdateYAMLFromParams<TRestDetectorLightCollectionProcess>(fNode);
}

double TRestDetectorLightCollectionProcess::GetCaptureProbability(double distance) const {
    if (!fCapture || distance > fMaxDistance) return 0;
    return std::max(0.0, fCapture(distance));
}

void TRestDetectorLightCollectionProcess::InitProcess() {
    const std::string user = "TRestDetectorLightCollectionProcess";
    if (!fManager) throw std::runtime_error(user + ": no TRestManager");

    auto* readouts =
        dynamic_cast<TRestDetectorReadoutManager*>(fManager->GetMetadataClass("TRestDetectorReadoutManager"));
    if (!readouts) throw std::runtime_error(user + ": no TRestDetectorReadoutManager found");
    fReadout = readouts->SelectReadout(fReadoutName, user);

    // The light response: scintillation profiles, and the attenuation of the attenuation model
    fLight = nullptr;
    const auto* medium = TRestDetectorMedium::Find(*fManager, fMediumName, "light", user);
    if (medium) {
        fLight = &medium->GetLightResponse(user);
        fMediumName = medium->GetName();
    } else if (fModel == "attenuation") {
        throw std::runtime_error(user + ": the attenuation model needs a TRestDetectorMedium with light");
    }
    fWLS = nullptr;
    if (!fWLSMediumName.empty()) {
        if (fModel != "localCapture") throw std::runtime_error(user + ": wlsMedium is for localCapture");
        const auto* wlsMedium =
            TRestDetectorMedium::Find(*fManager, fWLSMediumName, "wls", user, "wlsMedium");
        fWLS = &wlsMedium->GetWLSResponse(user);
    }
    // The stored metadata keeps the readout and medium used
    UpdateYAMLFromParams<TRestDetectorLightCollectionProcess>(fNode);

    // Seed 0: random seed. Otherwise re-seeded at each event in ProcessEvent
    fRandom.SetSeed(fSeed);

    RegisterObservable("NPhotoelectrons", fNPhotoelectrons);
    RegisterObservable("NHitsWithoutChannel", fNHitsWithoutChannel);
}

void TRestDetectorLightCollectionProcess::AddPhotoelectrons(const XYZVector& position, int photoelectrons,
                                                            double time, int particleID,
                                                            TRestDetectorHitsEvent& output) {
    if (photoelectrons <= 0) return;
    const auto* profile = fLight ? fLight->GetScintillationProfile(particleID) : nullptr;
    for (int p = 0; p < photoelectrons; p++) {
        double t = time;
        if (profile) t += profile->SampleTime(fRandom);
        if (fWLS) t += fWLS->GetReemissionTime(fRandom);
        output.AddHit(position, 1.0, t, TRestDetectorHitsEvent::REST_HitType::XYZ, particleID);
    }
    fNPhotoelectrons += photoelectrons;
}

int TRestDetectorLightCollectionProcess::Detect(int photons) {
    if (photons <= 0) return 0;
    if (fQuantumEfficiency >= 1) return photons;
    return Binomial(fRandom, photons, fQuantumEfficiency);
}

void TRestDetectorLightCollectionProcess::CaptureLocally(const XYZVector& position, int photons, double time,
                                                         int particleID, TRestDetectorHitsEvent& output) {
    const auto channels = fReadout->GetChannelsInRange(position, fMaxDistance);
    // Multinomial: each channel takes a binomial of the photons left, with its probability relative to
    // what is left of the total
    int left = photons;
    double probabilityLeft = 1;
    bool collected = false;
    for (const auto& channel : channels) {
        // With WLS fibres the capture is by the fibre, half of the light in each direction: a channel
        // takes the halves of the paths reaching its photosensor
        std::vector<TRestDetectorReadout::GuidedLightPath> paths;
        double probability = GetCaptureProbability(channel.distance);
        if (fWLS && probability > 0) {
            paths = fReadout->GetGuidedLightPaths(channel.daqID, position);
            if (paths.empty()) {
                throw std::runtime_error("TRestDetectorLightCollectionProcess: the channel " +
                                         std::to_string(channel.daqID) + " of readout " +
                                         fReadout->GetName() +
                                         " has no distance to its readout (needed by wlsMedium)");
            }
            probability *= std::min(1.0, paths.size() / 2.0);
        }
        if (probability <= 0) continue;
        collected = true;
        if (probability > probabilityLeft + 1e-9) {
            throw std::runtime_error(
                "TRestDetectorLightCollectionProcess: the capture probabilities of the channels around a "
                "position add up to more than 1 (captureProbability too large for the distance between "
                "channels)");
        }
        if (left <= 0) {
            probabilityLeft -= probability;
            continue;
        }
        const int captured =
            probability >= probabilityLeft ? left : Binomial(fRandom, left, probability / probabilityLeft);
        probabilityLeft -= probability;
        left -= captured;
        if (!fWLS) {
            AddPhotoelectrons(channel.position, Detect(captured), time, particleID, output);
            continue;
        }
        // Along the fibre to its read end: the captured photons shared equally among the paths
        int share = captured;
        for (size_t k = 0; k < paths.size(); k++) {
            const int onPath =
                k + 1 == paths.size() ? share : Binomial(fRandom, share, 1.0 / (paths.size() - k));
            share -= onPath;
            double survival = paths[k].reflected ? fWLS->GetEndReflectivity() : 1.0;
            const double attenuationLength = fWLS->GetAttenuationLength();
            if (attenuationLength > 0) survival *= std::exp(-paths[k].distance / attenuationLength);
            const double lightSpeed = fWLS->GetLightSpeed();
            const double arrival = lightSpeed > 0 ? time + paths[k].distance / lightSpeed : time;
            AddPhotoelectrons(channel.position, Detect(Binomial(fRandom, onPath, survival)), arrival,
                              particleID, output);
        }
    }
    if (!collected) fNHitsWithoutChannel++;
}

void TRestDetectorLightCollectionProcess::Attenuate(const XYZVector& position, int photons, double time,
                                                    int particleID, TRestDetectorHitsEvent& output) {
    const int channel = fReadout->GetChannelFromPosition(position.X(), position.Y(), position.Z());
    if (channel < 0) {
        fNHitsWithoutChannel++;
        return;
    }
    const double distance = fReadout->GetDistanceToReadout(channel, position);
    if (std::isnan(distance)) {
        throw std::runtime_error("TRestDetectorLightCollectionProcess: the channel " +
                                 std::to_string(channel) + " of readout " + fReadout->GetName() +
                                 " has no distance to its readout (e.g. no readoutVolume)");
    }
    const double attenuationLength = fLight->GetAttenuationLength();
    const double lightSpeed = fLight->GetEffectiveLightSpeed();
    const int arriving =
        attenuationLength > 0 ? Binomial(fRandom, photons, std::exp(-distance / attenuationLength)) : photons;
    AddPhotoelectrons(position, Detect(arriving), lightSpeed > 0 ? time + distance / lightSpeed : time,
                      particleID, output);
}

bool TRestDetectorLightCollectionProcess::ProcessEvent(const TRestEvent& input, TRestEvent& output) {
    const auto* inputEvent = dynamic_cast<const TRestDetectorHitsEvent*>(&input);
    auto* outputEvent = dynamic_cast<TRestDetectorHitsEvent*>(&output);
    if (!inputEvent || !outputEvent) {
        throw std::runtime_error(
            "TRestDetectorLightCollectionProcess: input/output is not a TRestDetectorHitsEvent");
    }

    outputEvent->Initialize();
    outputEvent->SetEventInfo(inputEvent);
    if (fSeed != 0) fRandom.SetSeed(TRestEvent::EventSeed(fSeed, *inputEvent));
    fNPhotoelectrons = 0;
    fNHitsWithoutChannel = 0;

    for (size_t n = 0; n < inputEvent->GetNumberOfHits(); n++) {
        const int photons = static_cast<int>(std::lround(inputEvent->GetEnergy(n)));
        if (photons <= 0) continue;
        if (fModel == "localCapture") {
            CaptureLocally(inputEvent->GetPosition(n), photons, inputEvent->GetTime(n),
                           inputEvent->GetParticleID(n), *outputEvent);
        } else {
            Attenuate(inputEvent->GetPosition(n), photons, inputEvent->GetTime(n),
                      inputEvent->GetParticleID(n), *outputEvent);
        }
    }

    RESTDebug << "TRestDetectorLightCollectionProcess: event " << inputEvent->GetID() << ", "
              << fNPhotoelectrons << " photoelectrons, " << fNHitsWithoutChannel << " hits without channel"
              << RESTendl;
    return outputEvent->GetNumberOfHits() > 0;
}
