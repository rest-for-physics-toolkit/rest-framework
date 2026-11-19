/// \class TRestDetectorPhotosensorProcess
/// Response of the photosensors of the readout channels, SiPMs or PMTs, to the photoelectrons of a
/// TRestDetectorHitsEvent: TRestDetectorLightCollectionProcess gives a hit per photoelectron (charge 1, the
/// time its photon reaches the photosensor), this process gives a hit per avalanche of the photosensor, with
/// the same channel (position), time and particle, adding the effects of the photosensor:
///
/// - dark counts: Poisson(darkCountRate * T) avalanches per channel, uniform in darkCountWindow (of length
///   T), on every channel of the readout (particle 0), at a position of the channel (its read end, or its
///   centre if the read end is not in the channel, e.g. the PMT of a veto);
/// - optical crosstalk (SiPM): each avalanche triggers another one at the same time with probability
///   crosstalkProbability, which may trigger another one (geometric number of crosstalk avalanches);
/// - afterpulses: each avalanche is followed with probability afterpulseProbability by another one, delayed
///   by Exp(afterpulseDelay) (SiPM, trapped carriers) or by Gauss(afterpulseDelay, afterpulseSpread) (PMT,
///   ions); afterpulses and crosstalk avalanches may give crosstalk and afterpulses in turn;
/// - saturation (SiPM): each avalanche fires a random cell of the `cells` of its channel; an avalanche on a
///   cell that fired less than cellRecoveryTime before (0: at any time of the event) is lost. For n
///   simultaneous photoelectrons the number of avalanches is cells * (1 - exp(-n / cells)) on average;
/// - the time jitter of the photosensor, Gaussian of timeJitter (transit time spread), and the charge of an
///   avalanche, 1 or Gauss(1, gainSpread) truncated at 0 (gain fluctuation, e.g. 0.3-0.4 for a PMT, 0.1 for
///   a SiPM), in units of the mean single photoelectron charge.
///
/// The avalanches of a channel are processed in time order, from the input photoelectrons (an input hit of
/// charge q is round(q) photoelectrons) and the dark counts. Hits without a channel
/// (TRestDetectorReadout::GetChannelFromPosition) are dropped. The effects not wanted are left at 0 (all 0:
/// one hit of charge 1 per photoelectron, unchanged). The quantum efficiency is that of
/// TRestDetectorLightCollectionProcess (it decides which photons give photoelectrons); the pulse of an
/// avalanche and the mean gain go into TRestDetectorSignalShapingProcess. The saturation of the charge of a
/// PMT (non-linearity at high currents) is not modelled.
///
/// Observables: NAvalanches (output hits), NDarkCounts, NCrosstalk, NAfterpulses (avalanches of each kind
/// in the output) and NSaturationLosses (avalanches lost on busy cells).
///
/// Parameters: readout (empty: the data readout of the TRestDetectorReadoutManager), timeJitter (0),
/// gainSpread (0), darkCountRate (0, per channel, e.g. "100 kHz"), darkCountWindow (needed with
/// darkCountRate), crosstalkProbability (0), afterpulseProbability (0), afterpulseDelay, afterpulseSpread
/// (0: exponential delay), cells (0: no saturation), cellRecoveryTime (0), seed (0: random seed; otherwise
/// the generator is re-seeded at each event from the seed and the event and subevent IDs).
/// crosstalkProbability + afterpulseProbability must be below 1 (a finite number of avalanches).
///
/// Examples:
///
///    # SiPMs (LiquidO fibres)
///    photosensor:
///      class: "TRestDetectorPhotosensorProcess"
///      timeJitter: 0.2 ns
///      gainSpread: 0.1
///      darkCountRate: 100 kHz
///      darkCountWindow: [-4 us, 16 us]
///      crosstalkProbability: 0.05
///      afterpulseProbability: 0.02
///      afterpulseDelay: 50 ns
///      cells: 3600
///      cellRecoveryTime: 50 ns
///      observables: ["NAvalanches", "NDarkCounts", "NCrosstalk", "NAfterpulses", "NSaturationLosses"]
///
///    # PMT (vetoes): ion afterpulses about 1.5 us after the pulse
///    photosensor:
///      class: "TRestDetectorPhotosensorProcess"
///      readout: "Veto_Readout"
///      timeJitter: 0.5 ns
///      gainSpread: 0.35
///      darkCountRate: 1 kHz
///      darkCountWindow: [-4 us, 16 us]
///      afterpulseProbability: 0.01
///      afterpulseDelay: 1.5 us
///      afterpulseSpread: 0.2 us
///
/// The time jitter and the gain spread were parameters of TRestDetectorLightCollectionProcess in the
/// previous version of restv2 (timeJitter, gainSpread).

#include "TRestDetectorPhotosensorProcess.h"

#include <algorithm>
#include <cmath>
#include <memory>
#include <queue>
#include <stdexcept>
#include <unordered_map>
#include <vector>

#include "TRestDetectorHitsEvent.h"
#include "TRestDetectorReadout.h"
#include "TRestDetectorReadoutManager.h"
#include "TRestManager.h"

using ROOT::Math::XYZVector;

static const bool TRestDetectorPhotosensorProcess_FieldsRegistered = []() {
    auto& reg = TRestMetadataFieldRegistry::Instance();
    using P = TRestDetectorPhotosensorProcess;
    reg.RegisterField<P>("readout", &P::fReadoutName);
    reg.RegisterField<P>("timeJitter", &P::fTimeJitter);
    reg.RegisterField<P>("gainSpread", &P::fGainSpread);
    reg.RegisterField<P>("darkCountRate", &P::fDarkCountRate);
    reg.RegisterField<P>("darkCountWindow", &P::fDarkCountWindow);
    reg.RegisterField<P>("crosstalkProbability", &P::fCrosstalkProbability);
    reg.RegisterField<P>("afterpulseProbability", &P::fAfterpulseProbability);
    reg.RegisterField<P>("afterpulseDelay", &P::fAfterpulseDelay);
    reg.RegisterField<P>("afterpulseSpread", &P::fAfterpulseSpread);
    reg.RegisterField<P>("cells", &P::fCells);
    reg.RegisterField<P>("cellRecoveryTime", &P::fCellRecoveryTime);
    reg.RegisterField<P>("seed", &P::fSeed);
    return true;
}();

namespace {
const bool kRegistered = []() {
    MetadataClassRegistry::Instance().Register(
        "TRestDetectorPhotosensorProcess", [](const std::string& instanceName, const YAML::Node& params) {
            return std::make_unique<TRestDetectorPhotosensorProcess>(instanceName, params);
        });
    return true;
}();

/// Kind of an avalanche (its observable)
enum class Kind { Photoelectron, DarkCount, Crosstalk, Afterpulse };

/// An avalanche of a channel, before the saturation, the jitter and the gain
struct Avalanche {
    double time;
    XYZVector position;
    int particleID;
    Kind kind;
    /// Earliest first (std::priority_queue keeps the largest on top)
    bool operator<(const Avalanche& other) const { return time > other.time; }
};
}  // namespace

TRestDetectorPhotosensorProcess::TRestDetectorPhotosensorProcess() : TRestEventProcess() {
    fName = "TRestDetectorPhotosensorProcess";
}

TRestDetectorPhotosensorProcess::TRestDetectorPhotosensorProcess(const std::string& instanceName,
                                                                 const YAML::Node& node)
    : TRestEventProcess(instanceName, node) {
    LoadConfig();
}

TRestDetectorPhotosensorProcess::TRestDetectorPhotosensorProcess(const std::string& fileName,
                                                                 const std::string& sectionName)
    : TRestEventProcess(fileName, sectionName) {
    LoadConfig();
}

void TRestDetectorPhotosensorProcess::LoadConfig() {
    TRestEventProcess::LoadConfig();
    if (!fNode || fNode.IsNull()) {
        RESTError << "TRestDetectorPhotosensorProcess::LoadConfig - YAML node is missing" << RESTendl;
        return;
    }
    UpdateParamsFromYAML<TRestDetectorPhotosensorProcess>(fNode);
    const std::string prefix = "TRestDetectorPhotosensorProcess: ";
    if (fTimeJitter < 0 || fGainSpread < 0 || fAfterpulseSpread < 0 || fCellRecoveryTime < 0) {
        throw std::runtime_error(
            prefix + "timeJitter, gainSpread, afterpulseSpread and cellRecoveryTime must be >= 0");
    }
    if (fDarkCountRate < 0) throw std::runtime_error(prefix + "darkCountRate must be >= 0");
    if (fDarkCountRate > 0 && !(fDarkCountWindow.second > fDarkCountWindow.first)) {
        throw std::runtime_error(prefix +
                                 "darkCountRate needs a darkCountWindow [start, end] with end > start");
    }
    for (double probability : {fCrosstalkProbability, fAfterpulseProbability}) {
        if (probability < 0 || probability >= 1) {
            throw std::runtime_error(prefix +
                                     "crosstalkProbability and afterpulseProbability must be in [0, 1)");
        }
    }
    // Mean number of avalanches triggered by one: below 1, or the chain does not end
    if (fCrosstalkProbability + fAfterpulseProbability >= 1) {
        throw std::runtime_error(prefix + "crosstalkProbability + afterpulseProbability must be below 1");
    }
    if (fAfterpulseProbability > 0 && fAfterpulseDelay <= 0) {
        throw std::runtime_error(prefix + "afterpulseProbability needs afterpulseDelay > 0");
    }
    if (fCells < 0) throw std::runtime_error(prefix + "cells must be >= 0");
    UpdateYAMLFromParams<TRestDetectorPhotosensorProcess>(fNode);
}

void TRestDetectorPhotosensorProcess::InitProcess() {
    const std::string user = "TRestDetectorPhotosensorProcess";
    if (!fManager) throw std::runtime_error(user + ": no TRestManager");
    auto* readouts =
        dynamic_cast<TRestDetectorReadoutManager*>(fManager->GetMetadataClass("TRestDetectorReadoutManager"));
    if (!readouts) throw std::runtime_error(user + ": no TRestDetectorReadoutManager found");
    fReadout = readouts->SelectReadout(fReadoutName, user);
    // The stored metadata keeps the readout used
    UpdateYAMLFromParams<TRestDetectorPhotosensorProcess>(fNode);

    // A position of each channel for its dark counts: one that the readout gives back as the channel (as
    // TRestDetectorHitsToSignalProcess finds it), the read end or else the centre
    fChannelPositions.clear();
    if (fDarkCountRate > 0) {
        const auto inChannel = [this](int daqID, const XYZVector& p) {
            return std::isfinite(p.X()) && std::isfinite(p.Y()) && std::isfinite(p.Z()) &&
                   fReadout->GetChannelFromPosition(p.X(), p.Y(), p.Z()) == daqID;
        };
        for (int daqID : fReadout->GetDAQChannels()) {
            const XYZVector end = fReadout->GetReadoutPosition(daqID);
            const XYZVector centre = fReadout->GetPositionFromChannel(daqID);
            if (inChannel(daqID, end)) {
                fChannelPositions[daqID] = end;
            } else if (inChannel(daqID, centre)) {
                fChannelPositions[daqID] = centre;
            } else {
                throw std::runtime_error(user + ": no position of the channel " + std::to_string(daqID) +
                                         " of readout " + fReadout->GetName() + " for its dark counts");
            }
        }
        if (fChannelPositions.empty()) {
            throw std::runtime_error(user + ": the readout " + fReadout->GetName() +
                                     " has no channels for the dark counts");
        }
    }

    // Seed 0: random seed. Otherwise re-seeded at each event in ProcessEvent
    fRandom.SetSeed(fSeed);

    RegisterObservable("NAvalanches", fNAvalanches);
    RegisterObservable("NDarkCounts", fNDarkCounts);
    RegisterObservable("NCrosstalk", fNCrosstalk);
    RegisterObservable("NAfterpulses", fNAfterpulses);
    RegisterObservable("NSaturationLosses", fNSaturationLosses);
}

double TRestDetectorPhotosensorProcess::AfterpulseDelay() {
    if (fAfterpulseSpread > 0) return std::max(0.0, fRandom.Gaus(fAfterpulseDelay, fAfterpulseSpread));
    return fRandom.Exp(fAfterpulseDelay);
}

bool TRestDetectorPhotosensorProcess::ProcessEvent(const TRestEvent& input, TRestEvent& output) {
    const auto* inputEvent = dynamic_cast<const TRestDetectorHitsEvent*>(&input);
    auto* outputEvent = dynamic_cast<TRestDetectorHitsEvent*>(&output);
    if (!inputEvent || !outputEvent) {
        throw std::runtime_error(
            "TRestDetectorPhotosensorProcess: input/output is not a TRestDetectorHitsEvent");
    }

    outputEvent->Initialize();
    outputEvent->SetEventInfo(inputEvent);
    if (fSeed != 0) fRandom.SetSeed(TRestDetectorHitsEvent::EventSeed(fSeed, *inputEvent));
    fNAvalanches = 0;
    fNDarkCounts = 0;
    fNCrosstalk = 0;
    fNAfterpulses = 0;
    fNSaturationLosses = 0;

    // The avalanches of each channel, earliest first: the photoelectrons and the dark counts
    std::map<int, std::priority_queue<Avalanche>> channels;
    int lastChannel = -1;
    XYZVector lastPosition;
    for (size_t n = 0; n < inputEvent->GetNumberOfHits(); n++) {
        const XYZVector position = inputEvent->GetPosition(n);
        // The photoelectrons of a channel share its position (TRestDetectorLightCollectionProcess)
        if (n == 0 || position != lastPosition) {
            lastChannel = fReadout->GetChannelFromPosition(position.X(), position.Y(), position.Z());
            lastPosition = position;
        }
        if (lastChannel < 0) continue;
        const long photoelectrons = std::lround(inputEvent->GetEnergy(n));
        auto& queue = channels[lastChannel];
        for (long p = 0; p < photoelectrons; p++) {
            queue.push({inputEvent->GetTime(n), position, inputEvent->GetParticleID(n), Kind::Photoelectron});
        }
    }
    if (fDarkCountRate > 0) {
        const double start = fDarkCountWindow.first, length = fDarkCountWindow.second - start;
        for (const auto& [daqID, position] : fChannelPositions) {
            const int darkCounts = fRandom.Poisson(fDarkCountRate * length);
            if (darkCounts == 0) continue;
            auto& queue = channels[daqID];
            for (int d = 0; d < darkCounts; d++) {
                queue.push({start + length * fRandom.Rndm(), position, 0, Kind::DarkCount});
            }
        }
    }

    // Each channel in time order: saturation, then the crosstalk and afterpulses of each avalanche
    std::unordered_map<int, double> firedCells;  // cell: time of its last avalanche
    for (auto& [daqID, queue] : channels) {
        firedCells.clear();
        while (!queue.empty()) {
            const Avalanche avalanche = queue.top();
            queue.pop();
            if (fCells > 0) {
                const int cell = static_cast<int>(fRandom.Integer(fCells));
                const auto fired = firedCells.find(cell);
                if (fired != firedCells.end() &&
                    (fCellRecoveryTime <= 0 || avalanche.time - fired->second < fCellRecoveryTime)) {
                    fNSaturationLosses++;
                    continue;
                }
                firedCells[cell] = avalanche.time;
            }

            double time = avalanche.time;
            if (fTimeJitter > 0) time += fRandom.Gaus(0, fTimeJitter);
            const double charge = fGainSpread > 0 ? std::max(0.0, fRandom.Gaus(1, fGainSpread)) : 1.0;
            outputEvent->AddHit(avalanche.position, charge, time, TRestDetectorHitsEvent::REST_HitType::XYZ,
                                avalanche.particleID);
            fNAvalanches++;
            if (avalanche.kind == Kind::DarkCount) fNDarkCounts++;
            if (avalanche.kind == Kind::Crosstalk) fNCrosstalk++;
            if (avalanche.kind == Kind::Afterpulse) fNAfterpulses++;

            if (fCrosstalkProbability > 0 && fRandom.Rndm() < fCrosstalkProbability) {
                queue.push({avalanche.time, avalanche.position, avalanche.particleID, Kind::Crosstalk});
            }
            if (fAfterpulseProbability > 0 && fRandom.Rndm() < fAfterpulseProbability) {
                queue.push({avalanche.time + AfterpulseDelay(), avalanche.position, avalanche.particleID,
                            Kind::Afterpulse});
            }
        }
    }

    RESTDebug << "TRestDetectorPhotosensorProcess: event " << inputEvent->GetID() << ", " << fNAvalanches
              << " avalanches (" << fNDarkCounts << " dark counts, " << fNCrosstalk << " crosstalk, "
              << fNAfterpulses << " afterpulses), " << fNSaturationLosses << " lost by saturation"
              << RESTendl;
    return outputEvent->GetNumberOfHits() > 0;
}
