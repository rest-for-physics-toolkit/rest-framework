/// \class TRestDetectorHitsToSignalProcess
/// Converts a TRestDetectorHitsEvent into a TRestDetectorSignalEvent: the energy of each hit is added to the
/// signal of the readout channel at its position (TRestDetectorReadout::GetChannelFromPosition), at the
/// time of the hit. It is generic: the hits must already be where and when they are collected, e.g. the
/// electrons drifted to the readout plane at their arrival time (TRestDetectorHitsDriftProcess) or the
/// deposits of a scintillator read as a whole. Hits without channel are rejected and counted in the
/// observable NHitsRejected.
///
/// With sampling > 0 the times are rounded down to multiples of it (fewer points per signal); otherwise the
/// times are exact. Events without signals are rejected.
///
/// Differences with legacy: the drift (projection on the readout, drift time and drift volume cuts) is done
/// by TRestDetectorHitsDriftProcess, there is no light-speed time for veto channels (they go in their own
/// pipeline; the light has its own processes), no readout planes loop and the sampling is optional.
///
/// Parameters: readout (empty: the data readout of the TRestDetectorReadoutManager), sampling (0).
///
/// Example:
///
///    hitsToSignal:
///      class: "TRestDetectorHitsToSignalProcess"
///      observables:
///        - "NHitsRejected"

#include "TRestDetectorHitsToSignalProcess.h"

#include <cmath>
#include <map>
#include <stdexcept>

#include "TRestDetectorHitsEvent.h"
#include "TRestDetectorReadout.h"
#include "TRestDetectorReadoutManager.h"
#include "TRestDetectorSignalEvent.h"
#include "TRestManager.h"
#include "TRestRun.h"

using ROOT::Math::XYZVector;

static const bool TRestDetectorHitsToSignalProcess_FieldsRegistered = []() {
    auto& reg = TRestMetadataFieldRegistry::Instance();
    using P = TRestDetectorHitsToSignalProcess;
    reg.RegisterField<P>("readout", &P::fReadoutName);
    reg.RegisterField<P>("sampling", &P::fSampling);
    return true;
}();

namespace {
const bool kRegistered = []() {
    MetadataClassRegistry::Instance().Register(
        "TRestDetectorHitsToSignalProcess", [](const std::string& instanceName, const YAML::Node& params) {
            return std::make_unique<TRestDetectorHitsToSignalProcess>(instanceName, params);
        });
    return true;
}();
}  // namespace

TRestDetectorHitsToSignalProcess::TRestDetectorHitsToSignalProcess() : TRestEventProcess() {
    fName = "TRestDetectorHitsToSignalProcess";
}

TRestDetectorHitsToSignalProcess::TRestDetectorHitsToSignalProcess(const std::string& instanceName,
                                                                   const YAML::Node& node)
    : TRestEventProcess(instanceName, node) {
    LoadConfig();
}

TRestDetectorHitsToSignalProcess::TRestDetectorHitsToSignalProcess(const std::string& fileName,
                                                                   const std::string& sectionName)
    : TRestEventProcess(fileName, sectionName) {
    LoadConfig();
}

void TRestDetectorHitsToSignalProcess::LoadConfig() {
    TRestEventProcess::LoadConfig();
    if (!fNode || fNode.IsNull()) {
        RESTError << "TRestDetectorHitsToSignalProcess::LoadConfig - YAML node is missing" << RESTendl;
        return;
    }
    UpdateParamsFromYAML<TRestDetectorHitsToSignalProcess>(fNode);
    UpdateYAMLFromParams<TRestDetectorHitsToSignalProcess>(fNode);
}

void TRestDetectorHitsToSignalProcess::InitProcess() {
    if (fSampling < 0) throw std::runtime_error("TRestDetectorHitsToSignalProcess: sampling must be >= 0");
    if (!fManager) throw std::runtime_error("TRestDetectorHitsToSignalProcess: no TRestManager");

    auto* readouts =
        dynamic_cast<TRestDetectorReadoutManager*>(fManager->GetMetadataClass("TRestDetectorReadoutManager"));
    if (!readouts) {
        throw std::runtime_error("TRestDetectorHitsToSignalProcess: no TRestDetectorReadoutManager found");
    }
    fReadout = readouts->SelectReadout(fReadoutName, "TRestDetectorHitsToSignalProcess");
    // The stored metadata keeps the readout used
    UpdateYAMLFromParams<TRestDetectorHitsToSignalProcess>(fNode);

    RegisterObservable("NHitsRejected", fNHitsRejected);
}

bool TRestDetectorHitsToSignalProcess::ProcessEvent(const TRestEvent& input, TRestEvent& output) {
    const auto* hitsEvent = dynamic_cast<const TRestDetectorHitsEvent*>(&input);
    auto* signalEvent = dynamic_cast<TRestDetectorSignalEvent*>(&output);
    if (!hitsEvent || !signalEvent) {
        throw std::runtime_error(
            "TRestDetectorHitsToSignalProcess: input/output is not a TRestDetectorHitsEvent/"
            "TRestDetectorSignalEvent");
    }

    signalEvent->Initialize();
    signalEvent->SetEventInfo(hitsEvent);

    // Charge of each channel versus the time, sorted
    std::map<int, std::map<double, double>> charges;
    fNHitsRejected = 0;
    // Consecutive hits at the same position (e.g. the photoelectrons of a channel) share its channel
    XYZVector lastPosition;
    int lastChannel = -1;
    bool hasLast = false;
    for (size_t n = 0; n < hitsEvent->GetNumberOfHits(); n++) {
        const XYZVector position = hitsEvent->GetPosition(n);
        double time = hitsEvent->GetTime(n);
        if (!hasLast || position != lastPosition) {
            lastChannel = fReadout->GetChannelFromPosition(position.X(), position.Y(), position.Z());
            lastPosition = position;
            hasLast = true;
        }
        const int daqID = lastChannel;
        if (daqID < 0) {
            fNHitsRejected++;
            continue;
        }
        if (fSampling > 0) time = std::floor(time / fSampling) * fSampling;
        charges[daqID][time] += hitsEvent->GetEnergy(n);
    }

    for (const auto& [daqID, points] : charges) {
        signalEvent->AddSignal(daqID, std::vector<std::pair<double, double>>(points.begin(), points.end()));
    }

    RESTDebug << "TRestDetectorHitsToSignalProcess: " << signalEvent->GetNumberOfSignals() << " signals, "
              << fNHitsRejected << " hits rejected" << RESTendl;
    return signalEvent->GetNumberOfSignals() > 0;
}
