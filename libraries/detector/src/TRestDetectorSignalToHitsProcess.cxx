/// \class TRestDetectorSignalToHitsProcess
/// Creates the hits of an event from its signals. It is generic: the channel of each signal gives the
/// coordinates defined by the readout (TRestDetectorReadout::GetPositionFromChannel, e.g. only X or only Y
/// for strips, the centre of a scintillator; undefined ones are NaN) and the pulse gives the energy and the
/// time of the hits. The hit type (X, Y, XY, ...) is deduced from the coordinates that are defined. In a
/// detector with drift, the coordinate along the drift is obtained from the time afterwards, by
/// TRestDetectorHitsTimeProjectionProcess.
///
/// Methods (energy, time of each hit):
/// - onlyMax, tripleMax, tripleMaxAverage, gaussFit, landauFit, agetFit: one hit at the pulse maximum.
///   With smoothingPoints > 1, onlyMax, tripleMax and tripleMaxAverage take the maximum of the signal
///   smoothed with a moving average of that number of points (TRestPulseShapeAnalysis::GetSignalSmoothed;
///   the fits already average the noise), with the energy of the smoothed signal.
/// - integral: one hit with the charge integrated in integralRange (whole signal if empty) at the time
///   of the maximum in that range. For pulses without drift (scintillators, solid state, bolometers).
/// - qCenter: one hit at the charge-weighted mean time, energy = mean charge per point (as legacy).
/// - all: one hit per signal point.
/// - intWindow: one hit per integration window of intWindow, if the mean charge is over threshold.
///
/// Differences with legacy: no driftVelocity (TRestDetectorHitsTimeProjectionProcess).
///
/// Parameters: readout (empty: the data readout of the TRestDetectorReadoutManager), method
/// ("tripleMaxAverage"), smoothingPoints (5), intWindow (5), threshold (100), integralRange ([0, 0]).
///
/// Example:
///
///    signalToHits:
///      class: "TRestDetectorSignalToHitsProcess"
///      method: "tripleMaxAverage"

#include "TRestDetectorSignalToHitsProcess.h"

#include <cmath>
#include <stdexcept>

#include "TRestDetectorHitsEvent.h"
#include "TRestDetectorReadout.h"
#include "TRestDetectorReadoutManager.h"
#include "TRestDetectorSignalEvent.h"
#include "TRestManager.h"
#include "TRestPulseShapeAnalysis.h"
#include "TRestTools.h"

using namespace TRestPulseShapeAnalysis;

static const bool TRestDetectorSignalToHitsProcess_FieldsRegistered = []() {
    auto& reg = TRestMetadataFieldRegistry::Instance();
    using P = TRestDetectorSignalToHitsProcess;
    reg.RegisterField<P>("readout", &P::fReadoutName);
    reg.RegisterField<P>("method", &P::fMethod);
    reg.RegisterField<P>("smoothingPoints", &P::fSmoothingPoints);
    reg.RegisterField<P>("intWindow", &P::fIntWindow);
    reg.RegisterField<P>("threshold", &P::fThreshold);
    reg.RegisterField<P>("integralRange", &P::fIntegralRange);
    return true;
}();

namespace {
const bool kRegistered = []() {
    MetadataClassRegistry::Instance().Register(
        "TRestDetectorSignalToHitsProcess", [](const std::string& instanceName, const YAML::Node& params) {
            return std::make_unique<TRestDetectorSignalToHitsProcess>(instanceName, params);
        });
    return true;
}();
}  // namespace

TRestDetectorSignalToHitsProcess::TRestDetectorSignalToHitsProcess() : TRestEventProcess() {
    fName = "TRestDetectorSignalToHitsProcess";
}

TRestDetectorSignalToHitsProcess::TRestDetectorSignalToHitsProcess(const std::string& instanceName,
                                                                   const YAML::Node& node)
    : TRestEventProcess(instanceName, node) {
    LoadConfig();
}

TRestDetectorSignalToHitsProcess::TRestDetectorSignalToHitsProcess(const std::string& fileName,
                                                                   const std::string& sectionName)
    : TRestEventProcess(fileName, sectionName) {
    LoadConfig();
}

void TRestDetectorSignalToHitsProcess::LoadConfig() {
    TRestEventProcess::LoadConfig();

    if (!fNode || fNode.IsNull()) {
        RESTError << "TRestDetectorSignalToHitsProcess::LoadConfig - YAML node is missing" << RESTendl;
        return;
    }

    UpdateParamsFromYAML<TRestDetectorSignalToHitsProcess>(fNode);
    fMethod = TRestTools::NormalizeOption("method", fMethod, kMethods);
    if (fMethod == "intWindow" && !(fIntWindow > 0)) {
        throw std::runtime_error("TRestDetectorSignalToHitsProcess: method intWindow needs intWindow > 0");
    }
    UpdateYAMLFromParams<TRestDetectorSignalToHitsProcess>(fNode);
}

void TRestDetectorSignalToHitsProcess::InitProcess() {
    auto* readouts = fManager ? dynamic_cast<TRestDetectorReadoutManager*>(
                                    fManager->GetMetadataClass("TRestDetectorReadoutManager"))
                              : nullptr;
    if (!readouts) {
        throw std::runtime_error("TRestDetectorSignalToHitsProcess: no TRestDetectorReadoutManager found");
    }
    fReadout = readouts->SelectReadout(fReadoutName, "TRestDetectorSignalToHitsProcess");

    // The stored metadata keeps the readout used, not the default of the YAML
    UpdateYAMLFromParams<TRestDetectorSignalToHitsProcess>(fNode);
}

bool TRestDetectorSignalToHitsProcess::ProcessEvent(const TRestEvent& input, TRestEvent& output) {
    const auto* signalEvent = dynamic_cast<const TRestDetectorSignalEvent*>(&input);
    auto* hitsEvent = dynamic_cast<TRestDetectorHitsEvent*>(&output);
    if (!signalEvent || !hitsEvent) {
        throw std::runtime_error(
            "TRestDetectorSignalToHitsProcess: input/output is not a TRestDetectorSignalEvent/"
            "TRestDetectorHitsEvent");
    }

    hitsEvent->Initialize();
    hitsEvent->SetEventInfo(&input);

    for (int s = 0; s < signalEvent->GetNumberOfSignals(); s++) {
        const auto& signal = signalEvent->GetSignal(s);
        if (signal.GetNPoints() == 0) continue;
        const int daqID = signal.GetSignalID();

        auto channelPos = fReadout->GetPositionFromChannel(daqID);
        if (std::isnan(channelPos.X()) && std::isnan(channelPos.Y()) && std::isnan(channelPos.Z())) {
            RESTDebug << "TRestDetectorSignalToHitsProcess: channel " << daqID << " not found in readout"
                      << RESTendl;
            continue;
        }

        // (time, energy) pairs of the hits created by this signal
        std::vector<std::pair<double, double>> hits;
        auto graph = signal.GetGraph();
        if (fMethod == "intWindow") {
            for (const auto& hit : GetIntWindow(graph, fIntWindow)) {
                if (hit.second > fThreshold) hits.push_back(hit);
            }
        } else if (fMethod == "integral") {
            const bool wholeSignal = fIntegralRange.first >= fIntegralRange.second;
            // Negative times are valid (e.g. a trigger delay): a flag marks the maximum found
            double charge = 0, maxCharge = 0, maxTime = 0;
            bool found = false;
            for (int j = 0; j < signal.GetNPoints(); j++) {
                const auto [time, value] = signal.GetPoint(j);
                if (!wholeSignal && (time < fIntegralRange.first || time > fIntegralRange.second)) continue;
                charge += value;
                if (!found || value > maxCharge) {
                    maxCharge = value;
                    maxTime = time;
                    found = true;
                }
            }
            if (found && charge > 0) hits.emplace_back(maxTime, charge);
        } else if (fMethod == "qCenter") {
            double charge = 0, weightedTime = 0;
            for (int j = 0; j < signal.GetNPoints(); j++) {
                const auto [time, value] = signal.GetPoint(j);
                charge += value;
                weightedTime += time * value;
            }
            if (charge > 0) hits.emplace_back(weightedTime / charge, charge / signal.GetNPoints());
        } else if (fMethod == "all") {
            for (int j = 0; j < signal.GetNPoints(); j++) hits.push_back(signal.GetPoint(j));
        } else {
            const bool smooth = fSmoothingPoints > 1 && (fMethod == "onlyMax" || fMethod == "tripleMax" ||
                                                         fMethod == "tripleMaxAverage");
            if (smooth) {
                const auto smoothed = GetSignalSmoothed(signal.GetCharge(), fSmoothingPoints);
                for (int j = 0; j < graph.GetN(); j++) graph.SetPointY(j, smoothed[j]);
            }
            const auto hit = GetMax(graph, fMethod);
            // Failures (empty signal, failed agetFit) return (-1, -1), rejected by the amplitude: negative
            // times are valid. A failed gaussFit or landauFit gives the maximum point
            if (std::isfinite(hit.first) && hit.second > 0) hits.push_back(hit);
        }

        const auto type = TRestHitsData::GetHitType(channelPos.X(), channelPos.Y(), channelPos.Z());
        for (const auto& [time, energy] : hits) hitsEvent->AddHit(channelPos, energy, time, type);
    }

    return hitsEvent->GetNumberOfHits() > 0;
}
