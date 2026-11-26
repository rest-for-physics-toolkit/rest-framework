/// \class TRestRawSignalPeaksFinderProcess
/// Finds all the peaks of each signal of a TRestRawSignalEvent with TRestPulseShapeAnalysis::GetPeaks and
/// groups the peaks of the event in time coincidence windows, as TRestDetectorSignalPeaksFinderProcess
/// does for detector signals. The threshold of each signal is relative to its own noise: `sigmaThreshold`
/// times the sigma of its baseline (computed in `baselineRange` with `baselineOption`, and subtracted),
/// and at least `minThreshold`. So channels with different gains (e.g. the vetoes, a scintillator and a
/// photomultiplier each) are treated alike without calibrating them, and their amplitudes keep the whole
/// range of the ADC. The detector process, on detector signals (baseline subtracted, often zero
/// suppressed), has no noise to refer to and takes a threshold in signal units.
///
/// Times are in bins of the raw signals. The threshold applies to the signal smoothed with a moving average
/// of 2 * halfWindow + 1 points (its noise is smaller than the sigma of the baseline), and the amplitude of a
/// peak is the average of the points (not smoothed) within amplitudeHalfWindow of it. A pulse within the
/// baseline range raises the sigma of that signal; the robust option (median and IQR) is the least affected.
/// Replaces, with TRestDetectorSignalPeaksFinderProcess, the legacy TRestRawPeaksFinderProcess: its TPC
/// signals had a threshold of sigmaOverBaseline times the mean sigma of all the TPC signals, its veto
/// signals a threshold in ADC (thresholdOverBaseline, with halfWindow 2, maxGreaterEqual 0, refineMaximum
/// false and amplitudeHalfWindow 0), and the signals were split by their type in the readout; here the
/// signals of each detector are split beforehand with TRestRawSignalChannelSelectionProcess.
///
/// Observables, the peaks sorted by time and then by signal ID:
/// - PeaksSignalID, PeaksTime (bins), PeaksAmplitude (over the baseline), PeaksBaselineSigma (sigma of
///   the baseline of the signal), PeaksWindowIndex: one entry per peak.
/// - WindowTime (time of the first peak), WindowMultiplicity (number of peaks), WindowAmplitude (sum of
///   the peak amplitudes): one entry per window.
/// - NPeaks, NPeaksSignals (signals with at least one peak), NWindows, PeaksAmplitudeSum.
///
/// Parameters: baselineRange ([5, 55] bins), baselineOption (robust), sigmaThreshold (5), minThreshold (0),
/// distance (10 bins), halfWindow (5), maxGreaterEqual (3), refineMaximum (true), amplitudeHalfWindow (1),
/// window (20 bins).
///
/// Example (vetoes, peaks coincident with the trigger: Sum(vetoPeaks_PeaksTime > 185 &&
/// vetoPeaks_PeaksTime < 225) > 0):
///
///    vetoPeaks:
///      class: "TRestRawSignalPeaksFinderProcess"
///      baselineRange: [5, 55]
///      sigmaThreshold: 5
///      distance: 10
///      window: 20

#include "TRestRawSignalPeaksFinderProcess.h"

#include <algorithm>
#include <set>
#include <stdexcept>
#include <tuple>

#include "TRestPulseShapeAnalysis.h"
#include "TRestRawSignalEvent.h"
#include "TRestTools.h"

static const bool TRestRawSignalPeaksFinderProcess_FieldsRegistered = []() {
    auto& reg = TRestMetadataFieldRegistry::Instance();
    using P = TRestRawSignalPeaksFinderProcess;
    reg.RegisterField<P>("baselineRange", &P::fBaselineRange);
    reg.RegisterField<P>("baselineOption", &P::fBaselineOption);
    reg.RegisterField<P>("sigmaThreshold", &P::fSigmaThreshold);
    reg.RegisterField<P>("minThreshold", &P::fMinThreshold);
    reg.RegisterField<P>("distance", &P::fDistance);
    reg.RegisterField<P>("halfWindow", &P::fHalfWindow);
    reg.RegisterField<P>("maxGreaterEqual", &P::fMaxGreaterEqual);
    reg.RegisterField<P>("refineMaximum", &P::fRefineMaximum);
    reg.RegisterField<P>("amplitudeHalfWindow", &P::fAmplitudeHalfWindow);
    reg.RegisterField<P>("window", &P::fWindow);
    return true;
}();

namespace {
const bool kRegistered = []() {
    MetadataClassRegistry::Instance().Register(
        "TRestRawSignalPeaksFinderProcess", [](const std::string& instanceName, const YAML::Node& params) {
            return std::make_unique<TRestRawSignalPeaksFinderProcess>(instanceName, params);
        });
    return true;
}();
}  // namespace

TRestRawSignalPeaksFinderProcess::TRestRawSignalPeaksFinderProcess() : TRestEventProcess() {
    fName = "TRestRawSignalPeaksFinderProcess";
}

TRestRawSignalPeaksFinderProcess::TRestRawSignalPeaksFinderProcess(const std::string& instanceName,
                                                                   const YAML::Node& node)
    : TRestEventProcess(instanceName, node) {
    LoadConfig();
}

TRestRawSignalPeaksFinderProcess::TRestRawSignalPeaksFinderProcess(const std::string& fileName,
                                                                   const std::string& sectionName)
    : TRestEventProcess(fileName, sectionName) {
    LoadConfig();
}

void TRestRawSignalPeaksFinderProcess::LoadConfig() {
    TRestEventProcess::LoadConfig();

    if (!fNode || fNode.IsNull()) {
        RESTError << "TRestRawSignalPeaksFinderProcess::LoadConfig - YAML node is missing" << RESTendl;
        return;
    }

    UpdateParamsFromYAML<TRestRawSignalPeaksFinderProcess>(fNode);
    fBaselineOption = TRestTools::NormalizeOption("baselineOption", fBaselineOption,
                                                  TRestPulseShapeAnalysis::kBaselineOptions);
    UpdateYAMLFromParams<TRestRawSignalPeaksFinderProcess>(fNode);
}

void TRestRawSignalPeaksFinderProcess::InitProcess() {
    if (fBaselineRange.first < 0 || fBaselineRange.second <= fBaselineRange.first + 1) {
        throw std::runtime_error(
            "TRestRawSignalPeaksFinderProcess: baselineRange must be [first, second) with at least two bins");
    }
    if (fSigmaThreshold < 0 || fMinThreshold < 0 || fDistance < 0 || fWindow < 0 || fHalfWindow < 1 ||
        fMaxGreaterEqual < 0 || fAmplitudeHalfWindow < 0) {
        throw std::runtime_error(
            "TRestRawSignalPeaksFinderProcess: sigmaThreshold, minThreshold, distance, window, "
            "maxGreaterEqual and amplitudeHalfWindow must be >= 0 and halfWindow >= 1");
    }

    RegisterObservable("PeaksSignalID", fObs.PeaksSignalID);
    RegisterObservable("PeaksTime", fObs.PeaksTime);
    RegisterObservable("PeaksAmplitude", fObs.PeaksAmplitude);
    RegisterObservable("PeaksBaselineSigma", fObs.PeaksBaselineSigma);
    RegisterObservable("PeaksWindowIndex", fObs.PeaksWindowIndex);

    RegisterObservable("WindowTime", fObs.WindowTime);
    RegisterObservable("WindowMultiplicity", fObs.WindowMultiplicity);
    RegisterObservable("WindowAmplitude", fObs.WindowAmplitude);

    RegisterObservable("NPeaks", fObs.NPeaks);
    RegisterObservable("NPeaksSignals", fObs.NPeaksSignals);
    RegisterObservable("NWindows", fObs.NWindows);
    RegisterObservable("PeaksAmplitudeSum", fObs.PeaksAmplitudeSum);
}

bool TRestRawSignalPeaksFinderProcess::ProcessEvent(const TRestEvent& input, TRestEvent& output) {
    output.CopyFrom(&input);

    const auto* signalEvent = dynamic_cast<const TRestRawSignalEvent*>(&input);
    if (!signalEvent) {
        throw std::runtime_error("TRestRawSignalPeaksFinderProcess: input is not a TRestRawSignalEvent");
    }

    fObs.clear();

    // (time, signal ID, amplitude, sigma of the baseline) of all the peaks of the event
    std::vector<std::tuple<double, int, double, double>> peaks;
    std::set<int> peakSignals;
    for (int s = 0; s < signalEvent->GetNumberOfSignals(); s++) {
        const auto& signal = signalEvent->GetSignal(s);
        double baseline = 0, sigma = 0;
        const auto data = TRestPulseShapeAnalysis::CalculateBaselineAndSigma(
            signal.GetData(), fBaselineRange.first, fBaselineRange.second, baseline, sigma, fBaselineOption);
        if (data.empty()) continue;

        std::vector<std::pair<double, double>> points(data.size());
        for (size_t bin = 0; bin < data.size(); bin++) points[bin] = {static_cast<double>(bin), data[bin]};
        const double threshold = std::max(fSigmaThreshold * sigma, fMinThreshold);
        for (const auto& [time, amplitude] :
             TRestPulseShapeAnalysis::GetPeaks(points, threshold, fDistance, fHalfWindow, fMaxGreaterEqual,
                                               fRefineMaximum, fAmplitudeHalfWindow)) {
            peaks.emplace_back(time, signal.GetSignalID(), amplitude, sigma);
            peakSignals.insert(signal.GetSignalID());
        }
    }
    std::sort(peaks.begin(), peaks.end());

    for (const auto& [time, signalID, amplitude, sigma] : peaks) {
        // Peaks are sorted by time, so a peak opens a new window if it is too far from the current one
        if (fObs.WindowTime.empty() || time - fObs.WindowTime.back() > fWindow) {
            fObs.WindowTime.push_back(time);
            fObs.WindowMultiplicity.push_back(0);
            fObs.WindowAmplitude.push_back(0);
        }
        fObs.WindowMultiplicity.back()++;
        fObs.WindowAmplitude.back() += amplitude;

        fObs.PeaksSignalID.push_back(signalID);
        fObs.PeaksTime.push_back(time);
        fObs.PeaksAmplitude.push_back(amplitude);
        fObs.PeaksBaselineSigma.push_back(sigma);
        fObs.PeaksWindowIndex.push_back(fObs.WindowTime.size() - 1);
        fObs.PeaksAmplitudeSum += amplitude;
    }

    fObs.NPeaks = peaks.size();
    fObs.NPeaksSignals = peakSignals.size();
    fObs.NWindows = fObs.WindowTime.size();

    return true;
}
