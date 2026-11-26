#ifndef TRESTRAWSIGNALPEAKSFINDERPROCESS_H
#define TRESTRAWSIGNALPEAKSFINDERPROCESS_H

#include <string>
#include <utility>
#include <vector>

#include "TRestEventProcess.h"

/// \class TRestRawSignalPeaksFinderProcess
/// \brief Finds all the peaks of the signals of a TRestRawSignalEvent, over a threshold in sigmas of the
/// baseline of each signal, and groups them in time coincidence windows. The event is not modified. See
/// the class documentation in TRestRawSignalPeaksFinderProcess.cxx (observables and example).
///
/// \author Luis Obis <luis.antonio.obis@gmail.com>
/// \author JuanAn García <juanangp@unizar.es>
class TRestRawSignalPeaksFinderProcess : public TRestEventProcess {
   private:
    struct PeaksObservables {
        std::vector<int> PeaksSignalID;
        std::vector<double> PeaksTime;
        std::vector<double> PeaksAmplitude;
        std::vector<double> PeaksBaselineSigma;
        std::vector<int> PeaksWindowIndex;

        std::vector<double> WindowTime;
        std::vector<int> WindowMultiplicity;
        std::vector<double> WindowAmplitude;

        int NPeaks = 0;
        int NPeaksSignals = 0;
        int NWindows = 0;
        double PeaksAmplitudeSum = 0.0;

        void clear() {
            auto clear_vectors = [](auto&... vecs) { (vecs.clear(), ...); };
            clear_vectors(PeaksSignalID, PeaksTime, PeaksAmplitude, PeaksBaselineSigma, PeaksWindowIndex,
                          WindowTime, WindowMultiplicity, WindowAmplitude);
            NPeaks = NPeaksSignals = NWindows = 0;
            PeaksAmplitudeSum = 0.0;
        }
    };

    PeaksObservables fObs;

   public:
    /// Bins [first, second) where the baseline and its sigma of each signal are computed
    std::pair<int, int> fBaselineRange = {5, 55};

    /// Baseline estimator: "" (mean and standard deviation), robust (median and IQR) or outliers
    std::string fBaselineOption = "robust";

    /// Threshold of the smoothed signal, in sigmas of the baseline of the signal
    double fSigmaThreshold = 5;

    /// Minimum threshold (ADC over the baseline), e.g. for signals without noise (simulations)
    double fMinThreshold = 0;

    /// Minimum distance between two peaks of the same signal (bins)
    double fDistance = 10;

    /// Half of the number of points of the moving average used to find the peaks
    int fHalfWindow = 5;

    /// Maximum number of points of the moving average window greater or equal than a peak
    int fMaxGreaterEqual = 3;

    /// Moves the peak to the highest smoothed point within distance after it
    bool fRefineMaximum = true;

    /// The peak amplitude is the average of the points within this number of points of the peak
    int fAmplitudeHalfWindow = 1;

    /// A coincidence window starts at the earliest peak not yet in a window and contains all the peaks up to
    /// window bins after it
    double fWindow = 20;

    /// \brief Default constructor (no configuration).
    TRestRawSignalPeaksFinderProcess();
    /// \brief Builds the process from an in-memory YAML section (LoadConfig).
    TRestRawSignalPeaksFinderProcess(const std::string& instanceName, const YAML::Node& node);
    /// \brief Builds the process from the section sectionName of a YAML file (LoadConfig).
    TRestRawSignalPeaksFinderProcess(const std::string& fileName, const std::string& sectionName);

    /// \brief Reads the parameters from the YAML node.
    void LoadConfig() override;
    /// \brief Checks the parameters and registers the observables.
    void InitProcess() override;
    /// \brief Finds the peaks of the signals and fills the observables; the event is copied unchanged.
    /// \return Always true.
    bool ProcessEvent(const TRestEvent& input, TRestEvent& output) override;
    /// \brief Nothing to do at the end of the processing.
    void EndProcess() override {}

    /// \brief Input event type, "TRestRawSignalEvent".
    std::string GetInputEvent() const override { return "TRestRawSignalEvent"; }
    /// \brief Output event type, "TRestRawSignalEvent".
    std::string GetOutputEvent() const override { return "TRestRawSignalEvent"; }

    /// \brief Returns the class name, "TRestRawSignalPeaksFinderProcess".
    std::string GetClassName() const override { return "TRestRawSignalPeaksFinderProcess"; }
};

#endif
