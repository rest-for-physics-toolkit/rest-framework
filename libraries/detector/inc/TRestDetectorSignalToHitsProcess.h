#ifndef TRESTDETECTORSIGNALTOHITSPROCESS_H
#define TRESTDETECTORSIGNALTOHITSPROCESS_H

#include <string>
#include <utility>
#include <vector>

#include "TRestEventProcess.h"

class TRestDetectorReadout;

/// \class TRestDetectorSignalToHitsProcess
/// \brief Transforms the signals of a TRestDetectorSignalEvent into the hits of a TRestDetectorHitsEvent
/// using the readout (channel position) and the pulses (energy and time); generic, no drift.
///
/// The parameters (YAML keys), the observables and an example are described in the class documentation
/// (TRestDetectorSignalToHitsProcess.cxx).
///
/// \author Javier Galán <javier.galan@unizar.es>
/// \author JuanAn García <juanangp@unizar.es>
class TRestDetectorSignalToHitsProcess : public TRestEventProcess {
   private:
    TRestDetectorReadout* fReadout = nullptr;

   public:
    /// Methods accepted by fMethod. onlyMax, tripleMax, tripleMaxAverage, gaussFit, landauFit and agetFit
    /// are delegated to TRestPulseShapeAnalysis::GetMax
    inline static const std::vector<std::string> kMethods = {
        "onlyMax", "tripleMax", "tripleMaxAverage", "gaussFit", "landauFit",
        "agetFit", "integral",  "qCenter",          "all",      "intWindow"};

    /// Name of the readout (in TRestDetectorReadoutManager) that decodes the signal channels. It can be
    /// omitted if the manager contains a single readout.
    std::string fReadoutName = "";

    /// Method used to extract the hits from each signal: onlyMax, tripleMax, tripleMaxAverage, gaussFit,
    /// landauFit, agetFit, integral, qCenter (one hit per signal), all (one hit per point) or intWindow
    /// (one hit per integration window).
    std::string fMethod = "tripleMaxAverage";

    /// Points of the moving average applied to the signal before taking its maximum with onlyMax,
    /// tripleMax and tripleMaxAverage (odd, 1 or less: no smoothing), so the noise does not spread the
    /// times; the energy is then the maximum of the smoothed signal
    int fSmoothingPoints = 5;

    /// Time window to integrate when method is intWindow
    TRestWithUnits fIntWindow = 5;

    /// Minimum integral of a window to create a hit when method is intWindow
    double fThreshold = 100;

    /// Time range integrated when method is integral. An empty range (first >= second) integrates the
    /// whole signal.
    std::pair<double, double> fIntegralRange = {0, 0};

    /// \brief Default constructor (no configuration).
    TRestDetectorSignalToHitsProcess();
    /// \brief Builds the process from an in-memory YAML section (LoadConfig).
    TRestDetectorSignalToHitsProcess(const std::string& instanceName, const YAML::Node& node);
    /// \brief Builds the process from the section sectionName of a YAML file (LoadConfig).
    TRestDetectorSignalToHitsProcess(const std::string& fileName, const std::string& sectionName);

    /// \brief Reads the parameters from the YAML node.
    void LoadConfig() override;
    /// \brief Gets the readout (written to the metadata).
    void InitProcess() override;
    /// \brief Creates the hits of the signals with the chosen method.
    /// \return false if no hit is created.
    bool ProcessEvent(const TRestEvent& input, TRestEvent& output) override;
    /// \brief Nothing to do at the end of the processing.
    void EndProcess() override {}

    /// \brief Input event type, "TRestDetectorSignalEvent".
    std::string GetInputEvent() const override { return "TRestDetectorSignalEvent"; }
    /// \brief Output event type, "TRestDetectorHitsEvent".
    std::string GetOutputEvent() const override { return "TRestDetectorHitsEvent"; }

    /// \brief Returns the class name, "TRestDetectorSignalToHitsProcess".
    std::string GetClassName() const override { return "TRestDetectorSignalToHitsProcess"; }
};

#endif
