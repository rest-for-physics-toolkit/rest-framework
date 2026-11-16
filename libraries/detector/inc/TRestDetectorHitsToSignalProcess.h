#ifndef TRESTDETECTORHITSTOSIGNALPROCESS_H
#define TRESTDETECTORHITSTOSIGNALPROCESS_H

#include <string>

#include "TRestEventProcess.h"

class TRestDetectorReadout;

/// \class TRestDetectorHitsToSignalProcess
/// \brief Collects the hits of a TRestDetectorHitsEvent in the readout channels: a TRestDetectorSignalEvent
/// with the charge (energy) of each channel versus the time of its hits (generic: no drift).
///
/// The parameters (YAML keys), the observables and an example are described in the class documentation
/// (TRestDetectorHitsToSignalProcess.cxx).
///
/// \author Javier Galán <javier.galan@unizar.es>
/// \author JuanAn García <juanangp@unizar.es>
class TRestDetectorHitsToSignalProcess : public TRestEventProcess {
   private:
    const TRestDetectorReadout* fReadout = nullptr;
    /// Hits without channel in the event
    int fNHitsRejected = 0;  //!

   public:
    /// Name of the readout (in TRestDetectorReadoutManager). It can be omitted if the manager contains a
    /// single readout.
    std::string fReadoutName = "";
    /// Time bin of the signals: the times are rounded down to multiples of it (0: exact times)
    TRestWithUnits fSampling = 0;  // us

    /// \brief Default constructor (no configuration).
    TRestDetectorHitsToSignalProcess();
    /// \brief Builds the process from an in-memory YAML section (LoadConfig).
    TRestDetectorHitsToSignalProcess(const std::string& instanceName, const YAML::Node& node);
    /// \brief Builds the process from the section sectionName of a YAML file (LoadConfig).
    TRestDetectorHitsToSignalProcess(const std::string& fileName, const std::string& sectionName);

    /// \brief Reads the parameters from the YAML node.
    void LoadConfig() override;
    /// \brief Gets the readout, checks sampling (>= 0) and registers the observable NHitsRejected.
    void InitProcess() override;
    /// \brief Adds the energy of each hit to the signal of its channel at its time (see the class
    /// description).
    /// \return false if the event has no signals.
    bool ProcessEvent(const TRestEvent& input, TRestEvent& output) override;
    /// \brief Nothing to do at the end of the processing.
    void EndProcess() override {}

    /// \brief Input event type, "TRestDetectorHitsEvent".
    std::string GetInputEvent() const override { return "TRestDetectorHitsEvent"; }
    /// \brief Output event type, "TRestDetectorSignalEvent".
    std::string GetOutputEvent() const override { return "TRestDetectorSignalEvent"; }

    /// \brief Returns the class name, "TRestDetectorHitsToSignalProcess".
    std::string GetClassName() const override { return "TRestDetectorHitsToSignalProcess"; }
};

#endif
