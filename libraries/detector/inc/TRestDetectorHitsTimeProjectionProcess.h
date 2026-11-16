#ifndef TRESTDETECTORHITSTIMEPROJECTIONPROCESS_H
#define TRESTDETECTORHITSTIMEPROJECTIONPROCESS_H

#include <string>

#include "TRestEventProcess.h"

class TRestDetectorReadout;

/// \class TRestDetectorHitsTimeProjectionProcess
/// \brief Reconstructs the coordinate along the drift of the hits of a TRestDetectorHitsEvent from their
/// time: the distance to the readout is driftVelocity * time (the inverse of TRestDetectorHitsDriftProcess).
///
/// The parameters (YAML keys), the observables and an example are described in the class documentation
/// (TRestDetectorHitsTimeProjectionProcess.cxx).
///
/// \author Javier Galán <javier.galan@unizar.es>
/// \author JuanAn García <juanangp@unizar.es>
class TRestDetectorHitsTimeProjectionProcess : public TRestEventProcess {
   private:
    const TRestDetectorReadout* fReadout = nullptr;
    /// Hits whose readout defines no distance (kept as they are) in the event
    int fNHitsWithoutDrift = 0;  //!

   public:
    /// Name of the readout (in TRestDetectorReadoutManager) that gives the drift geometry (plane and
    /// normal). It can be omitted if the manager contains a single readout.
    std::string fReadoutName = "";
    /// Name of the TRestDetectorMedium that gives the drift velocity (at its drift field) when driftVelocity
    /// <= 0. Empty: the only medium with `charge`.
    std::string fMediumName = "";
    /// Drift velocity (position units per time unit). If > 0 it overrides the one of the medium (e.g. an
    /// effective velocity tuned with data); written with the one used.
    TRestWithUnits fDriftVelocity = 0;  // mm/us

    /// \brief Default constructor (no configuration).
    TRestDetectorHitsTimeProjectionProcess();
    /// \brief Builds the process from an in-memory YAML section (LoadConfig).
    TRestDetectorHitsTimeProjectionProcess(const std::string& instanceName, const YAML::Node& node);
    /// \brief Builds the process from the section sectionName of a YAML file (LoadConfig).
    TRestDetectorHitsTimeProjectionProcess(const std::string& fileName, const std::string& sectionName);

    /// \brief Reads the parameters from the YAML node.
    void LoadConfig() override;
    /// \brief Gets the readout and, if driftVelocity <= 0, the drift velocity of the medium (written to the
    /// metadata); registers the observable NHitsWithoutDrift. Throws if there is no drift velocity.
    void InitProcess() override;
    /// \brief Sets the drift coordinate of each hit from its time (see the class description).
    /// \return false if the event has no hits.
    bool ProcessEvent(const TRestEvent& input, TRestEvent& output) override;
    /// \brief Nothing to do at the end of the processing.
    void EndProcess() override {}

    /// \brief Input event type, "TRestDetectorHitsEvent".
    std::string GetInputEvent() const override { return "TRestDetectorHitsEvent"; }
    /// \brief Output event type, "TRestDetectorHitsEvent".
    std::string GetOutputEvent() const override { return "TRestDetectorHitsEvent"; }

    /// \brief Returns the class name, "TRestDetectorHitsTimeProjectionProcess".
    std::string GetClassName() const override { return "TRestDetectorHitsTimeProjectionProcess"; }
};

#endif
