#ifndef TRESTDETECTORLIGHTTOSIGNALPROCESS_H
#define TRESTDETECTORLIGHTTOSIGNALPROCESS_H

#include <string>

#include "TRandom3.h"
#include "TRestEventProcess.h"

class TRestDetectorLightResponse;
class TRestDetectorWLSResponse;
class TRestDetectorReadout;

/// \class TRestDetectorLightToSignalProcess
/// \brief Photoelectrons of the channels (hits of TRestDetectorLightCollectionProcess) to signals: each
/// photoelectron is added to the signal of its channel at the time of its hit plus the emission time of the
/// scintillation of its particle (profile of the medium), the re-emission time of a wavelength shifter and
/// the time jitter of the photosensor.
///
/// The parameters (YAML keys), the observables and an example are described in the class documentation
/// (TRestDetectorLightToSignalProcess.cxx).
///
/// \author JuanAn García <juanangp@unizar.es>
class TRestDetectorLightToSignalProcess : public TRestEventProcess {
   private:
    const TRestDetectorReadout* fReadout = nullptr;
    const TRestDetectorLightResponse* fLight = nullptr;
    const TRestDetectorWLSResponse* fWLS = nullptr;
    TRandom3 fRandom;  //!

    /// Observables
    double fNPhotoelectrons = 0;   //!
    int fNHitsWithoutChannel = 0;  //!

   public:
    /// Name of the readout (in TRestDetectorReadoutManager). It can be omitted if the manager contains a
    /// single readout.
    std::string fReadoutName = "";
    /// Name of the TRestDetectorMedium with the scintillation profiles (`light`). It can be omitted if the
    /// manager contains a single medium.
    std::string fMediumName = "";
    /// Name of the TRestDetectorMedium with the wavelength shifting response (`wls`, e.g. the WLS fibres)
    /// that delays each photoelectron by its re-emission time. Empty: no wavelength shifter.
    std::string fWLSMediumName = "";
    /// Standard deviation of the Gaussian time jitter of the photosensor (transit time spread, 0: none)
    TRestWithUnits fTimeJitter = 0;  // us
    /// Relative spread (sigma / mean) of the charge of a photoelectron, the gain fluctuation of the
    /// photosensor (0: every photoelectron has charge 1)
    double fGainSpread = 0;
    /// Time bin of the signals: the times are rounded down to multiples of it (0: exact times)
    TRestWithUnits fSampling = 0;  // us
    /// Seed of the random generator (0: random seed). Otherwise the generator is re-seeded at each event
    /// with TRestDetectorHitsEvent::EventSeed (seed, event and subevent IDs): independent of the threads.
    unsigned long fSeed = 0;

    /// \brief Default constructor (no configuration).
    TRestDetectorLightToSignalProcess();
    /// \brief Builds the process from an in-memory YAML section (LoadConfig).
    TRestDetectorLightToSignalProcess(const std::string& instanceName, const YAML::Node& node);
    /// \brief Builds the process from the section sectionName of a YAML file (LoadConfig).
    TRestDetectorLightToSignalProcess(const std::string& fileName, const std::string& sectionName);

    /// \brief Reads the parameters from the YAML node.
    void LoadConfig() override;
    /// \brief Gets the readout, the light response of the medium and the wls response of wlsMedium, checks
    /// timeJitter, gainSpread and sampling (>= 0) and registers the observables NPhotoelectrons and
    /// NHitsWithoutChannel.
    void InitProcess() override;
    /// \brief Adds each photoelectron to the signal of its channel at its time (see the class description).
    /// \return false if the event has no signals.
    bool ProcessEvent(const TRestEvent& input, TRestEvent& output) override;
    /// \brief Nothing to do at the end of the processing.
    void EndProcess() override {}

    /// \brief Input event type, "TRestDetectorHitsEvent".
    std::string GetInputEvent() const override { return "TRestDetectorHitsEvent"; }
    /// \brief Output event type, "TRestDetectorSignalEvent".
    std::string GetOutputEvent() const override { return "TRestDetectorSignalEvent"; }

    /// \brief Returns the class name, "TRestDetectorLightToSignalProcess".
    std::string GetClassName() const override { return "TRestDetectorLightToSignalProcess"; }
};

#endif
