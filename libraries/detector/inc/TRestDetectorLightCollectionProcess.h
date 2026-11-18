#ifndef TRESTDETECTORLIGHTCOLLECTIONPROCESS_H
#define TRESTDETECTORLIGHTCOLLECTIONPROCESS_H

#include <functional>
#include <string>

#include "Math/Vector3D.h"
#include "TRandom3.h"
#include "TRestEventProcess.h"

class TRestDetectorHitsEvent;
class TRestDetectorLightResponse;
class TRestDetectorWLSResponse;
class TRestDetectorReadout;

/// \class TRestDetectorLightCollectionProcess
/// \brief Collection of the scintillation photons of the hits (TRestDetectorScintillationProcess) by the
/// readout channels and their detection: one hit per photoelectron, with its time (scintillation profile of
/// the particle, transport and re-emission in the fibres, jitter of the photosensor) and its charge (gain
/// spread of the photosensor), to be collected by TRestDetectorHitsToSignalProcess. A model for each kind of
/// scintillator: `localCapture` (opaque scintillator, the channels near the photons
/// capture them with a probability that depends on the distance, e.g. the fibres of LiquidO) or
/// `attenuation` (transparent scintillator, the photons reach the readout of their channel attenuated and
/// delayed with the distance, e.g. a veto with its PMT).
///
/// The parameters (YAML keys), the observables and an example are described in the class documentation
/// (TRestDetectorLightCollectionProcess.cxx).
///
/// \author Luis Obis <luis.antonio.obis@gmail.com>
/// \author JuanAn García <juanangp@unizar.es>
class TRestDetectorLightCollectionProcess : public TRestEventProcess {
   private:
    const TRestDetectorReadout* fReadout = nullptr;
    const TRestDetectorLightResponse* fLight = nullptr;
    const TRestDetectorWLSResponse* fWLS = nullptr;
    /// Capture probability versus the distance (mm), built from captureProbability (localCapture)
    std::function<double(double)> fCapture;  //!
    TRandom3 fRandom;                        //!

    /// Observables
    double fNPhotoelectrons = 0;   //!
    int fNHitsWithoutChannel = 0;  //!

    /// \brief Photoelectrons of the photons reaching a photosensor (Binomial(photons, quantumEfficiency)).
    int Detect(int photons);
    /// \brief Adds a hit per photoelectron at a channel position: the arrival time of its photon plus its
    /// emission time (scintillation profile of the particle), the re-emission of the wavelength shifter and
    /// the jitter of the photosensor, with charge 1 or Gauss(1, gainSpread).
    void AddPhotoelectrons(const ROOT::Math::XYZVector& position, int photoelectrons, double time,
                           int particleID, TRestDetectorHitsEvent& output);

    /// \brief The photoelectrons of the photons of a hit captured by the channels in range (localCapture),
    /// with the particle of the hit.
    void CaptureLocally(const ROOT::Math::XYZVector& position, int photons, double time, int particleID,
                        TRestDetectorHitsEvent& output);
    /// \brief The photoelectrons of the photons of a hit reaching the readout of its channel (attenuation),
    /// with the particle of the hit.
    void Attenuate(const ROOT::Math::XYZVector& position, int photons, double time, int particleID,
                   TRestDetectorHitsEvent& output);

   public:
    /// `model`: "localCapture" or "attenuation"
    std::string fModel = "localCapture";

    /// Name of the readout (in TRestDetectorReadoutManager) that collects the light. It can be omitted if
    /// the manager contains a single readout.
    std::string fReadoutName = "";

    /// Name of the TRestDetectorMedium with the light response: the scintillation time profiles and, for
    /// attenuation, attenuationLength and effectiveLightSpeed. Empty: the only medium with `light` (none for
    /// localCapture: no emission time).
    std::string fMediumName = "";

    /// `captureProbability` (localCapture): probability that a photon produced at a distance r (mm) from a
    /// channel is captured by it and reaches its photosensor, a ROOT formula of x = r (e.g.
    /// "0.03*exp(-x*x/200)") or a graph or histogram of a ROOT file, "file.root:name" (e.g. from an optical
    /// simulation of Geant4). With wlsMedium, the probability that the fibre of the channel captures it (the
    /// guiding to the photosensor is then given by the `wls` response)
    std::string fCaptureProbability = "";

    /// `wlsMedium` (localCapture): name of the TRestDetectorMedium of the fibres, whose `wls` response gives
    /// the attenuation and delay of the captured light along the fibre (half in each direction, the end not
    /// read reflecting endReflectivity of it) to its read end and the re-emission time. Empty: none.
    std::string fWLSMediumName = "";

    /// `maxDistance` (localCapture): channels farther than this collect no light
    TRestWithUnits fMaxDistance = 0;

    /// `quantumEfficiency`: probability that a photon reaching the photosensor (SiPM, PMT) of a channel gives
    /// a photoelectron
    double fQuantumEfficiency = 1;

    /// `timeJitter`: standard deviation of the Gaussian time jitter of the photosensor (transit time spread,
    /// 0: none)
    TRestWithUnits fTimeJitter = 0;  // us

    /// `gainSpread`: relative spread (sigma / mean) of the charge of a photoelectron, the gain fluctuation of
    /// the photosensor (0: charge 1)
    double fGainSpread = 0;

    /// Seed of the random generator (0: random seed). Otherwise the generator is re-seeded at each event
    /// with TRestDetectorHitsEvent::EventSeed (seed, event and subevent IDs): independent of the threads.
    unsigned long fSeed = 0;

    /// \brief Default constructor (no configuration).
    TRestDetectorLightCollectionProcess();
    /// \brief Builds the process from an in-memory YAML section (LoadConfig).
    TRestDetectorLightCollectionProcess(const std::string& instanceName, const YAML::Node& node);
    /// \brief Builds the process from the section sectionName of a YAML file (LoadConfig).
    TRestDetectorLightCollectionProcess(const std::string& fileName, const std::string& sectionName);

    /// \brief Reads the parameters; throws if the model is unknown, quantumEfficiency is not in [0, 1],
    /// timeJitter or gainSpread is negative, or for localCapture if maxDistance is not positive or
    /// captureProbability is not a valid formula or graph.
    void LoadConfig() override;
    /// \brief Gets the readout, the light response of the medium (needed for attenuation) and the wls
    /// response of wlsMedium, and registers the observables NPhotoelectrons and NHitsWithoutChannel.
    void InitProcess() override;
    /// \brief Distributes the photons of each hit among the channels collecting them and converts them into
    /// photoelectrons, a hit each.
    /// \return false if the event has no photoelectrons.
    bool ProcessEvent(const TRestEvent& input, TRestEvent& output) override;
    /// \brief Nothing to do at the end of the processing.
    void EndProcess() override {}

    /// \brief Capture probability at a distance (mm) from a channel (localCapture; 0 beyond maxDistance).
    double GetCaptureProbability(double distance) const;

    /// \brief Input event type, "TRestDetectorHitsEvent".
    std::string GetInputEvent() const override { return "TRestDetectorHitsEvent"; }
    /// \brief Output event type, "TRestDetectorHitsEvent".
    std::string GetOutputEvent() const override { return "TRestDetectorHitsEvent"; }

    /// \brief Returns the class name, "TRestDetectorLightCollectionProcess".
    std::string GetClassName() const override { return "TRestDetectorLightCollectionProcess"; }
};

#endif
