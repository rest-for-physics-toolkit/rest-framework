#ifndef TRESTDETECTORPHOTOSENSORPROCESS_H
#define TRESTDETECTORPHOTOSENSORPROCESS_H

#include <map>
#include <string>
#include <utility>

#include "Math/Vector3D.h"
#include "TRandom3.h"
#include "TRestEventProcess.h"

class TRestDetectorHitsEvent;
class TRestDetectorReadout;

/// \class TRestDetectorPhotosensorProcess
/// \brief Response of the photosensors (SiPM, PMT) to the photoelectrons of a TRestDetectorHitsEvent (a hit
/// per photoelectron, TRestDetectorLightCollectionProcess): dark counts, optical crosstalk, afterpulses,
/// saturation of the cells of a SiPM, the time jitter and the gain spread. A hit per avalanche, to be
/// collected by TRestDetectorHitsToSignalProcess.
///
/// The parameters (YAML keys), the observables and an example are described in the class documentation
/// (TRestDetectorPhotosensorProcess.cxx).
///
/// \author JuanAn García <juanangp@unizar.es>
class TRestDetectorPhotosensorProcess : public TRestEventProcess {
   private:
    const TRestDetectorReadout* fReadout = nullptr;
    /// Position of each channel given to its dark counts (dark counts only)
    std::map<int, ROOT::Math::XYZVector> fChannelPositions;  //!
    TRandom3 fRandom;                                        //!

    /// Observables
    int fNAvalanches = 0;        //!
    int fNDarkCounts = 0;        //!
    int fNCrosstalk = 0;         //!
    int fNAfterpulses = 0;       //!
    int fNSaturationLosses = 0;  //!

    /// \brief Delay of an afterpulse: Exp(afterpulseDelay), or Gauss(afterpulseDelay, afterpulseSpread)
    /// truncated at 0.
    double AfterpulseDelay();

   public:
    /// Name of the readout (in TRestDetectorReadoutManager) of the photosensors. It can be omitted if the
    /// manager contains a single readout.
    std::string fReadoutName = "";

    /// `timeJitter`: standard deviation of the Gaussian time jitter of the photosensor (transit time spread,
    /// 0: none)
    TRestWithUnits fTimeJitter = 0;  // us

    /// `gainSpread`: relative spread (sigma / mean) of the charge of an avalanche, the gain fluctuation of
    /// the photosensor (0: charge 1)
    double fGainSpread = 0;

    /// `darkCountRate`: rate of dark counts of each channel (e.g. 100 kHz for a SiPM, 1 kHz for a PMT; 0:
    /// none)
    TRestWithUnits fDarkCountRate = 0;  // MHz

    /// `darkCountWindow`: time range of the dark counts, e.g. [-2 us, 18 us] (the times of the photoelectrons
    /// are relative to the first energy deposit of the event)
    std::pair<TRestWithUnits, TRestWithUnits> fDarkCountWindow = {0, 0};  // us

    /// `crosstalkProbability`: probability that an avalanche triggers another one in a neighbouring cell at
    /// the same time (SiPM optical crosstalk; 0 for a PMT)
    double fCrosstalkProbability = 0;

    /// `afterpulseProbability`: probability that an avalanche is followed by an afterpulse
    double fAfterpulseProbability = 0;

    /// `afterpulseDelay`: mean delay of the afterpulses (exponential, e.g. tens of ns for a SiPM), or the
    /// position of the peak with afterpulseSpread (e.g. a few us for the ions of a PMT)
    TRestWithUnits fAfterpulseDelay = 0;  // us

    /// `afterpulseSpread`: standard deviation of a Gaussian delay of the afterpulses (0: exponential)
    TRestWithUnits fAfterpulseSpread = 0;  // us

    /// `cells`: number of cells (microcells) of the SiPM of a channel (0: no saturation, e.g. a PMT)
    int fCells = 0;

    /// `cellRecoveryTime`: time after an avalanche during which its cell gives no other (0: a cell fires
    /// once per event)
    TRestWithUnits fCellRecoveryTime = 0;  // us

    /// Seed of the random generator (0: random seed). Otherwise the generator is re-seeded at each event
    /// with TRestEvent::EventSeed (seed, event and subevent IDs): independent of the threads.
    unsigned long fSeed = 0;

    /// \brief Default constructor (no configuration).
    TRestDetectorPhotosensorProcess();
    /// \brief Builds the process from an in-memory YAML section (LoadConfig).
    TRestDetectorPhotosensorProcess(const std::string& instanceName, const YAML::Node& node);
    /// \brief Builds the process from the section sectionName of a YAML file (LoadConfig).
    TRestDetectorPhotosensorProcess(const std::string& fileName, const std::string& sectionName);

    /// \brief Reads the parameters; throws if a value is out of range (negative times, rates or spreads,
    /// probabilities not in [0, 1), crosstalkProbability + afterpulseProbability >= 1, darkCountWindow
    /// empty with dark counts, afterpulseDelay not positive with afterpulses, cells negative).
    void LoadConfig() override;
    /// \brief Gets the readout and, with dark counts, the position of each channel; registers the
    /// observables.
    void InitProcess() override;
    /// \brief Adds the dark counts, the crosstalk and the afterpulses of each channel, removes the
    /// avalanches of busy cells and gives each avalanche its jitter and charge.
    /// \return false if the event has no avalanches.
    bool ProcessEvent(const TRestEvent& input, TRestEvent& output) override;
    /// \brief Nothing to do at the end of the processing.
    void EndProcess() override {}

    /// \brief Input event type, "TRestDetectorHitsEvent".
    std::string GetInputEvent() const override { return "TRestDetectorHitsEvent"; }
    /// \brief Output event type, "TRestDetectorHitsEvent".
    std::string GetOutputEvent() const override { return "TRestDetectorHitsEvent"; }

    /// \brief Returns the class name, "TRestDetectorPhotosensorProcess".
    std::string GetClassName() const override { return "TRestDetectorPhotosensorProcess"; }
};

#endif
