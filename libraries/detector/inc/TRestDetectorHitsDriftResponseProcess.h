#ifndef TRESTDETECTORHITSDRIFTRESPONSEPROCESS_H
#define TRESTDETECTORHITSDRIFTRESPONSEPROCESS_H

#include <memory>
#include <string>
#include <vector>

#include "TRandom3.h"
#include "TRestDetectorDriftResponse.h"
#include "TRestEventProcess.h"

class TRestDetectorChargeResponse;
class TRestDetectorReadout;

/// \class TRestDetectorHitsDriftResponseProcess
/// \brief Converts the energy deposits of a TRestDetectorHitsEvent into ionization electrons and drifts them
/// to the boundary plane above the readout with a drift response computed by restGarfield in a field map
/// (TRestDetectorDriftResponse): the drift in a non-uniform field, for special simulations, instead of
/// TRestDetectorHitsDriftProcess (uniform field).
///
/// The parameters (YAML keys) and an example are described in the class documentation
/// (TRestDetectorHitsDriftResponseProcess.cxx).
///
/// \author JuanAn García <juanangp@unizar.es>
class TRestDetectorHitsDriftResponseProcess : public TRestEventProcess {
   private:
    const TRestDetectorReadout* fReadout = nullptr;
    const TRestDetectorChargeResponse* fCharge = nullptr;
    std::shared_ptr<const TRestDetectorDriftResponse> fResponse;  //!
    TRandom3 fRandom;                                             //!

   public:
    /// Values accepted by fFluctuations
    inline static const std::vector<std::string> kFluctuations = {"fano", "poisson", "none"};
    /// `responseFile`: file of the drift response (restGarfield, `drift` key)
    std::string fResponseFile = "";
    /// Name of the readout (in TRestDetectorReadoutManager): its normal and the projection of the hits. It
    /// can be omitted if the manager contains a single readout.
    std::string fReadoutName = "";
    /// Name of the TRestDetectorMedium (with `charge`: wValue, fanoFactor). Empty: the only medium with
    /// `charge`.
    std::string fMediumName = "";
    /// Fluctuation of the number of electrons: fano (gaussian with the medium Fano factor), poisson or none.
    std::string fFluctuations = "fano";
    /// Maximum number of electrons per event (0: no limit). Above it, each electron carries the energy of
    /// several (an effective W-value), as TRestDetectorHitsDriftProcess.
    int fMaxHits = 1000;
    /// Seed of the random generator (0: random seed). Otherwise the generator is re-seeded at each event
    /// with TRestEvent::EventSeed (seed, event and subevent IDs): independent of the threads.
    unsigned long fSeed = 0;

    /// \brief Default constructor (no configuration).
    TRestDetectorHitsDriftResponseProcess();
    /// \brief Builds the process from an in-memory YAML section (LoadConfig).
    TRestDetectorHitsDriftResponseProcess(const std::string& instanceName, const YAML::Node& node);
    /// \brief Builds the process from the section sectionName of a YAML file (LoadConfig).
    TRestDetectorHitsDriftResponseProcess(const std::string& fileName, const std::string& sectionName);

    /// \brief Reads the parameters from the YAML node.
    void LoadConfig() override;
    /// \brief Gets the readout, the charge response of the medium (wValue > 0, and fanoFactor > 0 with fano
    /// fluctuations) and reads the drift response of responseFile. Throws if one of them is missing.
    void InitProcess() override;
    /// \brief Converts the hits into drifted electrons (see the class description).
    /// \return false if no electron survives.
    bool ProcessEvent(const TRestEvent& input, TRestEvent& output) override;
    /// \brief Nothing to do at the end of the processing.
    void EndProcess() override {}

    /// \brief Input event type, "TRestDetectorHitsEvent".
    std::string GetInputEvent() const override { return "TRestDetectorHitsEvent"; }
    /// \brief Output event type, "TRestDetectorHitsEvent".
    std::string GetOutputEvent() const override { return "TRestDetectorHitsEvent"; }

    /// \brief Returns the class name, "TRestDetectorHitsDriftResponseProcess".
    std::string GetClassName() const override { return "TRestDetectorHitsDriftResponseProcess"; }
};

#endif
