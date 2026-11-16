/// \class TRestGeant4ToDetectorHitsProcess
/// Transfers the hits (position, energy, time) of a TRestGeant4Event with positive energy to a
/// TRestDetectorHitsEvent as XYZ hits, with the Particle ID (PDG encoding) of the track of each hit
/// (fHitParticleID, e.g. for the scintillation profile of each particle). Only the hits of the volumes
/// selected as in TRestGeant4QuenchingProcess are transferred
/// (TRestGeant4GeometryInfo::GetVolumesMatchingExpression: physical volumes matching each regular expression
/// or, if none, those of the matching logical volumes);
/// without `volumes`, those of the sensitive volumes of the simulation (`volumes: [".*"]` or `["^."]`
/// for all of them).
///
/// The times of the hits are relative to the first energy deposit of the Geant4 event in any volume
/// (TRestGeant4Event::GetFirstDepositTime, the same origin as the times of TRestGeant4AnalysisProcess and for
/// every detector of the event): the times of a radioactive decay, ~1e13 us after the generation of the
/// event, become those of the detector (a trigger at the first deposit, drift and scintillation times of ns
/// or us). The time of that first deposit is the observable TimeOffset (the event time is not changed: it is
/// the generation time of the event).
///
/// Difference with legacy: there is no `type` per volume (legacy VETO hits), the hits of different
/// detectors are transferred by separate processes or pipelines with their own `volumes`.
///
/// Parameters: volumes ([], the sensitive volumes), list of regular expressions (e.g. "gas" for any name
/// containing gas, "^gas" for the names starting with it).
///
/// Example:
///
///    geant4ToHits:
///      class: "TRestGeant4ToDetectorHitsProcess"
///      volumes:
///        - "gasVolume"
///        - "^scintillatorVolume"

#include "TRestGeant4ToDetectorHitsProcess.h"

#include <stdexcept>

#include "TRestDetectorHitsEvent.h"
#include "TRestGeant4Event.h"
#include "TRestGeant4Metadata.h"
#include "TRestRun.h"

static const bool TRestGeant4ToDetectorHitsProcess_FieldsRegistered = []() {
    auto& reg = TRestMetadataFieldRegistry::Instance();
    using P = TRestGeant4ToDetectorHitsProcess;
    reg.RegisterField<P>("volumes", &P::fVolumes);
    return true;
}();

namespace {
const bool kRegistered = []() {
    MetadataClassRegistry::Instance().Register(
        "TRestGeant4ToDetectorHitsProcess", [](const std::string& instanceName, const YAML::Node& params) {
            return std::make_unique<TRestGeant4ToDetectorHitsProcess>(instanceName, params);
        });
    return true;
}();
}  // namespace

TRestGeant4ToDetectorHitsProcess::TRestGeant4ToDetectorHitsProcess() : TRestEventProcess() {
    fName = "TRestGeant4ToDetectorHitsProcess";
}

TRestGeant4ToDetectorHitsProcess::TRestGeant4ToDetectorHitsProcess(const std::string& instanceName,
                                                                   const YAML::Node& node)
    : TRestEventProcess(instanceName, node) {
    LoadConfig();
}

TRestGeant4ToDetectorHitsProcess::TRestGeant4ToDetectorHitsProcess(const std::string& fileName,
                                                                   const std::string& sectionName)
    : TRestEventProcess(fileName, sectionName) {
    LoadConfig();
}

void TRestGeant4ToDetectorHitsProcess::LoadConfig() {
    TRestEventProcess::LoadConfig();

    if (!fNode || fNode.IsNull()) {
        RESTError << "TRestGeant4ToDetectorHitsProcess::LoadConfig - YAML node is missing" << RESTendl;
        return;
    }

    UpdateParamsFromYAML<TRestGeant4ToDetectorHitsProcess>(fNode);
    UpdateYAMLFromParams<TRestGeant4ToDetectorHitsProcess>(fNode);
}

void TRestGeant4ToDetectorHitsProcess::InitProcess() {
    fVolumeIDs.clear();

    const auto* metadata =
        fRunInfo ? dynamic_cast<TRestGeant4Metadata*>(fRunInfo->GetMetadataClass("TRestGeant4Metadata"))
                 : nullptr;
    if (!metadata) {
        throw std::runtime_error(
            "TRestGeant4ToDetectorHitsProcess: no TRestGeant4Metadata found, needed to select volumes");
    }
    const auto& geometryInfo = metadata->GetGeant4GeometryInfo();

    if (fVolumes.empty()) {
        const auto& sensitive = metadata->GetSensitiveVolumes();
        fVolumeIDs = geometryInfo.GetVolumeIDsFromNames({sensitive.begin(), sensitive.end()});
        if (fVolumeIDs.empty()) {
            throw std::runtime_error(
                "TRestGeant4ToDetectorHitsProcess: no 'volumes' and no sensitive volume found");
        }
    }

    // Volume selection as TRestGeant4QuenchingProcess: the names matching each expression, compared
    // with the physical and filtered names of the hit volumes (resolved here once per volume ID)
    for (const auto& expression : fVolumes) {
        const auto ids = geometryInfo.GetVolumeIDsMatchingExpression(expression);
        if (ids.empty()) {
            RESTWarning << "TRestGeant4ToDetectorHitsProcess: no volume matches '" << expression << "'"
                        << RESTendl;
        }
        fVolumeIDs.insert(ids.begin(), ids.end());
    }
    for (const auto id : fVolumeIDs) {
        RESTDebug << "TRestGeant4ToDetectorHitsProcess: volume " << geometryInfo.GetVolumeFromID(id)
                  << " (id " << id << ") selected" << RESTendl;
    }

    if (fVolumeIDs.empty()) {
        throw std::runtime_error("TRestGeant4ToDetectorHitsProcess: no volume matches 'volumes'");
    }

    RegisterObservable("TimeOffset", fTimeOffset);
}

bool TRestGeant4ToDetectorHitsProcess::ProcessEvent(const TRestEvent& input, TRestEvent& output) {
    const auto* g4Event = dynamic_cast<const TRestGeant4Event*>(&input);
    auto* hitsEvent = dynamic_cast<TRestDetectorHitsEvent*>(&output);
    if (!g4Event || !hitsEvent) {
        throw std::runtime_error(
            "TRestGeant4ToDetectorHitsProcess: input/output is not a TRestGeant4Event/"
            "TRestDetectorHitsEvent");
    }

    hitsEvent->Initialize();
    hitsEvent->SetEventInfo(&input);

    const auto& data = g4Event->fEventData;
    const auto& hits = data.hitsStorage;
    // Times relative to the first deposit of the event
    fTimeOffset = g4Event->GetFirstDepositTime();
    // The hits of the track t are [trackStartIndices[t], trackStartIndices[t] + trackNHits[t])
    for (size_t t = 0; t < data.trackParticleIDs.size(); t++) {
        const size_t start = data.trackStartIndices.at(t);
        const size_t end = start + data.trackNHits.at(t);
        for (size_t i = start; i < end; i++) {
            if (hits.energy[i] <= 0) continue;
            if (fVolumeIDs.count(data.hitVolumeID.at(i)) == 0) continue;
            hitsEvent->AddHit(hits.x[i], hits.y[i], hits.z[i], hits.energy[i], hits.time[i] - fTimeOffset,
                              TRestDetectorHitsEvent::REST_HitType::XYZ, data.trackParticleIDs[t]);
        }
    }

    RESTDebug << "TRestGeant4ToDetectorHitsProcess: " << hitsEvent->GetNumberOfHits() << " hits added"
              << RESTendl;

    return hitsEvent->GetNumberOfHits() > 0;
}
