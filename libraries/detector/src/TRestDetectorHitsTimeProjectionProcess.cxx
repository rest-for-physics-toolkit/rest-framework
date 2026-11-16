/// \class TRestDetectorHitsTimeProjectionProcess
/// Time projection, as in a time projection chamber (TPC): reconstructs the coordinate along the drift of
/// each hit from its time, the inverse of TRestDetectorHitsDriftProcess: the hit is placed where its distance
/// to the readout (TRestDetectorReadout::GetDistanceToReadout, e.g. the Micromegas plane, positive towards
/// its cathode as set by its normal) is driftVelocity * time, along the axis of the readout normal (e.g. Y
/// for a Micromegas with normal [0, -1, 0]). The other coordinates are kept, also the undefined ones (NaN,
/// e.g. Y of an X strip): an X hit of TRestDetectorSignalToHitsProcess becomes an XZ hit for a readout
/// normal along Z. The energy, time and particle of the hits are kept. Hits whose readout defines no distance
/// (e.g. a readout without drift) are kept as they are and counted in the observable NHitsWithoutDrift.
///
/// The drift velocity is the one of the TRestDetectorMedium with `charge` at its drift field, unless
/// driftVelocity is given (e.g. an effective velocity tuned with data); the one used is written to the
/// metadata. The time origin is that of the hits: a trigger delay must be subtracted beforehand.
///
/// Differences with legacy: the drift velocity and the drift coordinate were part of
/// TRestDetectorSignalToHitsProcess (driftVelocity, plane position + normal * distance).
///
/// Parameters: readout (empty: the data readout of the TRestDetectorReadoutManager), driftVelocity (0: the
/// one of the medium), medium (name; empty: the only medium with `charge`; only if driftVelocity <= 0).
///
/// Example:
///
///    signalToHits:
///      class: "TRestDetectorSignalToHitsProcess"
///      method: "tripleMaxAverage"
///    timeProjection:
///      class: "TRestDetectorHitsTimeProjectionProcess"
///      driftVelocity: 16.5 mm/us

#include "TRestDetectorHitsTimeProjectionProcess.h"

#include <cmath>
#include <stdexcept>

#include "TRestDetectorHitsEvent.h"
#include "TRestDetectorMedium.h"
#include "TRestDetectorReadout.h"
#include "TRestDetectorReadoutManager.h"
#include "TRestManager.h"

using ROOT::Math::XYZVector;

static const bool TRestDetectorHitsTimeProjectionProcess_FieldsRegistered = []() {
    auto& reg = TRestMetadataFieldRegistry::Instance();
    using P = TRestDetectorHitsTimeProjectionProcess;
    reg.RegisterField<P>("readout", &P::fReadoutName);
    reg.RegisterField<P>("medium", &P::fMediumName);
    reg.RegisterField<P>("driftVelocity", &P::fDriftVelocity);
    return true;
}();

namespace {
const bool kRegistered = []() {
    MetadataClassRegistry::Instance().Register(
        "TRestDetectorHitsTimeProjectionProcess",
        [](const std::string& instanceName, const YAML::Node& params) {
            return std::make_unique<TRestDetectorHitsTimeProjectionProcess>(instanceName, params);
        });
    return true;
}();
}  // namespace

TRestDetectorHitsTimeProjectionProcess::TRestDetectorHitsTimeProjectionProcess() : TRestEventProcess() {
    fName = "TRestDetectorHitsTimeProjectionProcess";
}

TRestDetectorHitsTimeProjectionProcess::TRestDetectorHitsTimeProjectionProcess(
    const std::string& instanceName, const YAML::Node& node)
    : TRestEventProcess(instanceName, node) {
    LoadConfig();
}

TRestDetectorHitsTimeProjectionProcess::TRestDetectorHitsTimeProjectionProcess(const std::string& fileName,
                                                                               const std::string& sectionName)
    : TRestEventProcess(fileName, sectionName) {
    LoadConfig();
}

void TRestDetectorHitsTimeProjectionProcess::LoadConfig() {
    TRestEventProcess::LoadConfig();
    if (!fNode || fNode.IsNull()) {
        RESTError << "TRestDetectorHitsTimeProjectionProcess::LoadConfig - YAML node is missing" << RESTendl;
        return;
    }
    UpdateParamsFromYAML<TRestDetectorHitsTimeProjectionProcess>(fNode);
    UpdateYAMLFromParams<TRestDetectorHitsTimeProjectionProcess>(fNode);
}

void TRestDetectorHitsTimeProjectionProcess::InitProcess() {
    const std::string user = "TRestDetectorHitsTimeProjectionProcess";
    if (!fManager) throw std::runtime_error(user + ": no TRestManager");

    auto* readouts =
        dynamic_cast<TRestDetectorReadoutManager*>(fManager->GetMetadataClass("TRestDetectorReadoutManager"));
    if (!readouts) throw std::runtime_error(user + ": no TRestDetectorReadoutManager found");
    fReadout = readouts->SelectReadout(fReadoutName, user);

    if (fDriftVelocity <= 0) {
        const auto* medium = TRestDetectorMedium::Find(*fManager, fMediumName, "charge", user);
        if (!medium) {
            throw std::runtime_error(user +
                                     ": set 'driftVelocity' or define a TRestDetectorMedium with charge");
        }
        fDriftVelocity = medium->GetChargeResponse(user).GetDriftVelocity();
        if (fDriftVelocity <= 0) {
            throw std::runtime_error(user + ": the medium '" + medium->GetName() +
                                     "' gives no drift velocity");
        }
        fMediumName = medium->GetName();
    }
    // The stored metadata keeps the readout and the drift velocity used
    UpdateYAMLFromParams<TRestDetectorHitsTimeProjectionProcess>(fNode);

    RegisterObservable("NHitsWithoutDrift", fNHitsWithoutDrift);
}

bool TRestDetectorHitsTimeProjectionProcess::ProcessEvent(const TRestEvent& input, TRestEvent& output) {
    const auto* inputEvent = dynamic_cast<const TRestDetectorHitsEvent*>(&input);
    auto* outputEvent = dynamic_cast<TRestDetectorHitsEvent*>(&output);
    if (!inputEvent || !outputEvent) {
        throw std::runtime_error(
            "TRestDetectorHitsTimeProjectionProcess: input/output is not a TRestDetectorHitsEvent");
    }

    outputEvent->Initialize();
    outputEvent->SetEventInfo(inputEvent);
    fNHitsWithoutDrift = 0;

    for (size_t n = 0; n < inputEvent->GetNumberOfHits(); n++) {
        XYZVector position = inputEvent->GetPosition(n);
        const double time = inputEvent->GetTime(n);
        // The distance to the readout is a projection on its normal: the undefined coordinates (NaN) are
        // evaluated at 0, they do not change it along an axis normal
        const XYZVector defined(std::isnan(position.X()) ? 0 : position.X(),
                                std::isnan(position.Y()) ? 0 : position.Y(),
                                std::isnan(position.Z()) ? 0 : position.Z());
        const double distance = fReadout->GetDistanceToReadout(defined);
        const XYZVector normal = fReadout->GetReadoutNormal(defined);
        if (std::isnan(distance) || normal.R() == 0) {
            fNHitsWithoutDrift++;
            outputEvent->AddHit(position, inputEvent->GetEnergy(n), time, inputEvent->GetType(n),
                                inputEvent->GetParticleID(n));
            continue;
        }
        // Axis of the normal (the readout normals are along an axis), moved by the missing distance
        const double nx = std::abs(normal.X()), ny = std::abs(normal.Y()), nz = std::abs(normal.Z());
        const double shift = fDriftVelocity * time - distance;
        if (nx >= ny && nx >= nz) {
            position.SetX(defined.X() + shift / normal.X());
        } else if (ny >= nz) {
            position.SetY(defined.Y() + shift / normal.Y());
        } else {
            position.SetZ(defined.Z() + shift / normal.Z());
        }
        outputEvent->AddHit(position, inputEvent->GetEnergy(n), time,
                            TRestHitsData::GetHitType(position.X(), position.Y(), position.Z()),
                            inputEvent->GetParticleID(n));
    }

    return outputEvent->GetNumberOfHits() > 0;
}
