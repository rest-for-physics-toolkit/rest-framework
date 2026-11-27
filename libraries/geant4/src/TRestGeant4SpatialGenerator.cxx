#include "TRestGeant4SpatialGenerator.h"

#include <Math/AxisAngle.h>

#include <algorithm>
#include <atomic>
#include <cmath>
#include <limits>
#include <mutex>
#include <stdexcept>

#include "TRestGeant4CRY.h"
#include "TRestGeant4PrimaryGeneratorTypes.h"
#include "TRestTools.h"

using namespace std;
using namespace TRestGeant4PrimaryGeneratorTypes;
using ROOT::Math::XYZVector;
using TRestTools::ToLower;

static const bool TRestGeant4SpatialGenerator_FieldsRegistered = []() {
    auto& reg = TRestMetadataFieldRegistry::Instance();
    reg.RegisterField<TRestGeant4SpatialGenerator>("type", &TRestGeant4SpatialGenerator::fType);

    using Point = TRestGeant4SpatialGeneratorPoint;
    reg.RegisterField<Point>("position", &Point::fPosition);

    using Shape = TRestGeant4SpatialGeneratorShape;
    reg.RegisterField<Shape>("from", &Shape::fFrom);
    reg.RegisterField<Shape>("shape", &Shape::fShape);
    reg.RegisterField<Shape>("size", &Shape::fSize);
    reg.RegisterField<Shape>("position", &Shape::fPosition);
    reg.RegisterField<Shape>("rotationAxis", &Shape::fRotationAxis);
    reg.RegisterField<Shape>("rotationValue", &Shape::fRotationValue);

    using Cosmic = TRestGeant4SpatialGeneratorCosmic;
    reg.RegisterField<Cosmic>("radius", &Cosmic::fRadius);
    reg.RegisterField<Cosmic>("flux", &Cosmic::fFlux);

    using CRY = TRestGeant4SpatialGeneratorCRY;
    reg.RegisterField<CRY>("position", &CRY::fPosition);
    reg.RegisterField<CRY>("verticalAxis", &CRY::fVerticalAxis);
    reg.RegisterField<CRY>("showersUsed", &CRY::fShowersUsed);
    reg.RegisterField<CRY>("exposureTime", &CRY::fExposureTime);
    return true;
}();

namespace {

XYZVector ToVector(const array<TRestWithUnits, 3>& values) { return {values[0], values[1], values[2]}; }

}  // namespace

// --- TRestGeant4SpatialGenerator ---

unique_ptr<TRestGeant4SpatialGenerator> TRestGeant4SpatialGenerator::Create(const string& name,
                                                                            const YAML::Node& node) {
    if (!node["type"]) throw runtime_error("TRestGeant4SpatialGenerator: the generator 'type' is missing");
    const auto type = ToLower(node["type"].as<string>());
    if (type == "point") return make_unique<TRestGeant4SpatialGeneratorPoint>(name, node);
    if (type == "volume") return make_unique<TRestGeant4SpatialGeneratorVolume>(name, node);
    if (type == "surface") return make_unique<TRestGeant4SpatialGeneratorSurface>(name, node);
    if (type == "cosmic") return make_unique<TRestGeant4SpatialGeneratorCosmic>(name, node);
    if (type == "cry") return make_unique<TRestGeant4SpatialGeneratorCRY>(name, node);
    throw runtime_error("TRestGeant4SpatialGenerator: unknown generator type '" + type +
                        "', it must be point, volume, surface, cosmic or cry");
}

// --- TRestGeant4SpatialGeneratorPoint ---

TRestGeant4SpatialGeneratorPoint::TRestGeant4SpatialGeneratorPoint(const string& name, const YAML::Node& node)
    : TRestGeant4SpatialGenerator(name, node) {
    LoadConfig();
}

void TRestGeant4SpatialGeneratorPoint::LoadConfig() {
    UpdateParamsFromYAML<TRestGeant4SpatialGeneratorPoint>(fNode);
    UpdateYAMLFromParams<TRestGeant4SpatialGeneratorPoint>(fNode);
}

XYZVector TRestGeant4SpatialGeneratorPoint::Sample(TRandom&) { return GetCenter(); }

XYZVector TRestGeant4SpatialGeneratorPoint::GetCenter() const { return ToVector(fPosition); }

// --- TRestGeant4SpatialGeneratorShape ---

XYZVector TRestGeant4SpatialGeneratorShape::GetCenter() const { return ToVector(fPosition); }

void TRestGeant4SpatialGeneratorShape::SetGdmlPlacement(const XYZVector& center, const XYZVector& axis,
                                                        double angle) {
    fPosition = {center.X(), center.Y(), center.Z()};
    fRotationAxis = {axis.X(), axis.Y(), axis.Z()};
    fRotationValue = angle;
}

void TRestGeant4SpatialGeneratorShape::CheckShape(const map<string, size_t>& shapes) {
    // As legacy, a generator from a geometry volume has no shape
    if (!fFrom.empty()) {
        fShape = "";
        return;
    }
    fShape = ToLower(fShape);
    const auto shape = shapes.find(fShape);
    if (shape == shapes.end()) {
        string names;
        for (const auto& [name, nSizes] : shapes) names += (names.empty() ? "" : ", ") + name;
        throw runtime_error(GetClassName() + ": 'from' (a geometry volume) or 'shape' (" + names +
                            ") is needed, '" + fShape + "' is not valid");
    }
    for (size_t i = 0; i < shape->second; i++) {
        if (fSize[i] <= 0) {
            throw runtime_error(GetClassName() + ": the " + fShape + " needs " + to_string(shape->second) +
                                " positive 'size' values");
        }
    }
    if (fRotationValue != 0 && ToVector(fRotationAxis).R() == 0) {
        throw runtime_error(GetClassName() + ": 'rotationAxis' must not be zero");
    }
}

XYZVector TRestGeant4SpatialGeneratorShape::ToWorld(const XYZVector& local) const {
    if (fRotationValue == 0) return local + GetCenter();
    const ROOT::Math::AxisAngle rotation(ToVector(fRotationAxis).Unit(), fRotationValue);
    return rotation * local + GetCenter();
}

// --- TRestGeant4SpatialGeneratorVolume ---

TRestGeant4SpatialGeneratorVolume::TRestGeant4SpatialGeneratorVolume(const string& name,
                                                                     const YAML::Node& node)
    : TRestGeant4SpatialGeneratorShape(name, node) {
    LoadConfig();
}

void TRestGeant4SpatialGeneratorVolume::LoadConfig() {
    UpdateParamsFromYAML<TRestGeant4SpatialGeneratorVolume>(fNode);
    CheckShape({{"box", 3}, {"sphere", 1}, {"cylinder", 2}});
    UpdateYAMLFromParams<TRestGeant4SpatialGeneratorVolume>(fNode);
}

XYZVector TRestGeant4SpatialGeneratorVolume::Sample(TRandom& random) {
    if (IsGdml()) throw runtime_error(GetClassName() + ": a GDML volume needs the Geant4 geometry");
    const double phi = random.Uniform(0, 2 * M_PI);
    XYZVector local;
    if (fShape == "box") {
        local = {fSize[0] * (random.Uniform() - 0.5), fSize[1] * (random.Uniform() - 0.5),
                 fSize[2] * (random.Uniform() - 0.5)};
    } else if (fShape == "sphere") {
        local = fSize[0] * cbrt(random.Uniform()) * IsotropicDirection(random);
    } else {
        const double r = fSize[0] * sqrt(random.Uniform());
        local = {r * cos(phi), r * sin(phi), fSize[1] * (random.Uniform() - 0.5)};
    }
    return ToWorld(local);
}

// --- TRestGeant4SpatialGeneratorSurface ---

TRestGeant4SpatialGeneratorSurface::TRestGeant4SpatialGeneratorSurface(const string& name,
                                                                       const YAML::Node& node)
    : TRestGeant4SpatialGeneratorShape(name, node) {
    LoadConfig();
}

void TRestGeant4SpatialGeneratorSurface::LoadConfig() {
    UpdateParamsFromYAML<TRestGeant4SpatialGeneratorSurface>(fNode);
    CheckShape({{"sphere", 1}, {"cylinder", 2}, {"circle", 1}, {"wall", 2}});
    UpdateYAMLFromParams<TRestGeant4SpatialGeneratorSurface>(fNode);
}

XYZVector TRestGeant4SpatialGeneratorSurface::Sample(TRandom& random) {
    if (IsGdml()) throw runtime_error(GetClassName() + ": a GDML surface needs the Geant4 geometry");
    const double phi = random.Uniform(0, 2 * M_PI);
    XYZVector local;
    if (fShape == "sphere") {
        local = fSize[0] * IsotropicDirection(random);
    } else if (fShape == "cylinder") {
        local = {fSize[0] * cos(phi), fSize[0] * sin(phi), fSize[1] * (random.Uniform() - 0.5)};
    } else if (fShape == "circle") {
        const double r = fSize[0] * sqrt(random.Uniform());
        local = {r * cos(phi), r * sin(phi), 0};
    } else {
        local = {fSize[0] * (random.Uniform() - 0.5), fSize[1] * (random.Uniform() - 0.5), 0};
    }
    return ToWorld(local);
}

// --- TRestGeant4SpatialGeneratorCosmic ---

TRestGeant4SpatialGeneratorCosmic::TRestGeant4SpatialGeneratorCosmic(const string& name,
                                                                     const YAML::Node& node)
    : TRestGeant4SpatialGenerator(name, node) {
    LoadConfig();
}

void TRestGeant4SpatialGeneratorCosmic::LoadConfig() {
    UpdateParamsFromYAML<TRestGeant4SpatialGeneratorCosmic>(fNode);
    if (fRadius < 0) throw runtime_error(GetClassName() + ": 'radius' must be positive");
    if (fFlux < 0) throw runtime_error(GetClassName() + ": 'flux' must be positive");
    UpdateYAMLFromParams<TRestGeant4SpatialGeneratorCosmic>(fNode);
}

XYZVector TRestGeant4SpatialGeneratorCosmic::Sample(TRandom&) {
    throw runtime_error(GetClassName() + ": the cosmic origin depends on the direction");
}

XYZVector TRestGeant4SpatialGeneratorCosmic::SampleFromDirection(TRandom& random,
                                                                 const XYZVector& direction) const {
    if (fRadius <= 0) throw runtime_error(GetClassName() + ": the radius is not set");
    // Two unit vectors normal to the direction
    const XYZVector reference = abs(direction.Z()) < 0.9 ? XYZVector(0, 0, 1) : XYZVector(1, 0, 0);
    const XYZVector v = direction.Cross(reference).Unit();
    const XYZVector w = direction.Cross(v);

    const double r = fRadius * sqrt(random.Uniform());
    const double phi = random.Uniform(0, 2 * M_PI);
    const XYZVector onDisk = r * (cos(phi) * v + sin(phi) * w);
    return onDisk - sqrt(fRadius * fRadius - r * r) * direction;
}

void TRestGeant4SpatialGeneratorCosmic::SetDefaultRadius(double radius) {
    if (fRadius > 0) return;
    fRadius = radius;
    UpdateYAMLFromParams<TRestGeant4SpatialGeneratorCosmic>(fNode);
}

void TRestGeant4SpatialGeneratorCosmic::SetFlux(double flux) {
    fFlux = flux;
    UpdateYAMLFromParams<TRestGeant4SpatialGeneratorCosmic>(fNode);
}

// --- TRestGeant4SpatialGeneratorCRY ---

namespace {
TRestGeant4SpatialGeneratorCRY::EngineFactory& CRYEngineFactory() {
    static TRestGeant4SpatialGeneratorCRY::EngineFactory factory;
    return factory;
}

/// Mass (keV) of the particles of CRY by their Geant4 name, for their speed; others taken at the light speed
double CRYParticleMass(const string& name) {
    if (name == "e-" || name == "e+") return 510.99895;
    if (name == "mu-" || name == "mu+") return 105658.3755;
    if (name == "pi-" || name == "pi+") return 139570.39;
    if (name == "pi0") return 134976.8;
    if (name == "kaon-" || name == "kaon+") return 493677;
    if (name == "neutron") return 939565.42;
    if (name == "proton") return 938272.088;
    return 0;
}

/// Moves a particle to where its line enters the box of half size `half` centered at the origin, its time
/// delayed by the flight (nothing if it starts inside). False if its line misses the box ahead of it
bool MoveToBox(TRestGeant4ParticleState& particle, const XYZVector& half) {
    const XYZVector origin = particle.GetOrigin();
    const XYZVector direction = particle.GetDirection().Unit();
    const double o[3] = {origin.X(), origin.Y(), origin.Z()};
    const double d[3] = {direction.X(), direction.Y(), direction.Z()};
    const double h[3] = {half.X(), half.Y(), half.Z()};
    double enter = -numeric_limits<double>::infinity(), leave = numeric_limits<double>::infinity();
    for (int i = 0; i < 3; i++) {
        if (d[i] == 0) {
            if (abs(o[i]) >= h[i]) return false;
            continue;
        }
        const double t1 = (-h[i] - o[i]) / d[i], t2 = (h[i] - o[i]) / d[i];
        enter = max(enter, min(t1, t2));
        leave = min(leave, max(t1, t2));
    }
    if (leave <= max(enter, 0.0)) return false;
    if (enter <= 0) return true;  // inside the world
    // 1 um inside the world, not on its surface
    const double distance = enter + min(1e-3, (leave - enter) / 2);
    const double mass = CRYParticleMass(particle.GetParticleName());
    const double kinetic = particle.GetEnergy();
    const double beta = mass > 0 ? sqrt(kinetic * (kinetic + 2 * mass)) / (kinetic + mass) : 1;
    constexpr double lightSpeed = 299792.458;  // mm/us
    particle.SetOrigin(origin + distance * direction);
    particle.SetTime(particle.GetTime() + distance / (beta * lightSpeed));
    return true;
}
}  // namespace

namespace {
/// Adds to an atomic double (fetch_add of doubles is C++20)
void AtomicAdd(std::atomic<double>& value, double delta) {
    double current = value.load();
    while (!value.compare_exchange_weak(current, current + delta)) {
    }
}
}  // namespace

/// The CRY generators of the copies (threads), each one its own CRY with the seed `seed` + its index, and
/// their total showers and exposure, shared by the copies
struct TRestGeant4SpatialGeneratorCRY::State {
    std::mutex mutex;  // the list of engines
    std::string setup;
    long seed = 0;
    std::vector<std::shared_ptr<Engine>> engines;
    std::atomic<int64_t> showers{0};
    std::atomic<double> exposure{0};  // us
};

/// The CRY generator of a copy (thread) and its exposure already added to the total
struct TRestGeant4SpatialGeneratorCRY::Engine {
    std::unique_ptr<TRestGeant4CRYEngine> engine;
    double exposure = 0;
};

TRestGeant4SpatialGeneratorCRY::TRestGeant4SpatialGeneratorCRY(const string& name, const YAML::Node& node)
    : TRestGeant4SpatialGenerator(name, node) {
    LoadConfig();
}

void TRestGeant4SpatialGeneratorCRY::LoadConfig() {
    UpdateParamsFromYAML<TRestGeant4SpatialGeneratorCRY>(fNode);
    if (ToVector(fVerticalAxis).R() == 0) throw runtime_error(GetClassName() + ": 'verticalAxis' is zero");
    UpdateYAMLFromParams<TRestGeant4SpatialGeneratorCRY>(fNode);
}

XYZVector TRestGeant4SpatialGeneratorCRY::Sample(TRandom&) {
    throw runtime_error(GetClassName() + ": the positions are the ones of the showers");
}

XYZVector TRestGeant4SpatialGeneratorCRY::GetCenter() const { return ToVector(fPosition); }

void TRestGeant4SpatialGeneratorCRY::SetEngineFactory(EngineFactory factory) {
    CRYEngineFactory() = std::move(factory);
}

bool TRestGeant4SpatialGeneratorCRY::HasEngine() { return static_cast<bool>(CRYEngineFactory()); }

void TRestGeant4SpatialGeneratorCRY::SetCRY(const TRestGeant4CRY& cry) {
    fCRYSetup = cry.GetCRYSetup();
    fSubBoxLength = cry.fSubBoxLength;
    fMaxExposureTime = cry.fMaxExposureTime;
}

void TRestGeant4SpatialGeneratorCRY::Start(long seed) {
    if (fState) return;
    if (!HasCRY()) throw runtime_error(GetClassName() + ": no TRestGeant4CRY section (the CRY options)");
    if (!CRYEngineFactory()) {
        throw runtime_error(GetClassName() +
                            ": restG4 was built without CRY (REST_CRY, CRY_PATH or $CRYHOME)");
    }
    auto state = make_shared<State>();
    state->setup = fCRYSetup;
    state->seed = seed;
    fState = state;
}

std::unique_ptr<TRestGeant4SpatialGenerator> TRestGeant4SpatialGeneratorCRY::Clone() const {
    auto copy = make_unique<TRestGeant4SpatialGeneratorCRY>(*this);
    copy->fEngine.reset();  // its own CRY, created at its first shower
    return copy;
}

bool TRestGeant4SpatialGeneratorCRY::NextShower(vector<TRestGeant4ParticleState>& particles) {
    if (!fState) throw runtime_error(GetClassName() + ": CRY is not started");
    if (!fEngine) {
        // The CRY of this copy, created in its thread (the random generator of CRY is one per thread)
        lock_guard<mutex> lock(fState->mutex);
        fEngine = make_shared<Engine>();
        const long seed = fState->seed + static_cast<long>(fState->engines.size());
        fEngine->engine = CRYEngineFactory()(fState->setup, seed);
        fState->engines.push_back(fEngine);
    }
    particles.clear();
    const XYZVector center = GetCenter();
    const bool world = fWorldHalfSize.R() > 0;
    // The showers without a particle crossing the world are skipped (they count in the exposure)
    while (particles.empty()) {
        if (fMaxExposureTime > 0 && fState->exposure.load() >= fMaxExposureTime) return false;
        vector<TRestGeant4ParticleState> shower = fEngine->engine->Generate();
        fState->showers++;
        const double exposure = fEngine->engine->GetTimeSimulated();
        AtomicAdd(fState->exposure, exposure - fEngine->exposure);
        fEngine->exposure = exposure;
        for (auto& particle : shower) {
            particle.SetOrigin(RotateToWorld(particle.GetOrigin()) + center);
            particle.SetDirection(RotateToWorld(particle.GetDirection()));
            if (!world || MoveToBox(particle, fWorldHalfSize)) particles.push_back(particle);
        }
        // Without the world, a shower is used as it is (an empty one too: CRY returns at least a particle)
        if (!world) break;
    }
    return true;
}

XYZVector TRestGeant4SpatialGeneratorCRY::RotateToWorld(const XYZVector& v) const {
    const XYZVector up = ToVector(fVerticalAxis).Unit();
    const XYZVector z(0, 0, 1);
    const XYZVector axis = z.Cross(up);
    const double angle = acos(clamp(z.Dot(up), -1.0, 1.0));
    const ROOT::Math::AxisAngle rotation(axis.R() > 0 ? axis.Unit() : XYZVector(1, 0, 0), angle);
    return rotation * v;
}

vector<XYZVector> TRestGeant4SpatialGeneratorCRY::GetPlaneCorners() const {
    const double half = fSubBoxLength / 2;
    vector<XYZVector> corners;
    for (const double x : {-half, half}) {
        for (const double y : {-half, half}) corners.push_back(RotateToWorld({x, y, 0}) + GetCenter());
    }
    return corners;
}

void TRestGeant4SpatialGeneratorCRY::UpdateExposure() {
    if (!fState) return;
    SetExposure(fState->showers.load(), fState->exposure.load());
}

void TRestGeant4SpatialGeneratorCRY::SetExposure(int64_t showers, double exposureTime) {
    fShowersUsed = showers;
    fExposureTime = exposureTime;
    UpdateYAMLFromParams<TRestGeant4SpatialGeneratorCRY>(fNode);
}
