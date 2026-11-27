#include "TRestGeant4PrimaryGeneratorInfo.h"

#include <stdexcept>

using namespace std;

static const bool TRestGeant4PrimaryGeneratorInfo_FieldsRegistered = []() {
    auto& reg = TRestMetadataFieldRegistry::Instance();
    reg.RegisterField<TRestGeant4PrimaryGeneratorInfo>(
        "spatialDensityFunction", &TRestGeant4PrimaryGeneratorInfo::fSpatialDensityFunction);
    reg.RegisterField<TRestGeant4PrimaryGeneratorInfo>("source", &TRestGeant4PrimaryGeneratorInfo::fSources);
    reg.RegisterField<TRestGeant4PrimaryGeneratorInfo>("activity",
                                                       &TRestGeant4PrimaryGeneratorInfo::fActivity);
    return true;
}();

TRestGeant4PrimaryGeneratorInfo::TRestGeant4PrimaryGeneratorInfo(const string& name, const YAML::Node& node)
    : TRestMetadata(name, node) {
    LoadConfig();
}

TRestGeant4PrimaryGeneratorInfo::TRestGeant4PrimaryGeneratorInfo(const TRestGeant4PrimaryGeneratorInfo& other)
    : TRestMetadata(other) {
    CopyGenerator(other);
}

TRestGeant4PrimaryGeneratorInfo& TRestGeant4PrimaryGeneratorInfo::operator=(
    const TRestGeant4PrimaryGeneratorInfo& other) {
    if (this != &other) {
        TRestMetadata::operator=(other);
        CopyGenerator(other);
    }
    return *this;
}

void TRestGeant4PrimaryGeneratorInfo::CopyGenerator(const TRestGeant4PrimaryGeneratorInfo& other) {
    // Its own section: writing the copy must not modify the section of the original (nor its sources)
    fNode = YAML::Clone(other.fNode);
    fSpatialDensityFunction = other.fSpatialDensityFunction;
    fSources.clear();
    for (const auto& source : other.fSources) fSources.push_back(source.Clone());
    // Not rebuilt from the section: SetGdmlPlacement and SetCosmicRadius modify it
    fSpatialGenerator = other.fSpatialGenerator ? other.fSpatialGenerator->Clone() : nullptr;
    // It writes to the section of the copy (e.g. the cry exposure), not to the one of the original
    if (fSpatialGenerator) fSpatialGenerator->SetYAMLNode(fNode);
}

void TRestGeant4PrimaryGeneratorInfo::LoadConfig() {
    // The sources (`source`, one or a list) are a registered field
    fSources.clear();
    UpdateParamsFromYAML<TRestGeant4PrimaryGeneratorInfo>(fNode);

    // The spatial generator reads its fields (type, from, shape...) from the generator section itself
    fSpatialGenerator = TRestGeant4SpatialGenerator::Create("generator", fNode);

    if (IsCRY()) {
        if (!fSources.empty()) {
            throw runtime_error(GetClassName() + ": the cry generator has no 'source' (CRY gives them)");
        }
    } else if (fSources.empty()) {
        throw runtime_error(GetClassName() + ": no 'source' in the generator");
    }

    if (!fSpatialDensityFunction.empty() &&
        !dynamic_cast<const TRestGeant4SpatialGeneratorShape*>(fSpatialGenerator.get())) {
        throw runtime_error(GetClassName() +
                            ": 'spatialDensityFunction' needs a volume or surface generator");
    }
    if (fActivity < 0) throw runtime_error(GetClassName() + ": 'activity' must be positive");
    if (fActivity > 0 && (IsCosmic() || IsCRY())) {
        throw runtime_error(GetClassName() +
                            ": 'activity' is not used by the cosmic and cry generators (their exposure is "
                            "computed from their flux)");
    }
    if (fActivity > 0 && GetSurfaceRate() > 0) {
        throw runtime_error(GetClassName() +
                            ": 'activity' and 'flux' given, only one of them can set the rate");
    }
    if (IsCosmic() && (fSources.size() != 1 || fSources[0].IsDecay0())) {
        throw runtime_error(GetClassName() +
                            ": the cosmic generator needs a single source, not a Decay0 one");
    }

    UpdateYAMLFromParams<TRestGeant4PrimaryGeneratorInfo>(fNode);
}

std::set<std::string> TRestGeant4PrimaryGeneratorInfo::GetSharedNodeKeys() const {
    return TRestMetadataFieldRegistry::Instance().GetFieldKeys(fSpatialGenerator.get());
}

void TRestGeant4PrimaryGeneratorInfo::LoadDistributions(long seed) {
    for (auto& source : fSources) source.LoadDistributions();
    // The cosmic generator (a single source): the flux of its source, to normalize to a time
    if (auto* cosmic = dynamic_cast<TRestGeant4SpatialGeneratorCosmic*>(fSpatialGenerator.get())) {
        if (fSources.size() == 1) cosmic->SetFlux(fSources[0].GetFlux());
    }
    if (auto* cry = dynamic_cast<TRestGeant4SpatialGeneratorCRY*>(fSpatialGenerator.get())) cry->Start(seed);
}

bool TRestGeant4PrimaryGeneratorInfo::HasDecay0Sources() const {
    for (const auto& source : fSources) {
        if (source.IsDecay0()) return true;
    }
    return false;
}

bool TRestGeant4PrimaryGeneratorInfo::IsCosmic() const {
    return dynamic_cast<const TRestGeant4SpatialGeneratorCosmic*>(fSpatialGenerator.get()) != nullptr;
}

bool TRestGeant4PrimaryGeneratorInfo::IsCRY() const {
    return dynamic_cast<const TRestGeant4SpatialGeneratorCRY*>(fSpatialGenerator.get()) != nullptr;
}

void TRestGeant4PrimaryGeneratorInfo::SetCRY(const TRestGeant4CRY& cry) {
    if (auto* generator = dynamic_cast<TRestGeant4SpatialGeneratorCRY*>(fSpatialGenerator.get())) {
        generator->SetCRY(cry);
    }
}

vector<ROOT::Math::XYZVector> TRestGeant4PrimaryGeneratorInfo::GetCRYPlaneCorners() const {
    const auto* cry = dynamic_cast<const TRestGeant4SpatialGeneratorCRY*>(fSpatialGenerator.get());
    return cry ? cry->GetPlaneCorners() : vector<ROOT::Math::XYZVector>();
}

void TRestGeant4PrimaryGeneratorInfo::SetCRYWorld(const ROOT::Math::XYZVector& halfSize) {
    if (auto* cry = dynamic_cast<TRestGeant4SpatialGeneratorCRY*>(fSpatialGenerator.get())) {
        cry->SetWorld(halfSize);
    }
}

void TRestGeant4PrimaryGeneratorInfo::UpdateCRYExposure() {
    if (auto* cry = dynamic_cast<TRestGeant4SpatialGeneratorCRY*>(fSpatialGenerator.get())) {
        cry->UpdateExposure();
    }
}

int64_t TRestGeant4PrimaryGeneratorInfo::GetCRYShowersUsed() const {
    const auto* cry = dynamic_cast<const TRestGeant4SpatialGeneratorCRY*>(fSpatialGenerator.get());
    return cry ? cry->fShowersUsed : 0;
}

double TRestGeant4PrimaryGeneratorInfo::GetCRYExposureTime() const {
    const auto* cry = dynamic_cast<const TRestGeant4SpatialGeneratorCRY*>(fSpatialGenerator.get());
    return cry ? double(cry->fExposureTime) : 0;
}

void TRestGeant4PrimaryGeneratorInfo::SetCRYExposure(int64_t showers, double exposureTime) {
    if (auto* cry = dynamic_cast<TRestGeant4SpatialGeneratorCRY*>(fSpatialGenerator.get())) {
        cry->SetExposure(showers, exposureTime);
    }
}

bool TRestGeant4PrimaryGeneratorInfo::IsGdml() const {
    return fSpatialGenerator && fSpatialGenerator->IsGdml();
}

bool TRestGeant4PrimaryGeneratorInfo::IsGdmlSurface() const {
    return IsGdml() && dynamic_cast<const TRestGeant4SpatialGeneratorSurface*>(fSpatialGenerator.get());
}

string TRestGeant4PrimaryGeneratorInfo::GetGdmlVolume() const {
    const auto* shape = dynamic_cast<const TRestGeant4SpatialGeneratorShape*>(fSpatialGenerator.get());
    return shape ? shape->fFrom : "";
}

void TRestGeant4PrimaryGeneratorInfo::SetGdmlPlacement(const ROOT::Math::XYZVector& center,
                                                       const ROOT::Math::XYZVector& axis, double angle) {
    if (auto* shape = dynamic_cast<TRestGeant4SpatialGeneratorShape*>(fSpatialGenerator.get())) {
        shape->SetGdmlPlacement(center, axis, angle);
    }
}

void TRestGeant4PrimaryGeneratorInfo::SetCosmicRadius(double radius) {
    if (auto* cosmic = dynamic_cast<TRestGeant4SpatialGeneratorCosmic*>(fSpatialGenerator.get())) {
        cosmic->SetDefaultRadius(radius);
    }
}

double TRestGeant4PrimaryGeneratorInfo::GetCosmicSurfaceCm2() const {
    const auto* cosmic = dynamic_cast<const TRestGeant4SpatialGeneratorCosmic*>(fSpatialGenerator.get());
    return cosmic ? cosmic->GetSurfaceCm2() : 0;
}

double TRestGeant4PrimaryGeneratorInfo::GetCosmicFlux() const {
    const auto* cosmic = dynamic_cast<const TRestGeant4SpatialGeneratorCosmic*>(fSpatialGenerator.get());
    return cosmic ? cosmic->fFlux : 0;
}

double TRestGeant4PrimaryGeneratorInfo::GetSurfaceRate() const {
    const auto* surface = dynamic_cast<const TRestGeant4SpatialGeneratorSurface*>(fSpatialGenerator.get());
    return surface ? surface->fFlux * surface->GetFluxSurfaceCm2() : 0;
}

void TRestGeant4PrimaryGeneratorInfo::Print() const {
    if (!fSpatialGenerator) return;
    RESTMetadata << "Generator: " << fSpatialGenerator->GetClassName() << RESTendl;
    if (IsGdml()) RESTMetadata << "Generator volume: " << GetGdmlVolume() << RESTendl;
    if (IsCosmic()) {
        RESTMetadata << "Cosmic generation surface: " << GetCosmicSurfaceCm2() << " cm2, flux "
                     << GetCosmicFlux() << " counts/s/cm2" << RESTendl;
    }
    if (GetSurfaceRate() > 0) {
        RESTMetadata << "Surface flux: " << GetSurfaceRate() << " primaries/s" << RESTendl;
    }
    if (IsCRY()) {
        const auto* cry = dynamic_cast<const TRestGeant4SpatialGeneratorCRY*>(fSpatialGenerator.get());
        RESTMetadata << "CRY: " << cry->GetCRYSetup() << RESTendl;
    }
    if (!fSpatialDensityFunction.empty()) {
        RESTMetadata << "Density: " << fSpatialDensityFunction << RESTendl;
    }
    for (const auto& source : fSources) {
        if (source.IsDecay0()) {
            RESTMetadata << "Source: Decay0 file " << source.fUse << RESTendl;
        } else {
            RESTMetadata << "Source: " << source.fParticleName << ", "
                         << source.GetEnergyDistribution()->GetClassName() << ", "
                         << source.GetAngularDistribution()->GetClassName() << RESTendl;
        }
    }
}
