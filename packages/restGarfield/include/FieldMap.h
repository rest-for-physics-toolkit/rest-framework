#ifndef RESTGARFIELD_FIELDMAP_H
#define RESTGARFIELD_FIELDMAP_H

/// \file FieldMap.h
/// \brief Electric field maps of restGarfield: the ASCII export of CST Studio Suite as a Garfield++
/// ComponentGrid, and several of them (e.g. a fine map of the amplification structure and a coarser one above
/// it) as one component, with periodicity.

#include <array>
#include <memory>
#include <string>
#include <vector>

#include "Garfield/Component.hh"

namespace Garfield {
class ComponentGrid;
class Medium;
}  // namespace Garfield

/// \brief Region of a node of a field map: in the gas, in a conductor or out of the solved volume (the field
/// exactly 0), or in a dielectric (dielectric map); out of every map
enum class CSTRegion : unsigned char { kGas = 0, kConductor = 1, kDielectric = 2, kOutside = 3 };

/// \brief Loads the electric field exported by CST Studio Suite as ASCII (2D/3D field result, ASCIIExport
/// with file type "ascii", a regular grid: Mode "FixedWidth" or "FixedNumber") into a Garfield++
/// ComponentGrid, without potential.
///
/// The file has a header line with the columns and their units, e.g.
/// `x [mm]  y [mm]  z [mm]  Ex [V/m]  Ey [V/m]  Ez [V/m]`, then a line of dashes and one line per node
/// `x y z Ex Ey Ez` (CST 2026). The units of the header (any length and field units of REST_Units, e.g. mm,
/// um, m and V/m, V/cm) convert the values to the Garfield units (cm, V/cm); the grid is deduced from the
/// nodes (ComponentGrid "XYZ"). The coordinates are moved by `offset` (mm), the position in the detector
/// frame of the origin of the CST frame. The nodes with the field exactly 0 (Ex = Ey = Ez = 0) are out of the
/// gas: CST writes 0 out of the solved volume (e.g. the corners of the exported box around a cylindrical gas
/// volume) and inside the conductors. They are inactive nodes of the grid: in a cell with one of them,
/// ComponentGrid gives no medium (status -5) and a drift or an avalanche ends there, like at a wall. The
/// Garfield medium of the grid is not set (ComponentGrid::SetMedium).
///
/// The dielectrics (e.g. the kapton of a microbulk Micromegas) have a field and are not distinguished from
/// the gas by it: `dielectricFile` is the export of the electric displacement D on nodes of the same grid
/// (header `x [um] ... Dx [C/m^2] ...`; all of them or a part, e.g. only the layers with dielectrics). A
/// node where |D| > 1.5 eps0 |E| (relative permittivity above 1.5; 1 in the gas) is inactive as well:
/// the electrons that reach a dielectric are lost (its charging up is not simulated).
///
/// \param fileName ASCII export of the electric field by CST.
/// \param offset Position (mm) of the origin of the CST frame in the detector frame.
/// \param dielectricFile ASCII export of the displacement field D on nodes of the same grid (empty: none).
/// \return The component, owned by the caller.
/// \throws std::runtime_error if a file cannot be read, has no header with the units of the coordinates and
/// the field (C/m^2 for D), a node line cannot be read, ComponentGrid rejects the grid (e.g. not regular), or
/// a node of the dielectric file is not a node of the field map.
std::unique_ptr<Garfield::ComponentGrid> LoadCSTFieldMap(const std::string& fileName,
                                                         const std::array<double, 3>& offset = {0, 0, 0},
                                                         const std::string& dielectricFile = "");

/// \brief Several field maps as one Garfield++ component, e.g. a fine map of the amplification structure
/// (mesh, holes, dielectric) and a coarser map of the drift field above it, which cannot be one regular
/// grid. At each point, the field and the medium of the first map (in the order of AddMap) whose mesh
/// contains it; outside every map, status -6 and no medium.
///
/// Periodicity (Component::EnablePeriodicityX/Y/Z or EnableMirrorPeriodicityX/Y/Z, e.g. along x and y for the
/// unit cell of a readout): the period is the extent of the box enclosing all the maps along that axis; a
/// point is moved into it before choosing the map, and with mirror periodicity the field component along
/// that axis changes sign in the mirrored cells (as ComponentGrid). The maps themselves are not periodic.
class CSTFieldMap : public Garfield::Component {
   public:
    CSTFieldMap();
    /// \brief Out of line: ComponentGrid is only declared here
    ~CSTFieldMap() override;

    /// \brief Adds a map (e.g. from LoadCSTFieldMap), after the maps added before (they take precedence where
    /// they overlap), with the medium of SetMedium.
    /// \throws std::runtime_error if the map has no mesh.
    void AddMap(std::unique_ptr<Garfield::ComponentGrid> map);
    /// \brief Sets the Garfield medium of all the maps (also those added later).
    void SetMedium(Garfield::Medium* medium);
    /// \brief Number of maps
    size_t GetNumberOfMaps() const { return fMaps.size(); }

    void ElectricField(const double x, const double y, const double z, double& ex, double& ey, double& ez,
                       Garfield::Medium*& m, int& status) override;
    void ElectricField(const double x, const double y, const double z, double& ex, double& ey, double& ez,
                       double& v, Garfield::Medium*& m, int& status) override;
    Garfield::Medium* GetMedium(const double x, const double y, const double z) override;
    /// \brief Where an electron that stopped at a point (cm) ended: the region of the nodes of the map around
    /// it (its cell and the cells next to it: Garfield++ ends an electron just past the boundary of the gas),
    /// dielectric if one of them is in a dielectric, else conductor if one of them has the field 0, else gas;
    /// kOutside out of every map.
    CSTRegion GetRegionNear(double x, double y, double z) const;
    /// \brief No potential: false.
    bool GetVoltageRange(double& vmin, double& vmax) override;
    /// \brief The box enclosing the maps (cm), infinite along the periodic axes.
    bool GetBoundingBox(double& xmin, double& ymin, double& zmin, double& xmax, double& ymax,
                        double& zmax) override;
    /// \brief The box enclosing the maps (cm): the period along the periodic axes.
    bool GetElementaryCell(double& xmin, double& ymin, double& zmin, double& xmax, double& ymax,
                           double& zmax) override;

   protected:
    void Reset() override;
    void UpdatePeriodicity() override;

   private:
    /// The map containing a point (cm) moved into the cell along the periodic axes: its index (-1: none),
    /// the point in the cell and whether it is mirrored along each axis
    int Locate(double x, double y, double z, std::array<double, 3>& cell,
               std::array<bool, 3>& mirrored) const;

    std::vector<std::unique_ptr<Garfield::ComponentGrid>> fMaps;
    /// Limits (cm) of each map: low and high along x, y, z
    std::vector<std::array<double, 6>> fLimits;
    /// Box enclosing the maps (cm)
    std::array<double, 3> fLow = {0, 0, 0}, fHigh = {0, 0, 0};
    Garfield::Medium* fMedium = nullptr;
};

/// \brief One map of LoadCSTFieldMaps: the export of the electric field and, optionally, that of the
/// displacement field D on the same nodes (the dielectrics, see LoadCSTFieldMap).
struct CSTFieldMapFiles {
    std::string fField;
    std::string fDielectric;
};

/// \brief Loads several CST exports (LoadCSTFieldMap, all with the same offset) into one CSTFieldMap, in the
/// order given (a map takes precedence over the following ones where they overlap: the fine map first).
/// \throws std::runtime_error if there are no maps or one cannot be loaded.
std::unique_ptr<CSTFieldMap> LoadCSTFieldMaps(const std::vector<CSTFieldMapFiles>& maps,
                                              const std::array<double, 3>& offset = {0, 0, 0});

#endif
