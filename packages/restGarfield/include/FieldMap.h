#ifndef RESTGARFIELD_FIELDMAP_H
#define RESTGARFIELD_FIELDMAP_H

/// \file FieldMap.h
/// \brief Electric field maps of restGarfield: the ASCII export of CST Studio Suite as a Garfield++
/// ComponentGrid.

#include <array>
#include <memory>
#include <string>

namespace Garfield {
class ComponentGrid;
}

/// \brief Loads the electric field exported by CST Studio Suite as ASCII (2D/3D field result, ASCIIExport
/// with file type "ascii", a regular grid: Mode "FixedWidth" or "FixedNumber") into a Garfield++
/// ComponentGrid, without potential.
///
/// The file has a header line with the columns and their units, e.g.
/// `x [mm]  y [mm]  z [mm]  Ex [V/m]  Ey [V/m]  Ez [V/m]`, then a line of dashes and one line per node
/// `x y z Ex Ey Ez` (CST 2026). The units of the header (any length and field units of REST_Units, e.g. mm,
/// m and V/m, V/cm) convert the values to the Garfield units (cm, V/cm); the grid is deduced from the nodes
/// (ComponentGrid "XYZ"). The coordinates are moved by `offset` (mm), the position in the detector frame of
/// the origin of the CST frame. The nodes with the field exactly 0 (Ex = Ey = Ez = 0) are out of the gas:
/// CST writes 0 out of the solved volume (e.g. the corners of the exported box around a cylindrical gas
/// volume) and inside the conductors. They are inactive nodes of the grid: in a cell with one of them,
/// ComponentGrid gives no medium (status -5) and a drift ends there, like at a wall. The Garfield medium of
/// the grid is not set (ComponentGrid::SetMedium).
///
/// \param fileName ASCII export of CST.
/// \param offset Position (mm) of the origin of the CST frame in the detector frame.
/// \return The component, owned by the caller.
/// \throws std::runtime_error if the file cannot be read, has no header with the units of the coordinates
/// and the field, a node line cannot be read, or ComponentGrid rejects the grid (e.g. not regular).
std::unique_ptr<Garfield::ComponentGrid> LoadCSTFieldMap(const std::string& fileName,
                                                         const std::array<double, 3>& offset = {0, 0, 0});

#endif
