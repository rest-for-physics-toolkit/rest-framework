#ifndef TRESTFIBREREADOUT_H
#define TRESTFIBREREADOUT_H

#include <array>
#include <string>
#include <vector>

#include "TRestDetectorReadout.h"

/// \class TRestFibreFamily
/// \brief A family of parallel fibres of a TRestFibreReadout (an entry of `fibres`): their axis and the grid
/// of their positions in the plane transverse to it.
class TRestFibreFamily : public TRestMetadata {
   public:
    /// `axis`: direction of the fibres, "X", "Y" or "Z"
    std::string fAxis = "Z";
    /// `pitch`: distance between neighbour fibres, the same along both transverse axes
    TRestWithUnits fPitch = 0;
    /// `offset`: position of the first fibre along the two transverse axes (in cyclic order after the axis:
    /// X: Y and Z; Y: Z and X; Z: X and Y), from the lower corner of the box of the readout
    std::array<TRestWithUnits, 2> fOffset = {0.0, 0.0};
    /// `volume`: logical volume of the fibre in the GDML of the simulation (the placements written by
    /// TRestFibreReadout::WriteGDMLPlacements)
    std::string fVolume = "";
    /// `readoutEnd`: end of the fibres read by a photosensor, "+" (the end at the upper side of the box along
    /// the axis), "-" (lower side), "both" (a channel per end), or alternating along alternateAxis, "+-"
    /// (the first fibre read at +, the next one at -, ...) or "-+"
    std::string fReadoutEnd = "+";
    /// `alternateAxis`: with readoutEnd "+-" or "-+", the transverse axis ("X", "Y" or "Z") along which
    /// consecutive fibres alternate their read end (whole rows or columns read at the same end)
    std::string fAlternateAxis = "";

    /// \brief Builds an empty family.
    TRestFibreFamily();
    /// \brief Builds the family from its YAML node (LoadConfig).
    TRestFibreFamily(const std::string& name, const YAML::Node& node);
    /// \brief Reads the keys; throws if the axis is not X, Y or Z, the pitch is not positive, readoutEnd is
    /// not +, -, both, +- or -+, or alternateAxis is not a transverse axis with +- and -+ (empty otherwise).
    void LoadConfig() override;
    /// \brief Nothing to initialise.
    void Initialize() override {}
    /// \brief Returns the class name, "TRestFibreFamily".
    std::string GetClassName() const override { return "TRestFibreFamily"; }

    /// \brief Axis of the fibres: 0 X, 1 Y, 2 Z.
    int GetAxis() const { return fAxis == "X" ? 0 : fAxis == "Y" ? 1 : 2; }
    /// \brief Read ends: {-1}, {+1} or {-1, +1}; alternating, the end of the first fibre ({+1} for "+-").
    std::vector<int> GetReadoutEnds() const {
        if (fReadoutEnd == "both") return {-1, 1};
        return {fReadoutEnd == "-" || fReadoutEnd == "-+" ? -1 : 1};
    }
    /// \brief Axis along which the read end alternates: 0 X, 1 Y, 2 Z, -1 if it does not alternate.
    int GetAlternateAxis() const {
        if (fReadoutEnd != "+-" && fReadoutEnd != "-+") return -1;
        return fAlternateAxis == "X" ? 0 : fAlternateAxis == "Y" ? 1 : 2;
    }
};

/// \class TRestFibreReadoutParameters
/// \brief The `readoutParameters` of a TRestFibreReadout: the box (size, position) and the fibre families
/// (`fibres`).
class TRestFibreReadoutParameters : public TRestMetadata {
   public:
    /// `size`: size of the box filled with fibres (the scintillator), along X, Y and Z
    std::array<TRestWithUnits, 3> fSize = {0.0, 0.0, 0.0};
    /// `position`: centre of the box in the global frame (the placement of the scintillator volume, not
    /// rotated)
    std::array<TRestWithUnits, 3> fPosition = {0.0, 0.0, 0.0};
    /// `fibres`: the fibre families
    std::vector<TRestFibreFamily> fFibres;

    /// \brief Builds empty parameters.
    TRestFibreReadoutParameters();
    /// \brief Builds the parameters from the `readoutParameters` node (LoadConfig).
    TRestFibreReadoutParameters(const std::string& name, const YAML::Node& node);
    /// \brief Reads the parameters; throws if the size is not positive or there are no fibres.
    void LoadConfig() override;
    /// \brief Nothing to initialise.
    void Initialize() override {}
    /// \brief Returns the class name, "TRestFibreReadoutParameters".
    std::string GetClassName() const override { return "TRestFibreReadoutParameters"; }
};

/// \class TRestFibreReadout
/// \brief Readout of a scintillator crossed by a lattice of fibres, each one read by a photosensor (e.g.
/// LiquidO, an opaque scintillator with wavelength shifting fibres read by SiPMs).
///
/// The parameters and an example are described in the class documentation (TRestFibreReadout.cxx).
///
/// \author JuanAn García <juanangp@unizar.es>
class TRestFibreReadout : public TRestDetectorReadout {
   public:
    TRestFibreReadoutParameters fReadoutParams;  ///< `readoutParameters`

   private:
    /// The grid of a family: axis, transverse axes, number of fibres along them, read ends, first physical
    /// ID (a block of nU * nV channels per read end) and transverse index (0: i along u, 1: j along v)
    /// along which the read end alternates (-1: it does not; the fibres with an odd index are read at the
    /// end opposite to ends[0])
    struct Grid {
        int axis = 2;
        int u = 0;
        int v = 1;
        int nU = 0;
        int nV = 0;
        std::vector<int> ends = {1};
        int firstID = 0;
        int alternate = -1;
    };
    std::vector<Grid> fGrids;  //!

    /// \brief Builds the grids of the families from the parameters (LoadConfig).
    void BuildGrids();
    /// \brief Lower corner of the box.
    double GetMinimum(int axis) const {
        return fReadoutParams.fPosition[axis] - fReadoutParams.fSize[axis] / 2;
    }
    /// \brief Transverse coordinate of fibre index along the transverse axis t (0: u, 1: v) of a family.
    double GetFibreCoordinate(size_t family, int t, int index) const;
    /// \brief Family, indices and read end (-1 or +1) of a physical ID; false if it is not a channel.
    bool GetFibre(int physicalID, size_t& family, int& i, int& j, int& end) const;
    /// \brief Coordinate along the axis of a family of its end (-1: lower side of the box, +1: upper).
    double GetEndCoordinate(size_t family, int end) const;
    /// \brief Read end (-1 or +1) of the fibre i, j of a grid for its read end index e.
    static int GetFibreEnd(const Grid& grid, size_t e, int i, int j);

   public:
    /// \brief Constructor.
    TRestFibreReadout();
    /// \brief Builds a fibre readout from an in-memory YAML node section.
    TRestFibreReadout(const std::string& instanceName, const YAML::Node& node);
    /// \brief Builds a fibre readout from file+section metadata references.
    TRestFibreReadout(const std::string& fileName, const std::string& sectionName);

    /// \brief Returns the class name, "TRestFibreReadout".
    std::string GetClassName() const override { return "TRestFibreReadout"; }
    /// \brief Loads the parameters and builds the grids of the families.
    void LoadConfig() override;

    /// \brief Builds the geometry (the fibres as thin tubes, to draw them) and the decoding: the decoding
    /// file if `decodingFile` is given, otherwise the DAQ ID of each fibre is its physical ID.
    void BuildGeometry() override;

    /// \brief Total number of fibres.
    int GetNumberOfFibres() const;
    /// \brief Total number of channels (read ends of the fibres).
    int GetNumberOfChannels() const;
    /// \brief The closest fibre (any family) to a position inside the box, measured in the plane transverse
    /// to the fibre, and of a fibre read at both ends the channel of the end closer to the position; -1
    /// outside the box or if it is not connected.
    int GetChannelFromPosition(double x, double y, double z) const override;
    /// \brief Position of a fibre: its two transverse coordinates, the coordinate along it NaN (e.g. a
    /// fibre along Y gives an XZ hit).
    ROOT::Math::XYZVector GetPositionFromChannel(int daqID) const override;
    /// \brief The box of the readout.
    std::pair<ROOT::Math::XYZVector, ROOT::Math::XYZVector> GetChannelsBoundingBox() const override;
    /// \brief Neighbour fibres of the same family and row (physical IDs - 1 and + 1 along the second
    /// transverse axis), -1 at the ends.
    std::pair<int, int> GetAdjacentChannels(int daqID) const override;
    /// \brief The channels of the fibres closer than maxDistance to the position in the plane transverse to
    /// them (the coordinate along the fibre inside the box), one per read end, with that distance, at the
    /// position of their read end (GetReadoutPosition).
    std::vector<ChannelInRange> GetChannelsInRange(const ROOT::Math::XYZVector& position,
                                                   double maxDistance) const override;
    /// \brief Position of the read end of the fibre of a channel (on the face of the box); NaN if it is not
    /// a channel.
    ROOT::Math::XYZVector GetReadoutPosition(int daqID) const override;
    /// \brief Distance along the fibre of a channel from a position to its read end (the light travels it
    /// inside the fibre, e.g. TRestDetectorLightCollectionProcess with a wls medium); NaN if it is not a
    /// channel.
    double GetDistanceToReadout(int daqID, const ROOT::Math::XYZVector& position) const override;
    /// The distance without a channel is that of the base class (NaN: the fibres have no drift)
    using TRestDetectorReadout::GetDistanceToReadout;
    /// \brief Paths of the light guided by the fibre of a channel from a position to its read end: the
    /// direct one (GetDistanceToReadout) and, if the fibre is read at one end only, the one reflected at
    /// the other end (2 L - d, L the length of the fibre); empty if it is not a channel.
    std::vector<GuidedLightPath> GetGuidedLightPaths(int daqID,
                                                     const ROOT::Math::XYZVector& position) const override;

    /// \brief Writes the placements of the fibres (`physvol` elements of the GDML, named fibre_<physical
    /// ID>, with the logical volume `volume` of their family) to a file, to be included as an external
    /// entity inside the volume of the scintillator (placed at `position`, not rotated).
    /// \return false if the file cannot be written or a family has no volume.
    bool WriteGDMLPlacements(const std::string& fileName) const;
};

#endif
