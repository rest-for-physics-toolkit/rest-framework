#ifndef TRESTDETECTORDRIFTRESPONSE_H
#define TRESTDETECTORDRIFTRESPONSE_H

#include <array>
#include <string>
#include <utility>
#include <vector>

#include "Math/Vector3D.h"
#include "TRestMetadata.h"

class TFile;
class TRandom;

/// \class TRestDetectorDriftResponse
/// \brief Drift of the electrons in a non-uniform field (e.g. a CST field map) as a response per start
/// point: the configuration of its computation (the `drift` key of a TRestDetectorTransportManager, computed
/// by restGarfield with Garfield++) and the response, written to and read from its own file without
/// Garfield++, used by TRestDetectorHitsDriftResponseProcess.
///
/// The parameters (YAML keys), the file format and the use of the response are described in the class
/// documentation (TRestDetectorDriftResponse.cxx).
///
/// \author JuanAn García <juanangp@unizar.es>
class TRestDetectorDriftResponse : public TRestMetadata {
   public:
    /// Name of the response in RESTMetadataStore and of its tree in the response file
    inline static const std::string kName = "DriftResponse";

    // Configuration (`drift` key of TRestDetectorTransportManager)
    /// `fieldMap`: ASCII export of the field by CST Studio Suite (restGarfield LoadCSTFieldMap)
    std::string fFieldMap = "";
    /// `fieldMapOffset`: position (mm) of the origin of the CST frame in the detector frame
    std::vector<TRestWithUnits> fFieldMapOffset = {0, 0, 0};
    /// `readouts`: YAML file with the TRestDetectorReadoutManager of the readout
    std::string fReadoutsFile = "";
    /// `readout`: name of the readout (empty: the data readout of the manager)
    std::string fReadoutName = "";
    /// `boundary`: distance (mm) from the readout to the plane where the drift ends (above the mesh, where
    /// the amplification response starts)
    TRestWithUnits fBoundary = 0;
    /// `transportTable`: key of the transport table of the medium file (empty: the only one)
    std::string fTransportTable = "";
    /// `xRange`, `yRange`, `zRange`: grid of start points (mm, detector frame)
    std::pair<TRestWithUnits, TRestWithUnits> fXRange = {0, 0};
    /// `yRange` of the grid (mm)
    std::pair<TRestWithUnits, TRestWithUnits> fYRange = {0, 0};
    /// `zRange` of the grid (mm)
    std::pair<TRestWithUnits, TRestWithUnits> fZRange = {0, 0};
    /// `gridNodes`: nodes of the grid along x, y and z (at least 2 each)
    std::vector<int> fGridNodes = {2, 2, 2};
    /// `electrons`: primary electrons drifted from each node
    int fElectrons = 1000;
    /// `stepSize`: distance step of the drift (mm, Garfield++ AvalancheMC)
    TRestWithUnits fStepSize = 1;
    /// `positionBins`: bins per axis of the 2D histogram of the arrival position of each node
    int fPositionBins = 20;
    /// `timeBins`: bins of the histogram of the drift time of each node
    int fTimeBins = 20;
    /// `responseFile`: file of the response
    std::string fResponseFile = "";
    /// `overwrite`: compute the response again if the file exists
    bool fOverwrite = false;

    /// \brief Response of a node of the grid
    struct Node {
        int fPrimaries = 0;  ///< Primary electrons drifted
        int fArrived = 0;    ///< Electrons that reached the boundary plane
        /// Range of the displacement (mm) along the axes u and v of the readout plane (PlaneAxes)
        double fUMin = 0, fUMax = 0, fVMin = 0, fVMax = 0;
        /// Range of the drift time (us)
        double fTMin = 0, fTMax = 0;
        /// Counts of the displacement, positionBins x positionBins (u bin + positionBins * v bin)
        std::vector<int> fPositionCounts;
        /// Counts of the drift time, timeBins
        std::vector<int> fTimeCounts;
    };

    // Response (filled by restGarfield or Read)
    /// Medium, its pressure (bar) and temperature (K), and the transport table used
    std::string fMedium = "";
    double fPressure = 0;     ///< Pressure (bar)
    double fTemperature = 0;  ///< Temperature (K)
    /// Nodes of the grid, x fastest (NodeIndex)
    std::vector<Node> fNodes;  //!

    /// \brief Builds an empty response.
    TRestDetectorDriftResponse();
    /// \brief Builds a configuration from a `drift` node (LoadConfig).
    TRestDetectorDriftResponse(const std::string& instanceName, const YAML::Node& node);

    /// \brief Returns the class name, "TRestDetectorDriftResponse".
    std::string GetClassName() const override { return "TRestDetectorDriftResponse"; }
    /// \brief Reads the configuration; throws if a value is out of range (fieldMapOffset not 3 values, a
    /// range not min < max, gridNodes not 3 values >= 2, electrons, stepSize, bins not > 0, boundary < 0).
    void LoadConfig() override;
    /// \brief Nothing to initialise.
    void Initialize() override {}

    /// \brief Number of nodes of the grid
    size_t GetNumberOfNodes() const;
    /// \brief Index of the node (i, j, k) of the grid, x fastest
    size_t NodeIndex(int i, int j, int k) const;
    /// \brief Position (mm) of a node of the grid
    ROOT::Math::XYZVector GetNodePosition(size_t index) const;
    /// \brief Whether a position (mm) is inside the grid
    bool Contains(const ROOT::Math::XYZVector& position) const;

    /// \brief Axes u and v of the readout plane of a normal: u is the x axis (the y axis for a normal close
    /// to x) without its component along the normal, v = normal x u. E.g. normal +z: u = x, v = y.
    static std::pair<ROOT::Math::XYZVector, ROOT::Math::XYZVector> PlaneAxes(
        const ROOT::Math::XYZVector& normal);

    /// \brief Sets the response of a node from the arrivals of its electrons (displacement u, v in mm and
    /// drift time in us of each one): the histograms over the range of the arrivals.
    void SetNode(size_t index, int primaries, const std::vector<std::array<double, 3>>& arrivals);

    /// \brief Fraction of the electrons from a position (mm) that reach the boundary plane: the fractions
    /// arrived / primaries of the 8 neighbouring nodes interpolated linearly. 0 outside the grid.
    double GetSurvival(const ROOT::Math::XYZVector& position) const;
    /// \brief Arrival of an electron from a position (mm) that reaches the boundary: one of the 8
    /// neighbouring nodes chosen with its interpolation weight times its survival, and its displacement (u,
    /// v, mm) and drift time (us) sampled from its histograms (uniform within the bin).
    /// \return false outside the grid or if no neighbouring node has arrivals.
    bool SampleArrival(const ROOT::Math::XYZVector& position, TRandom& random, double& u, double& v,
                       double& time) const;

    /// \brief Writes the response to a file (opened with RECREATE or UPDATE): its metadata (configuration,
    /// medium, conditions, grid, bins and generator) in RESTMetadataStore and a tree with a node per entry.
    void Write(TFile* file, const YAML::Node& generator) const;
    /// \brief Reads the response of a file. Throws if it has none or it is not consistent.
    static TRestDetectorDriftResponse Read(TFile* file);

   private:
    /// The 8 neighbouring nodes of a position and their interpolation weights; false outside the grid
    bool Neighbours(const ROOT::Math::XYZVector& position, std::array<size_t, 8>& nodes,
                    std::array<double, 8>& weights) const;
};

#endif
