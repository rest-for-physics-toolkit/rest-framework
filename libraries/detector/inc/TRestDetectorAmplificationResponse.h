#ifndef TRESTDETECTORAMPLIFICATIONRESPONSE_H
#define TRESTDETECTORAMPLIFICATIONRESPONSE_H

#include <array>
#include <string>
#include <vector>

#include "TRestMetadata.h"

class TFile;
class TRandom;

/// \class TRestDetectorAmplificationResponse
/// \brief Amplification of the electrons in the field map of a unit cell of the amplification structure
/// (e.g. a microbulk Micromegas exported by CST): the configuration of its computation (the `amplification`
/// key of a TRestDetectorTransportManager, computed by restGarfield with Garfield++ microscopic avalanches)
/// and the response (transparency, avalanche size and time), written to and read from its own file without
/// Garfield++, used by TRestDetectorAvalancheProcess (`responseFile`).
///
/// The parameters (YAML keys), the file format and the use of the response are described in the class
/// documentation (TRestDetectorAmplificationResponse.cxx).
///
/// \author JuanAn García <juanangp@unizar.es>
class TRestDetectorAmplificationResponse : public TRestMetadata {
   public:
    /// Name of the response in RESTMetadataStore
    inline static const std::string kName = "AmplificationResponse";

    // Configuration (`amplification` key of TRestDetectorTransportManager)
    /// `fieldMaps`: ASCII exports of the electric field by CST Studio Suite, the first taking precedence
    /// where they overlap (e.g. the fine map of the mesh and the coarse map above it)
    std::vector<std::string> fFieldMaps = {};
    /// `dielectricMaps`: exports of the displacement field D on the nodes of each field map, in the same
    /// order ("" or missing: none); their dielectric nodes stop the electrons
    std::vector<std::string> fDielectricMaps = {};
    /// `fieldMapOffset`: position (mm) of the origin of the CST frame in the frame of the response
    std::vector<TRestWithUnits> fFieldMapOffset = {0, 0, 0};
    /// `periodicity`: axes along which the cell repeats ("", "x", "y", "xy"...): the period is the extent of
    /// the maps along them
    std::string fPeriodicity = "xy";
    /// `mirror`: mirror periodicity instead of a simple translation
    bool fMirror = false;
    /// `start`: z (mm) of the plane where the electrons start, the boundary plane of the drift response
    TRestWithUnits fStart = 0;
    /// `anode`: z (mm) below which an electron is collected by the anode (the electrons drift towards -z)
    TRestWithUnits fAnode = 0;
    /// `electrons`: primary electrons, started at random over the cell on the start plane
    int fElectrons = 1000;
    /// `usePenning`: enable the Penning transfer of the gas (penningRate)
    bool fUsePenning = false;
    /// `penningRate`: Penning transfer probability
    double fPenningRate = 0;
    /// `sizeBins`: bins of the histogram of the avalanche size
    int fSizeBins = 100;
    /// `timeBins`: bins of the histogram of the time from the start plane to the anode
    int fTimeBins = 50;
    /// `responseFile`: file of the response
    std::string fResponseFile = "";
    /// `overwrite`: compute the response again if the file exists
    bool fOverwrite = false;

    // Response (filled by restGarfield or Read)
    /// Medium, its pressure (bar) and temperature (K), and the Penning rate used (0: none)
    std::string fMedium = "";
    double fPressure = 0;     ///< Pressure (bar)
    double fTemperature = 0;  ///< Temperature (K)
    double fPenning = 0;      ///< Penning rate used
    /// Primary electrons started and those whose avalanche reached the anode (at least one electron)
    int fPrimaries = 0;
    /// Primary electrons whose avalanche reached the anode
    int fCollected = 0;

    /// \brief Where the primary electrons whose avalanche did not reach the anode ended (its first electron
    /// that stopped)
    struct Losses {
        int fOut = 0;         ///< Out of the field maps or of the cell (e.g. back above the start plane)
        int fConductor = 0;   ///< On a conductor (nodes with the field 0, e.g. the mesh)
        int fDielectric = 0;  ///< On a dielectric (e.g. the kapton)
        int fAttached = 0;    ///< Attached in the gas
        int fOther = 0;       ///< Elsewhere (e.g. in the gas, below the transport cut)
    };
    /// Losses of the primary electrons (fPrimaries - fCollected in total)
    Losses fLosses;
    /// Mean and variance of the size of the collected avalanches (electrons at the anode)
    double fMeanSize = 0;
    /// Variance of the size of the collected avalanches
    double fSizeVariance = 0;
    /// Histogram of the size of the collected avalanches over [0, fSizeMax]
    double fSizeMax = 0;
    /// Counts of the avalanche size, sizeBins
    std::vector<int> fSizeCounts;  //!
    /// Range (us) of the time from the start plane to the anode (mean time of the electrons of an avalanche)
    double fTMin = 0, fTMax = 0;
    /// Counts of the time, timeBins
    std::vector<int> fTimeCounts;  //!

    /// \brief Builds an empty response.
    TRestDetectorAmplificationResponse();
    /// \brief Builds a configuration from an `amplification` node (LoadConfig).
    TRestDetectorAmplificationResponse(const std::string& instanceName, const YAML::Node& node);

    /// \brief Returns the class name, "TRestDetectorAmplificationResponse".
    std::string GetClassName() const override { return "TRestDetectorAmplificationResponse"; }
    /// \brief Reads the configuration; throws if a value is out of range (no fieldMaps, more dielectricMaps
    /// than fieldMaps, fieldMapOffset not 3 values, periodicity with other letters than x, y, z, start not
    /// above anode, electrons, bins not > 0, penningRate not in [0, 1]).
    void LoadConfig() override;
    /// \brief Nothing to initialise.
    void Initialize() override {}

    /// \brief Sets the response from the avalanches of `primaries` electrons: the size (electrons at the
    /// anode, > 0) and time (us) of each one that reached the anode, and where the others ended.
    void SetResult(int primaries, const std::vector<std::array<double, 2>>& avalanches, const Losses& losses);
    /// \brief SetResult without losses
    void SetResult(int primaries, const std::vector<std::array<double, 2>>& avalanches) {
        SetResult(primaries, avalanches, Losses{});
    }

    /// \brief Fraction of the primary electrons whose avalanche reaches the anode (0 without primaries)
    double GetTransparency() const;
    /// \brief Polya parameter theta of the collected avalanches from their moments, sigma^2 / G^2 = 1 / (1 +
    /// theta) (0 without variance)
    double GetTheta() const;
    /// \brief Size of an avalanche sampled from its histogram (uniform within the bin)
    double SampleSize(TRandom& random) const;
    /// \brief Time (us) from the start plane to the anode sampled from its histogram (uniform within the bin)
    double SampleTime(TRandom& random) const;

    /// \brief Writes the response to a file (opened with RECREATE or UPDATE): its metadata (configuration,
    /// medium, conditions, counts, histograms and generator) in RESTMetadataStore, and the histograms as TH1D
    /// `AmplificationSize` and `AmplificationTime` to draw them.
    void Write(TFile* file, const YAML::Node& generator) const;
    /// \brief Reads the response of a file. Throws if it has none or it is not consistent.
    static TRestDetectorAmplificationResponse Read(TFile* file);
};

#endif
