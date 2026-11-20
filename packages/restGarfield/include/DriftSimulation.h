#ifndef RESTGARFIELD_DRIFTSIMULATION_H
#define RESTGARFIELD_DRIFTSIMULATION_H

/// \file DriftSimulation.h
/// \brief Drift responses of restGarfield in field maps (TRestDetectorDriftResponse), computed with
/// Garfield++ AvalancheMC.

class TRestDetectorDriftResponse;
class TRestDetectorMedium;
class TRestDetectorReadout;
class TRestDetectorTransportTable;

/// \brief Computes the drift response of a configuration (see TRestDetectorDriftResponse): the field map
/// (LoadCSTFieldMap) with the Garfield gas of the transport table at the pressure of the medium
/// (CreateGasFromTable), the drift region cut at the boundary plane (the readout normal must be along x, y or
/// z), and `electrons` electrons drifted with AvalancheMC from each node of the grid. An electron arrives if
/// it ends on the boundary plane; its displacement from the node projected on that plane (along PlaneAxes of
/// the readout normal) and its drift time fill the histograms of the node. The electrons of a node on or
/// below the boundary plane all arrive at once, without displacement (their drift is the amplification's); a
/// node outside the field map has no arrivals. Fills the response (nodes, medium, pressure, temperature and
/// transport table); nothing is written.
///
/// \param medium Gas medium of the response (its pressure and temperature).
/// \param table Transport table of the medium (e.g. read from its medium file).
/// \param readout Readout of the detector: distance to the readout and normal.
/// \param response Configuration of the response, filled.
/// \param verbose Print the progress.
/// \throws std::runtime_error if the field map cannot be loaded, the readout has no distance or its normal
/// is not along an axis, or the field map is not above the boundary plane.
void ComputeDriftResponse(const TRestDetectorMedium& medium, const TRestDetectorTransportTable& table,
                          const TRestDetectorReadout& readout, TRestDetectorDriftResponse& response,
                          bool verbose);

/// \brief Computes the drift response of a configuration and writes it to its responseFile (skipped if the
/// file exists, unless overwrite): the readout from the TRestDetectorReadoutManager of the `readouts` YAML
/// file (`readout`, or its data readout) and the transport table `transportTable` of the medium file (or its
/// only table), then ComputeDriftResponse.
///
/// \param medium Gas medium of the response, with its medium file.
/// \param response Configuration of the response (the `drift` section of TRestDetectorTransportManager).
/// \param verbose Print the progress.
/// \throws std::runtime_error if there is no readouts file or readout, the medium file has no transport
/// table (or several and no transportTable), or ComputeDriftResponse fails.
void RunDriftResponse(const TRestDetectorMedium& medium, TRestDetectorDriftResponse& response, bool verbose);

#endif
