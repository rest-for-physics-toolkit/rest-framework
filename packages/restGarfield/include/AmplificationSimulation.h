#ifndef RESTGARFIELD_AMPLIFICATIONSIMULATION_H
#define RESTGARFIELD_AMPLIFICATIONSIMULATION_H

/// \file AmplificationSimulation.h
/// \brief Amplification responses of restGarfield in the field maps of a unit cell
/// (TRestDetectorAmplificationResponse), computed with Garfield++ microscopic avalanches.

class TRestDetectorAmplificationResponse;
class TRestDetectorMedium;

/// \brief Computes the amplification response of a configuration (see TRestDetectorAmplificationResponse):
/// the field maps of the cell (LoadCSTFieldMaps, with their dielectric maps and the periodicity) with the
/// Garfield gas of the composition of the medium at its pressure and temperature (CreateGas, with the Penning
/// transfer if usePenning), and `electrons` avalanches (AvalancheMicroscopic) each started by one electron at
/// a random position over the cell on the plane z = start. The size of an avalanche is its number of
/// electrons ending below z = anode; one with electrons there is collected, with the mean time of those
/// electrons. Fills the response (medium, conditions, counts, histograms, the map of the cell by the start
/// point and the electrons at the anode out of the cell of their primary); nothing is written.
///
/// \param medium Gas medium of the response (its components, pressure and temperature).
/// \param response Configuration of the response, filled.
/// \param verbose Print the progress.
/// \param jobs Processes (fork) sharing the avalanches (Garfield++ is not thread safe), each with its own
/// random seed drawn from the generator of Garfield++; 1: in this process.
/// \throws std::runtime_error if the medium is not a gas, a field map cannot be loaded, the start or anode
/// planes are not inside the maps, the Penning transfer cannot be enabled, or a process fails.
void ComputeAmplificationResponse(const TRestDetectorMedium& medium,
                                  TRestDetectorAmplificationResponse& response, bool verbose, int jobs = 1);

/// \brief Computes the amplification response of a configuration (ComputeAmplificationResponse) and writes
/// it to its responseFile (skipped if the file exists, unless overwrite).
///
/// \param medium Gas medium of the response.
/// \param response Configuration of the response (the `amplification` section of
/// TRestDetectorTransportManager).
/// \param verbose Print the progress.
/// \param jobs Processes sharing the avalanches (`--jobs` of restGarfield).
/// \throws std::runtime_error if ComputeAmplificationResponse fails or the file cannot be written.
void RunAmplificationResponse(const TRestDetectorMedium& medium, TRestDetectorAmplificationResponse& response,
                              bool verbose, int jobs = 1);

#endif
