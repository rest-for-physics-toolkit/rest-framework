#ifndef RESTGARFIELD_TRANSPORTSIMULATION_H
#define RESTGARFIELD_TRANSPORTSIMULATION_H

#include <filesystem>
#include <memory>

/// \file TransportSimulation.h
/// \brief Transport tables of restGarfield (drift velocity, diffusion, Townsend and attachment
/// coefficients), computed with Garfield++ (Magboltz for gases).

namespace Garfield {
class MediumGas;
class MediumMagboltz;
}  // namespace Garfield
class TRestDetectorMedium;
class TRestDetectorTransportTable;

/// \brief Exclusive lock (flock) of the directory of a medium file while it is read or written, so several
/// restGarfield processes can add tables to the same file (e.g. the media Ar99_iC4H10_1, Ar98_iC4H10_2...
/// of ArIso.root). It waits for the lock of another process and releases it when destroyed.
class MediumFileLock {
   public:
    explicit MediumFileLock(const std::filesystem::path& mediumFile);
    ~MediumFileLock();
    MediumFileLock(const MediumFileLock&) = delete;
    MediumFileLock& operator=(const MediumFileLock&) = delete;

   private:
    int fDescriptor = -1;
};

/// \brief Garfield++ gas (MediumMagboltz) of a gas medium at a pressure (bar) and temperature (K), without
/// gas table. Throws if the medium is not a gas or its components are not valid for Garfield.
/// \param medium Gas medium definition (1 to 6 components, by name and fraction).
/// \param pressure Pressure in bar.
/// \param temperature Temperature in K.
/// \return The Garfield gas, owned by the caller.
std::unique_ptr<Garfield::MediumMagboltz> CreateGas(const TRestDetectorMedium& medium, double pressure,
                                                    double temperature);

/// \brief Garfield++ gas whose transport tables are those of a TRestDetectorTransportTable (e.g. read from
/// the medium file), without running Magboltz: for the drift of electrons and ions with AvalancheMC or
/// DriftLineRKF in a field map. Only the tables are set (no cross sections: not for AvalancheMicroscopic).
///
/// The tables are those of `table` at `pressure`, scaled as TRestDetectorTransportTable does (E/p): its
/// field nodes at that pressure. The Garfield gas keeps its default pressure, the pressure of its tables,
/// so Garfield does not scale them again. The electron velocity, diffusion, Townsend and attachment
/// coefficients are set (Garfield keeps the logarithm of the coefficients, -30 for 0), and the ion mobility
/// and diffusion if the table has ions. Units: field V/cm, velocity cm/ns, diffusion sqrt(cm),
/// coefficients 1/cm, mobility cm2/(V ns).
///
/// \param medium Gas medium definition of the table (its components name the Garfield gas).
/// \param table Transport table with electrons.
/// \param pressure Pressure in bar.
/// \return The Garfield gas, owned by the caller.
/// \throws std::runtime_error if the medium is not a gas, its components are not valid for Garfield, or the
/// table has no electrons.
std::unique_ptr<Garfield::MediumGas> CreateGasFromTable(const TRestDetectorMedium& medium,
                                                        const TRestDetectorTransportTable& table,
                                                        double pressure);

/// \brief Computes the tables of a TRestDetectorTransportTable configuration with Garfield++: for each
/// condition, creates the Garfield medium of the type of the medium definition (MediumMagboltz, which
/// computes the gas table with Magboltz, MediumSilicon, ...), evaluates the transport parameters on the field
/// grid and adds the table to the medium file (TRestDetectorTransportTable).
///
/// One table per temperature (by default the one of the medium) at the pressure of the configuration. Only
/// the carriers (electron, hole, ion) that the Garfield medium provides are stored; electrons are required.
/// A table already in the medium file is skipped unless overwrite is set.
///
/// \param medium Medium definition: a gas, or a solid of a single component Si, GaAs, GaN, CdTe or diamond.
/// Its mediumFile receives the tables.
/// \param simulation Configuration of the transport tables (the `transport` section of
/// TRestDetectorTransportManager).
/// \param verbose Verbose Magboltz gas table generation (GenerateGasTable).
/// \param jobs Gases: the field nodes are split in up to `jobs` processes (Magboltz is not thread safe),
/// each one computing a contiguous part of the grid (at least 2 nodes); the table is the same.
/// \throws std::runtime_error if the medium is not supported, usePenning is set for a non-gas, there is no
/// electron transport or the medium file cannot be opened.
void RunTransportSimulation(const TRestDetectorMedium& medium, TRestDetectorTransportTable& simulation,
                            bool verbose, int jobs = 1);

#endif
