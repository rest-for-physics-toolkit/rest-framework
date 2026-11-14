#ifndef TRESTPARTICLES_H
#define TRESTPARTICLES_H

#include <string>

/// \brief Particle names and IDs of REST: the ID of a particle is its PDG encoding (Monte Carlo particle
/// numbering scheme, as Geant4 GetPDGEncoding), the one stored in the events (e.g. the tracks of
/// TRestGeant4Event, the hits of TRestDetectorHitsEvent; 0 = unknown). The names are those of Geant4 (fixed
/// here, so a configuration means the same with any Geant4 version):
///
/// - a table of the usual particles: gamma, e-, e+, mu-, mu+, tau-, tau+, the neutrinos (nu_e, anti_nu_e...),
///   pi+, pi-, pi0, kaon+, kaon-, kaon0, anti_kaon0, kaon0L, kaon0S, proton, anti_proton, neutron,
///   anti_neutron, lambda and the light nuclei deuteron, triton, He3 and alpha;
/// - the other nuclei as element symbol and mass number, e.g. Ar40 (PDG 100ZZZAAAI: 1000180400). A Geant4
///   excitation energy in brackets is accepted (Ar40[1460.859]: level I = 9, an excited state);
/// - a number is taken as a PDG encoding.
///
/// The exact names of the particles of a simulation (excited states with their energy) are those stored by
/// restG4 in TRestGeant4PhysicsInfo, which falls back to these functions.
///
/// \author JuanAn García <juanangp@unizar.es>
namespace REST_Particles {

/// \brief PDG encoding of a particle name (see the namespace documentation); 0 if the name is not known.
int GetParticleID(const std::string& name);

/// \brief Name of a PDG encoding: the table name, the nucleus (Ar40; an excited state Ar40*, where Geant4
/// writes its energy) or the number itself if unknown; "" for 0.
std::string GetParticleName(int id);

/// \brief Whether the PDG encoding is a nucleus (100ZZZAAAI, including deuteron, triton, He3 and alpha).
inline bool IsIon(int id) { return id >= 1000000000; }
/// \brief Atomic number Z of a nucleus (0 if not a nucleus).
inline int GetIonZ(int id) { return IsIon(id) ? (id / 10000) % 1000 : 0; }
/// \brief Mass number A of a nucleus (0 if not a nucleus).
inline int GetIonA(int id) { return IsIon(id) ? (id / 10) % 1000 : 0; }

/// \brief Whether a particle matches a name of a configuration: "ion" matches any nucleus, a nucleus name
/// its Z and A whatever its excitation level (Ar40 matches the excited Ar40), any other name its PDG
/// encoding. Throws std::invalid_argument if the name is not known.
bool Matches(int id, const std::string& name);

}  // namespace REST_Particles

#endif
