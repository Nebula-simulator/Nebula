#ifndef __PARTICLE_H_
#define __PARTICLE_H_

/// Represents a particle (currently, always an electron) in the simulation.
struct particle
{
	vec3 pos;        ///< Position
	vec3 dir;        ///< Direction (unnormalised).
	real kin_energy; ///< Kinetic energy
};

#endif // __PARTICLE_H_
