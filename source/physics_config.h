#ifndef __PHYSICS_CONFIG_H_
#define __PHYSICS_CONFIG_H_

#include <iostream>
#include "core/particle.h"
#include "common/util/table_1D.h"
#include "common/util/table_2D.h"
#include "common/util/table_3D.h"
#include "common/util/random.h"
#include "core/material.h"
#include "physics/kieft/inelastic.h"
#include "physics/kieft/elastic.h"
#include "physics/full_penn.h"
#include "physics/boundary_intersect.h"

/*
 * Physics definitions below
 */

// Penn inelastic model
template<bool gpu_flag>
using inelastic_scatter = nbl::scatter::full_penn<gpu_flag,
	true  // Generate secondaries
>;

// Kieft & Bosch inelastic model
//template<bool gpu_flag>
//using inelastic_scatter = nbl::scatter::kieft_inelastic<gpu_flag,
//	true, // Optical phonon loss
//	true, // Generate secondaries
//	true, // Random instantaneous momentum for secondary
//	true  // Momentum conservation
//>;

// Kieft & Bosch elastic model
template<bool gpu_flag>
using elastic_scatter = nbl::scatter::kieft_elastic<gpu_flag,
	true, // Acoustic phonon loss
	true  // Atomic recoil loss
>;

// Material boundary intersection
using intersect_t = boundary_intersect<
	true, // Quantum-mechanical (probabilistic) transmission
	true, // Interface refraction
	false // Kieft & Bosch empirical interface absorption
>;


// Putting it all together
template<bool gpu_flag>
using scatter_physics = scatter_list<
	inelastic_scatter<gpu_flag>,
	elastic_scatter<gpu_flag>
>;

#endif
