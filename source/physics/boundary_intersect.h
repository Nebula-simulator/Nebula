#ifndef __BOUNDARY_INTERSECT_H_
#define __BOUNDARY_INTERSECT_H_

#include "../core/triangle.h"
#include "../core/cpu_material_manager.h"

/**
 * \brief Material boundary intersection.
 *
 * See  (Kieft paper)
 *
 * \tparam quantum_transmission Consider probabilistic transmission
 * \tparam interface_refraction Consider refraction at interface
 * \tparam interface_absorption Empirical interface absorption (Kieft et al.
 *                                  doi:10.1088/0022-3727/41/21/215310)
 */
template<
	bool quantum_transmission = true,
	bool interface_refraction = true,
	bool interface_absorption = false>
struct boundary_intersect
{
	/**
	 * \brief Print diagnostic info
	 */
	static void print_info(std::ostream& stream)
	{
		stream << std::boolalpha <<
			" * Material boundary crossing model\n"
			"   Options:\n"
			"     - Quantum mechanical transmission: " << quantum_transmission << "\n"
			"     - Interface refraction: " << interface_refraction << "\n"
			"     - Empirical interface absorption: " << interface_absorption << "\n";
	}

	/**
	 * \brief Perform intersection event.
	 */
	template<typename particle_manager, typename material_manager, bool gpu_flag>
	PHYSICS void execute(material_manager& material_mgr,
		particle_manager& particle_mgr, typename particle_manager::particle_index_t particle_idx,
		nbl::util::random_generator<gpu_flag>& rng) const
	{
		using material_index_t = typename material_manager::material_index_t;

		// Get particle data from particle manager
		auto this_particle = particle_mgr[particle_idx];
		const triangle this_triangle = *(particle_mgr.get_last_triangle(particle_idx));

		// Get particle direction, normal of the intersected triangle
		auto normalised_dir = normalised(this_particle.dir);
		auto last_triangle_normal = normalised(this_triangle.get_normal());

		// Get angle between direction of motion and triangle
		const real cos_theta = dot_product(last_triangle_normal, normalised_dir);


		// determine the material index in the following way:
		//   material_idx_in represents the current material
		//   material_idx_out represents the material when passing through the interface
		material_index_t material_idx_in, material_idx_out;
		if (cos_theta > 0)
		{
			material_idx_in = this_triangle.material_in;
			material_idx_out = this_triangle.material_out;
		}
		else
		{
			material_idx_in = this_triangle.material_out;
			material_idx_out = this_triangle.material_in;
		}


		// manage special cases for electron detection, electron mirrors and terminators.
		//   DETECTOR always detects.
		//   DETECTOR_LT/GE50 detect under certain circumstances.
		//     If not detected, they pass through as if no intersection event has taken place.
		switch (material_idx_out) {
		case nbl::special_materials::DETECTOR:
			particle_mgr.detect(particle_idx);
			return;
		case nbl::special_materials::DETECTOR_LT50:
			if (this_particle.kin_energy < 50)
				particle_mgr.detect(particle_idx);
			return;
		case nbl::special_materials::DETECTOR_GE50:
			if (this_particle.kin_energy >= 50)
				particle_mgr.detect(particle_idx);
			return;
		case nbl::special_materials::TERMINATOR:
			particle_mgr.terminate(particle_idx);
			return;
		case nbl::special_materials::MIRROR:
			this_particle.dir = normalised_dir - 2*last_triangle_normal*cos_theta;
			particle_mgr[particle_idx] = this_particle;
			return;
		case nbl::special_materials::NOP:
			return;
		default:
			break;
		}


		// determine the change in energy `dU` (in eV) when passing through the interface
		// see thesis T.V. Eq. 3.136
		real dU = 0;
		if (material_mgr.is_physical(material_idx_in)) {
			dU -= material_mgr[material_idx_in].barrier;
		}
		if (material_mgr.is_physical(material_idx_out)) {
			dU += material_mgr[material_idx_out].barrier;
		}

		// determine transmission probability (only if energy suffices)
		// see thesis T.V. Eq. 3.145
		if (this_particle.kin_energy*cos_theta*cos_theta + dU > 0)
		{
			const real s = sqrtr(1 + dU / (this_particle.kin_energy*cos_theta*cos_theta));
			const real T = (quantum_transmission
				? 4 * s / ((1 + s)*(1 + s))
				: 1);
			if (rng.unit() < T)
			{
				if (interface_refraction)
				{
					// determine the angle of refraction
					// see thesis T.V. Eq. 3.139
					this_particle.dir = (normalised_dir - last_triangle_normal*cos_theta)
						+ s * last_triangle_normal * cos_theta;
				}

				// if there is transmission, then adjust the kinetic energy,
				// update the current material index and EXIT.
				this_particle.kin_energy += dU;

				particle_mgr.set_material_index(particle_idx, material_idx_out);
				particle_mgr[particle_idx] = this_particle;
				return;
			}
		}

		// surface absorption? (this is in accordance with Kieft & Bosch code)
		// note that the default behaviour has this feature disabled
		if (interface_absorption)
		{
			if (dU < 0 && rng.unit() < expr(1 + 0.5_r*this_particle.kin_energy / dU))
			{
				particle_mgr.terminate(particle_idx);
				return;
			}
		}

		// the only remaining case is total internal reflection
		this_particle.dir = normalised_dir - 2 * last_triangle_normal*cos_theta;
		particle_mgr[particle_idx] = this_particle;
	}
};

#endif // __BOUNDARY_INTERSECT_H_
