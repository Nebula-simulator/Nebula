#ifndef __KIEFT_INELASTIC_H_
#define __KIEFT_INELASTIC_H_

namespace nbl { namespace scatter {

/**
 * \brief Inelastic scattering, according to the Kieft & Bosch model.
 *
 *   - See doi:10.1088/0022-3727/41/21/215310 (Kieft paper)
 *   - See doi:10.4233/uuid:f214f594-a21f-4318-9f29-9776d60ab06c (Verduin thesis)
 *
 * \tparam gpu_flag               Is the code to be run on a GPU?
 * \tparam optical_phonon_loss    Assume optical phonon loss for energy loss less than band gap
 * \tparam generate_secondary     Generate secondary electrons
 * \tparam instantaneous_momentum Consider instantaneous momentum for SE
 * \tparam momentum_conservation  Obey conservation of momentum
 */
template<bool gpu_flag,
	bool optical_phonon_loss = true,
	bool generate_secondary = true,
	bool instantaneous_momentum = true,
	bool momentum_conservation = true>
class kieft_inelastic
{
public:
	/**
	 * \brief Indicate when this class generates secondary electrons
	 */
	constexpr static bool may_create_se = generate_secondary;

	/**
	 * \brief Print diagnostic info
	 */
	static void print_info(std::ostream& stream)
	{
		stream << std::boolalpha <<
			" * Kieft & Bosch inelastic model\n"
			"   Options:\n"
			"     - Optical phonon loss: " << optical_phonon_loss << "\n"
			"     - Generate secondary electrons: " << generate_secondary << "\n"
			"     - Instantaneous momentum transfer: " << instantaneous_momentum << "\n"
			"     - Momentum conservation: " << momentum_conservation << "\n";
	}

	/**
	 * \brief Get maximal energy, in eV
	 */
	PHYSICS real get_max_energy() const
	{
		return expr(_log_imfp_table.get_scalemax());
	}

	/**
	 * \brief Sample a random free path length
	 */
	inline PHYSICS real sample_path(particle const & this_particle, util::random_generator<gpu_flag> & rng) const
	{
		// Get inverse mean free path for this kinetic energy
		const real imfp = expr(_log_imfp_table.get(logr(this_particle.kin_energy)));

		// Draw a distance
		return rng.exponential(1 / imfp);
	}

	/**
	 * \brief Perform a scattering event
	 */
	template<typename particle_manager>
	inline PHYSICS void execute(
		particle_manager& particle_mgr,
		typename particle_manager::particle_index_t particle_idx,
		util::random_generator<gpu_flag>& rng) const
	{
		// Retrieve current particle from global memory
		auto this_particle = particle_mgr[particle_idx];

		// draw a random zero-momentum energy loss of the primary electron
		real omega0;
		{// see thesis T.V. Eq. 3.82.
			const real x = logr(this_particle.kin_energy);
			const real y = rng.unit();
			omega0 = expr(_log_icdf_table.get(x, y));
		}

		// draw a random binding energy of the secondary electron.
		real binding;
		{
			const real x = logr(omega0);
			const real y = rng.unit();
			binding = _ionisation_table.get_rounddown(x, y);
		}

		// draw a random total energy loss for the primary electron
		real omega;
		{// see thesis T.V. Eq. 3.85.
			real omega_max = 0.5_r*(this_particle.kin_energy + omega0 - _fermi); // upper limit of eq. 9 in Ashley, but corrected for the fermi energy
			real omega_min = omega0;
			real w0 = minr(omega0 - 1, maxr(0, binding) - _fermi);
			if (this_particle.kin_energy > 2*omega0)
			{
				// equation 10 in Ashley
				omega_min = 0.5_r*(this_particle.kin_energy + omega0
					- sqrtr(this_particle.kin_energy*(this_particle.kin_energy - 2*omega0)));
				w0 = omega0;
			}

			const real U = rng.unit();
			if ((w0 > 0) && (omega_min > w0) && (omega_min < omega_max)) {
				// For nonzero binding energy, sample omega according to equation 7 in Ashley,
				// using the lower and upper limits as defined above.
				// For inner-shell ionization (Ebind > 50 eV) we substitute the Fermi-energy corrected
				// binding energy for omegaprime (so that the differential cross section becomes inversely
				// proportional to both the total energy transfer and the kinetic energy of the secondary
				// electron).
				const real f_min = 1 / w0*logr((omega_min - w0) / omega_min);
				const real f_max = 1 / w0*logr((omega_max - w0) / omega_max);
				omega = -w0 / expm1r(w0*(f_min*(1 - U) + f_max*U));
			}
			else {
				// In some cases (typically only occuring for binding < 50 eV) we get omega_min > omega_max.
				// This is due to our Fermi energy correction in the definition of omega_max. Physically, this
				// means that momentum cannot be conserved because the primary electron cannot have a final
				// kinetic energy that is lower than the Fermi energy. In this (relatively rare) case we have
				// to ignore momentum conservation and probe omega according to a 1/(omega)^2 distribution
				// with omega0 and omega_max as lower and upper limits respectively.
				omega = omega0*omega_max / (omega0*(1 - U) + omega_max*U);
			}
		}

		// special cases if there is no binding energy:
		if (binding < 0) {
			if (_band_gap < 0) {
				// TODO for metals: excitation of a fermi sea electron
			}
			else if (omega0 > _band_gap) {
				// electron excitation across the band gap (see page 78 thesis T.V.)
				binding = _band_gap;
			}
			else {
				// sub-band gap energy loss in semiconductors and insulators (see page 78 thesis T.V.)
				// energy loss due to longitudinal optical phonon excitation is assumed
				// update energy and EXIT
				if (optical_phonon_loss)
				{
					this_particle.kin_energy -= omega0;
					particle_mgr[particle_idx] = this_particle;
				}
				return;
			}
		}
		binding = maxr(0, binding);

		// Normalise direction
		normalise(this_particle.dir);

		// determine random normal vector to determine the scattering direction
		vec3 normal_dir = normalised(make_normal_vec(this_particle.dir, rng.phi()));

		// Determine the inelastic scattering angles.
		// We have strictly followed the method of Ivanchenko (and thus Kieft and Bosch)
		// See for more details thesis T.V. page 80-85.
		const real _K = this_particle.kin_energy - _fermi + 2*binding; // see thesis T.V. Eq. 3.105
		const real dK = binding + omega;        // see thesis T.V. Eq. 3.106
		const real cos_theta = sqrtr(dK / _K);  // see thesis T.V. Eq. 3.100
		const real sin_theta = sqrtr(1 - saturater(cos_theta*cos_theta));

		// Determine initial secondary direction (see thesis T.V. Eq. 3.107)
		// The initial direction is determined by assuming that the secondary electron
		// is at rest
		vec3 secondary_dir = this_particle.dir*cos_theta + normal_dir*sin_theta;

		// Add (optional) direction to account for the (intrinsic) instantaneous momentum
		//  of the secondary electron.
		// See thesis T.V. Eq. 3.108
		if (instantaneous_momentum)
		{
			normalise(secondary_dir);
			secondary_dir += sqrtr(binding / dK) * rng.uniform_vector();
		}

		// ensure proper normalization of the secondary directional vector.
		normalise(secondary_dir);

		if (generate_secondary)
		{
			particle secondary_particle;
			secondary_particle.kin_energy = _fermi + omega - binding; // See thesis T.V. Eq. 3.86
			secondary_particle.pos = this_particle.pos;
			secondary_particle.dir = secondary_dir;

			particle_mgr.create_secondary(particle_idx, secondary_particle);
		}

		this_particle.kin_energy -= omega;

		if (momentum_conservation)
		{
			// primary direction determined by non-relativistic momentum-conservation, i.e.:
			//   sin(theta)*primary_dir_2 = primary_dir - cos(theta)*secondary_dir;
			// See thesis T.V. Eq. 3.111
			this_particle.dir -= cos_theta * secondary_dir;
		}

		// Store the scattered particle in memory
		particle_mgr[particle_idx] = this_particle;
	}


	/**
	 * \brief Create, given a material file.
	 */
	static CPU kieft_inelastic create(hdf5_file const & mat)
	{
		if (!mat.exists("kieft/inelastic"))
			throw std::runtime_error("Kieft inelastic model not found in "
				"material file " + mat.get_filename());

		kieft_inelastic inel;

		inel._fermi = static_cast<real>(mat.get_property_quantity("fermi") / units::eV);
		inel._band_gap = static_cast<real>(mat.get_property_quantity("band_gap", -1*units::eV) / units::eV);

		{
			inel._log_imfp_table = mat.fill_table1D<real>("/kieft/inelastic/imfp");
			auto K_range = mat.get_log_dimscale("/kieft/inelastic/imfp", 0,
				inel._log_imfp_table.width());
			inel._log_imfp_table.set_scale(
				(real)std::log(K_range.front()/units::eV), (real)std::log(K_range.back()/units::eV));
			inel._log_imfp_table.mem_scope([&](real* imfp_vector)
			{
				const real unit = real(mat.get_unit("/kieft/inelastic/imfp") * units::nm);
				for (size_t x = 0; x < K_range.size(); ++x)
				{
					imfp_vector[x] = std::log(imfp_vector[x] * unit);
				}
			});
		}

		{
			inel._log_icdf_table = mat.fill_table2D<real>("/kieft/inelastic/w0_icdf");
			auto K_range = mat.get_log_dimscale("/kieft/inelastic/w0_icdf", 0,
				inel._log_icdf_table.width());
			auto P_range = mat.get_lin_dimscale("/kieft/inelastic/w0_icdf", 1,
				inel._log_icdf_table.height());
			inel._log_icdf_table.set_scale(
				(real)std::log(K_range.front()/units::eV), (real)std::log(K_range.back()/units::eV),
				(real)P_range.front(), (real)P_range.back());
			inel._log_icdf_table.mem_scope([&](real** icdf_vector)
			{
				const real unit = real(mat.get_unit("/kieft/inelastic/w0_icdf") / units::eV);
				const auto fermi = mat.get_property_quantity("fermi");
				for (size_t x = 0; x < K_range.size(); ++x)
				{
					auto K = K_range[x];
					for (size_t y = 0; y < P_range.size(); ++y)
					{
						icdf_vector[x][y] = (real)std::log(std::max(0.0, std::min<double>(
							(K - fermi) / units::eV,
							icdf_vector[x][y] * unit
						)));
					}
				}
			});
		}

		{
			inel._ionisation_table = mat.fill_table2D<real>("/kieft/ionization/binding_icdf");
			auto K_range = mat.get_log_dimscale("/kieft/ionization/binding_icdf", 0,
				inel._ionisation_table.width());
			auto P_range = mat.get_lin_dimscale("/kieft/ionization/binding_icdf", 1,
				inel._ionisation_table.height());
			inel._ionisation_table.set_scale(
				(real)std::log(K_range.front()/units::eV), (real)std::log(K_range.back()/units::eV),
				(real)P_range.front(), (real)P_range.back());
			inel._ionisation_table.mem_scope([&](real** ionisation_vector)
			{
				// Get outer_shells vector, containing outer-shell energies sorted from low to high.
				std::vector<units::quantity<double>> outer_shells;
				{
					const auto outer_shell_table = mat.fill_table1D<double>("/kieft/ionization/outer_shells");
					const auto unit = mat.get_unit("/kieft/ionization/outer_shells");
					for (size_t i = 0; i < outer_shell_table.width(); ++i)
					{
						const units::quantity<double> value = outer_shell_table(i) * unit;
						if (value < 100*units::eV)
							outer_shells.push_back(value);
					}
					std::sort(outer_shells.begin(), outer_shells.end());
				}

				// Create the simulation table
				const auto unit = mat.get_unit("/kieft/ionization/binding_icdf");
				for (size_t x = 0; x < K_range.size(); ++x)
				{
					auto K = K_range[x];
					for (size_t y = 0; y < P_range.size(); ++y)
					{
						const real P = (real)P_range[y];

						// Magic KB number
						const units::quantity<double> margin = 10 * units::eV;

						units::quantity<double> binding = -1 * units::eV;
						if (K > 100*units::eV)
						{
							binding = inel._ionisation_table.get_rounddown(
								(real)std::log((K+margin)/units::eV), P) * unit;
							if (binding < 50*units::eV || !std::isfinite(binding.value))
								binding = -1 * units::eV;
						}
						if (binding < 0*units::eV)
						{
							// Find largest outer shell less than or equal to K
							auto outer_shell_iterator = std::lower_bound(outer_shells.rbegin(),
								outer_shells.rend(), K, std::greater<units::quantity<double>>{});
							if (outer_shell_iterator != outer_shells.rend())
								binding = *outer_shell_iterator;
						}

						ionisation_vector[x][y] = real(binding / units::eV);
					}
				}
			});
		}

		return inel;
	}

	/**
	 * \brief Clone from another instance.
	 */
	template<bool source_gpu_flag>
	static CPU kieft_inelastic create(kieft_inelastic<source_gpu_flag, optical_phonon_loss, generate_secondary, instantaneous_momentum, momentum_conservation> const & source)
	{
		kieft_inelastic target;

		target._fermi = source._fermi;
		target._band_gap = source._band_gap;

		target._log_imfp_table = util::table_1D<real, gpu_flag>::create(source._log_imfp_table);
		target._log_icdf_table = util::table_2D<real, gpu_flag>::create(source._log_icdf_table);
		target._ionisation_table = util::table_2D<real, gpu_flag>::create(source._ionisation_table);

		return target;
	}

	/**
	 * \brief Dealllocate data held by an instance of this class.
	 */
	static CPU void destroy(kieft_inelastic & inel)
	{
		util::table_1D<real, gpu_flag>::destroy(inel._log_imfp_table);
		util::table_2D<real, gpu_flag>::destroy(inel._log_icdf_table);
		util::table_2D<real, gpu_flag>::destroy(inel._ionisation_table);
	}

private:
	/**
	 * \brief Table storing the inverse mean free path as function of energy.
	 *
	 * Actually, stores `log(inverse mean free path / nm^-1)` as function of
	 * `log(kinetic energy / eV)`.
	 */
	util::table_1D<real, gpu_flag> _log_imfp_table;

	/**
	 * \brief Table storing the probability distribution for "zero-momentum
	 * energy loss", omega prime.
	 *
	 * "ICDF" is short for Inverse Cumulative Distribution Function, which is
	 * what this table stores.
	 *
	 * Specifically, stores `log(omega' / eV)` as function of
	 *   - x axis: `log(kinetic energy / eV)`
	 *   - y axis: cumulative probability (between 0 and 1)
	 */
	util::table_2D<real, gpu_flag> _log_icdf_table;

	/**
	 * \brief Ionization table, representing probability of ionizing a given
	 * inner or outer shell.
	 *
	 * Specifically, stores `binding energy / eV` as function of
	 *   - x axis: `log(kinetic energy / eV)`
	 *   - y axis: cumulative probability (between 0 and 1)
	 *
	 * Be sure to use .get_rounddown() for non-interpolated values
	 */
	util::table_2D<real, gpu_flag> _ionisation_table;

	real _fermi;    ///< Fermi energy (eV)
	real _band_gap; ///< Band gap (eV)

	template<bool, bool, bool, bool, bool>
	friend class kieft_inelastic;
};

}} // namespace nbl::scatter

#endif // __KIEFT_INELASTIC_H_
