#include "config/config.h"
#include "physics_config.h"
#include "core/material.h"
#include "common/cli_params.h"
#include "common/work_pool.h"
#include "common/time_log.h"
#include "io/output_stream.h"
#include "io/load_tri_file.h"
#include "io/load_pri_file.h"

#include "drivers/cpu/simple_cpu_driver.h"

#include "geometry/trilist.h"
#include "geometry/octree.h"

#include <iostream>
#include <iomanip>
#include <algorithm>
#include <numeric>
#include <thread>
#include <mutex>

// Main typedefs
using geometry_t = nbl::geometry::octree<false>;
using material_t = material<scatter_physics<false>>;

using driver = nbl::drivers::simple_cpu_driver<
	scatter_physics<false>,
	intersect_t,
	geometry_t
>;

int main(int argc, char** argv)
{
	// Print version information
	std::clog << "This is Nebula version "
		<< VERSION_MAJOR << "." << VERSION_MINOR << "." << VERSION_PATCH << "\n\n"
		"Physics models:\n";
	scatter_physics<false>::print_info(std::clog);
	intersect_t::print_info(std::clog);
	std::clog << "\n" << std::string(80, '-') << "\n\n";

	// Settings
	cli_params p("[options] <geometry.tri> <primaries.pri> [material0.mat] .. [materialN.mat]");
	p.add_option("energy-threshold", "Lowest energy to simulate", 0);
	p.add_option("seed", "Random seed", 0x78c7e39b);
	p.parse(argc, argv);
	const real energy_threshold = p.get_flag<real>("energy-threshold");
	const unsigned int seed = p.get_flag<unsigned int>("seed");


	// Setup time logging
	time_log timer;

	// Interpret command-line options
	std::vector<std::string> pos_flags = p.get_positional();
	if (pos_flags.size() < 3)
	{
		p.print_usage(std::clog);
		return 1;
	}

	std::mt19937 random_generator(seed);

	// Load geometry
	std::clog << "Loading geometry..." << std::endl;
	timer.start();
	std::vector<triangle> triangles = nbl::load_tri_file(pos_flags[0]);
	timer.stop("Loading triangles");

	if (triangles.empty())
	{
		std::clog << "Error: could not load triangles!\n";
		p.print_usage(std::clog);
		return 1;
	}
	// Sanity check with number of materials
	{
		int max_material = -1;
		for (triangle const & tri : triangles)
		{
			if (tri.material_in > max_material)
				max_material = tri.material_in;
			if (tri.material_out > max_material)
				max_material = tri.material_out;
		}

		if (max_material > int(pos_flags.size())-3)
		{
			std::clog << "Error: not enough materials provided for this geometry!\n"
				"  Expected " << (max_material+1) << " materials, " << (pos_flags.size()-2) << " provided.\n";
			return 1;
		}
		if (max_material < int(pos_flags.size())-3)
		{
			std::clog << "Warning: too many materials provided for this geometry!\n"
				"  Expected " << (max_material+1) << " materials, " << (pos_flags.size()-2) << " provided.\n";
		}
	}

	timer.start();
	geometry_t geometry = geometry_t::create(triangles);
	timer.stop("Building acceleration structure");


	// Load materials
	std::clog << "Loading materials..." << std::endl;
	timer.start();
	nbl::cpu_material_manager<material_t> materials;
	for (size_t parameter_idx = 2; parameter_idx < pos_flags.size(); ++parameter_idx)
	{
		nbl::hdf5_file material(pos_flags[parameter_idx]);
		materials.add(material);

		std::clog << "  Material " << (parameter_idx - 2) << ":\n"
			"    Name: " << material.get_property_string("name") << "\n"
			"    cstool version: " << material.get_property_string("cstool_version") << "\n";
	}
	timer.stop("Loading materials");


	// Load primaries
	std::clog << "Loading primary electrons..." << std::endl;
	timer.start();
	std::vector<particle> primaries; std::vector<int2> pixels;
	std::tie(primaries, pixels) = nbl::load_pri_file(pos_flags[1], geometry.AABB_min(), geometry.AABB_max(), materials.get_max_energy());
	timer.stop("Loading primary electrons");

	if (primaries.empty())
	{
		std::clog << "Error: could not load primary electrons!\n";
		p.print_usage(std::clog);
		return 1;
	}

	// The driver only accepts uint32 tags. So we make a map: simulation tag is
	// the index of the primary particle in the "primaries" / "pixels" array.
	std::vector<uint32_t> gpu_tags(primaries.size());
	std::iota(gpu_tags.begin(), gpu_tags.end(), 0); // Fill with 0, 1, ... tags.size()-1

	// This manages the work to be done (thread-safe).
	work_pool pool(primaries.data(), gpu_tags.data(), primaries.size());

	intersect_t inter;

	// Print debug data
	std::clog << "\n"
		<< "Loaded " << triangles.size() << " triangles.\n"
		<< "  min = {" << geometry.AABB_min().x << ", " << geometry.AABB_min().y << ", " << geometry.AABB_min().z << "}\n"
		<< "  max = {" << geometry.AABB_max().x << ", " << geometry.AABB_max().y << ", " << geometry.AABB_max().z << "}\n"
		<< "Loaded " << primaries.size() << " primaries.\n"
		<< "Loaded " << materials.size() << " materials.\n\n" << std::flush;

	// Prepare output file
	output_stream out_file("stdout");

	// Simulation loop
	auto sim_loop = [&pool, &out_file, &pixels,
		&geometry, &inter, &materials, energy_threshold](uint32_t seed)
	{
		driver d(
			inter, materials, geometry,
			energy_threshold, materials.get_max_energy(), seed);
		output_buffer buff(out_file, 1024*(7*sizeof(float) + 2*sizeof(int)));

		for (;;)
		{
			// Push new particles
			auto work_data = pool.get_work(1);

			if (std::get<2>(work_data) == 0)
				break;

			d.push(
				std::get<0>(work_data),  // particle*
				std::get<1>(work_data),  // tag*
				std::get<2>(work_data)); // number

			// Simulate a little
			d.simulate_to_end();

			// Flush output data
			d.flush_detected([&buff,&pixels](particle p, uint32_t t)
			{
				buff.add(std::array<float, 7>{
					p.pos.x, p.pos.y, p.pos.z,
					p.dir.x, p.dir.y, p.dir.z, p.kin_energy});
				buff.add(std::array<int, 2>{
					pixels[t].x, pixels[t].y});
			});
		}

		buff.flush();
	};

	// Simulation
	const auto n_threads = std::thread::hardware_concurrency();
	std::clog << "Creating " << n_threads << " CPU drivers" << std::endl;
	std::vector<std::thread> threads;

	timer.start();
	for (unsigned int i = 0; i < n_threads; ++i)
		threads.push_back(std::thread(sim_loop, random_generator()));

	// Progress indicator
	for (;;)
	{
		std::this_thread::sleep_for(std::chrono::seconds(1));
		auto primaries_to_go = pool.get_primaries_to_go();
		std::clog << " \rProgress "
			<< std::fixed << std::setprecision(2) << 100 * (1 - ((double)primaries_to_go / primaries.size())) << "%";
		if (primaries_to_go == 0)
			break;
	}

	for (auto& t : threads)
		t.join();

	timer.stop("Simulation");

	std::clog << "\n\n";
	timer.print(std::clog);
	return 0;
}
