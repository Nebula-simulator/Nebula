#include "config/config.h"
#include "physics_config.h"
#include "core/material.h"
#include "common/cli_params.h"
#include "common/work_pool.h"
#include "common/time_log.h"
#include "io/output_stream.h"
#include "io/load_tri_file.h"
#include "io/load_pri_file.h"

#include "drivers/gpu/gpu_driver.h"

#include "geometry/trilist.h"
#include "geometry/octree.h"

#include <iostream>
#include <iomanip>
#include <algorithm>
#include <numeric>
#include <thread>
#include <mutex>
#include <condition_variable>

// Main typedefs
using geometry_t = nbl::geometry::octree<true>;
using cpu_material_t = material<scatter_physics<false>>;
using gpu_material_t = material<scatter_physics<true>>;

using driver = nbl::drivers::gpu_driver<
	scatter_physics<true>,
	intersect_t,
	geometry_t
>;

// Get maximal energy accepted by all material files
real get_max_energy(std::vector<cpu_material_t> const & materials)
{
	if (materials.size() == 0)
		return 0;

	real max_energy = std::numeric_limits<real>::infinity();
	for (auto const & mat : materials)
	{
		const real e = mat.get_max_energy();
		if (e < max_energy)
			max_energy = e;
	}
	return max_energy;
}

struct worker_data
{
	worker_data(
		std::string const & outfilename,
		time_log& timer)
	: out_file(outfilename), timer(timer)
	{}

	std::condition_variable cv;
	std::mutex cv_m;

	nbl::geometry::octree_builder::linearized_octree geometry;
	std::vector<cpu_material_t> materials;
	work_pool primaries;
	std::vector<int2> pixels;

	time_log& timer;

	uint32_t prescan_size;
	real batch_factor;

	uint32_t capacity;
	uint32_t frame_size;
	uint32_t batch_size;
	real energy_threshold;

	std::vector<uint32_t> running_count;

	output_stream out_file;

	enum class status_t
	{
		init,
		geometry_loaded,
		materials_loaded,
		primaries_loaded,
		prescan_done
	};
	status_t status = status_t::init;
};

void worker_thread(worker_data& data,
	int gpu_id,
	typename driver::seed_t seed,
	bool do_prescan)
{
	cudaSetDevice(gpu_id);

	// Upload geometry
	std::unique_lock<std::mutex> lk(data.cv_m);
	data.cv.wait(lk, [&data]() { return data.status >= worker_data::status_t::geometry_loaded; });
	lk.unlock();

	geometry_t geometry = geometry_t::create(data.geometry);


	// Upload materials
	{
		std::unique_lock<std::mutex> lk(data.cv_m);
		data.cv.wait(lk, [&data]{ return data.status >= worker_data::status_t::materials_loaded; });
		lk.unlock();
	}

	std::vector<gpu_material_t> materials;
	for (auto const & cpumat : data.materials)
		materials.push_back(gpu_material_t(cpumat));


	// Create driver
	intersect_t inter;
	driver d(data.capacity,
		geometry, inter, materials,
		data.energy_threshold, seed);


	// Do prescan, if desired
	if (do_prescan)
	{
		{
			std::unique_lock<std::mutex> lk(data.cv_m);
			data.cv.wait(lk, [&data]{ return data.status >= worker_data::status_t::primaries_loaded; });
			lk.unlock();
		}


		data.timer.start();
		std::vector<std::pair<uint32_t, uint32_t>> prescan_stats; // Holds (running_count, detected_count)
		// Push first batch
		{
			auto work_data = data.primaries.get_work(data.prescan_size);
			auto particles_pushed = d.push(
				std::get<0>(work_data),  // particle*
				std::get<1>(work_data),  // tag*
				std::get<2>(work_data)); // number
			prescan_stats.push_back({ particles_pushed, 0 });
		}
		// Execute prescan
		while (prescan_stats.back().first > 0)
		{
			d.do_iteration();

			// TODO: this can be optimised to just one function with one loop.
			prescan_stats.push_back({ d.get_running_count(), d.get_detected_count() });
			std::clog << " \rExecuting pre-scan"
				<< " | running: " << prescan_stats.back().first
				<< " | detected: " << prescan_stats.back().second;
		}
		data.timer.stop("Prescan");

		// Find frame_size and batch_size based on the prescan stats.
		// frame_size is the iteration number where running_count was maximal.
		const uint32_t frame_size = 1 + (uint32_t)std::distance(prescan_stats.begin(),
			std::max_element(prescan_stats.begin(), prescan_stats.end(),
			[](std::pair<uint32_t, uint32_t> p1, std::pair<uint32_t, uint32_t> p2) -> bool
			{ return p1.first < p2.first; }));
		// Batch size
		uint32_t batch_size;
		{
			real accumulator = 0;
			for (uint32_t i = 2*frame_size; i < prescan_stats.size(); i += frame_size)
				accumulator += prescan_stats[i].first / real(data.prescan_size);
			accumulator += 2*prescan_stats[frame_size].first / real(data.prescan_size);
			accumulator += 2*prescan_stats[frame_size].second / real(data.prescan_size);
			batch_size = uint32_t(data.batch_factor*data.capacity / accumulator);
		}
		std::clog << "\nframe_size = " << frame_size << " | batch_size = " << batch_size << std::endl;


		// Notify other threads
		{
			std::lock_guard<std::mutex> lg(data.cv_m);
			data.frame_size = frame_size;
			data.batch_size = batch_size;
			data.status = worker_data::status_t::prescan_done;
		}
		data.cv.notify_all();
	}
	else
	{
		std::unique_lock<std::mutex> lk(data.cv_m);
		data.cv.wait(lk, [&data]{ return data.status >= worker_data::status_t::prescan_done; });
		lk.unlock();
	}


	// Start simulation
	d.allocate_input_buffers(data.batch_size);
	output_buffer buff(data.out_file, 1024*(7*sizeof(float) + 2*sizeof(int)));

	for (;;)
	{
		// Copy the simulation to buffer,
		// push new data into the simulation
		d.buffer_detected();
		d.push_to_simulation();
		cudaDeviceSynchronize();

		// Execute frame
		for (uint32_t i = 0; i < data.frame_size; ++i)
			d.do_iteration();

		// Push new batch asynchronously
		d.push_to_buffer(data.primaries);

		// Output detected electrons from buffer
		auto running_count = d.flush_buffered([&buff, &data](particle p, uint32_t t)
		{
			buff.add(std::array<float, 7>{
				p.pos.x, p.pos.y, p.pos.z,
				p.dir.x, p.dir.y, p.dir.z, p.kin_energy});
			buff.add(std::array<int, 2>{
				data.pixels[t].x, data.pixels[t].y});
		});
		data.running_count[gpu_id] = running_count;

		cudaDeviceSynchronize();

		if (running_count == 0 && data.primaries.done())
			break;
	}
	buff.flush();
}


int main(int argc, char** argv)
{
	// Print version information
	std::clog << "This is Nebula version "
		<< VERSION_MAJOR << "." << VERSION_MINOR << "." << VERSION_PATCH << "\n\n"
		"Physics models:\n";
	scatter_physics<true>::print_info(std::clog);
	intersect_t::print_info(std::clog);
	std::clog << "\n" << std::string(80, '-') << "\n\n";


	// Settings
	cli_params p("[options] <geometry.tri> <primaries.pri> [material0.mat] .. [materialN.mat]");
	p.add_option("energy-threshold", "Lowest energy to simulate", 0);
	p.add_option("capacity", "Electron capacity on the GPU", 1000000);
	p.add_option("prescan-size", "Number of electrons to use for prescan", 1000);
	p.add_option("batch-factor", "Multiplication factor for electron batch size", 0.9);
	p.add_option("seed", "Random seed", 0x78c7e39b);
	p.add_option("sort-primaries", "Sort primary electrons before simulation", false);
	p.parse(argc, argv);
	const real energy_threshold = p.get_flag<real>("energy-threshold");
	const uint32_t capacity = p.get_flag<uint32_t>("capacity");
	const uint32_t prescan_size = p.get_flag<uint32_t>("prescan-size");
	const real batch_factor = p.get_flag<real>("batch-factor");
	const unsigned int seed = p.get_flag<unsigned int>("seed");
	const bool sort_primaries = p.get_flag<bool>("sort-primaries");

	// Setup time logging
	time_log timer;

	// Interpret command-line options
	std::vector<std::string> pos_flags = p.get_positional();
	if (pos_flags.size() < 3 || capacity <= 0 || prescan_size <= 0 || batch_factor <= 0)
	{
		p.print_usage(std::clog);
		return 1;
	}

	std::mt19937 random_generator(seed);


	// Setup worker_data structure
	worker_data data("stdout", timer);
	data.capacity = capacity;
	data.energy_threshold = energy_threshold;
	data.prescan_size = prescan_size;
	data.batch_factor = batch_factor;


	// Launch worker threads: one for each GPU
	std::vector<std::thread> threads;
	int n_gpus;
	cudaGetDeviceCount(&n_gpus);
	std::clog << "Found " << n_gpus << " GPUs:\n";
	for (int i = 0; i < n_gpus; ++i)
	{
		cudaDeviceProp props;
		cudaGetDeviceProperties(&props, i);
		std::clog << "  - [" << i << "] " << props.name << '\n';
		threads.push_back(std::thread(worker_thread, std::ref(data), i, random_generator(), i == 0));
	}
	std::clog << std::endl;

	data.running_count.resize(n_gpus);


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
	data.geometry = nbl::geometry::octree_builder::linearized_octree(
		nbl::geometry::octree_builder::octree_root(triangles));
	timer.stop("Building acceleration structure");

	{
		std::lock_guard<std::mutex> lg(data.cv_m);
		data.status = worker_data::status_t::geometry_loaded;
	}
	data.cv.notify_all();


	// Load materials
	std::clog << "Loading materials..." << std::endl;
	timer.start();
	for (size_t parameter_idx = 2; parameter_idx < pos_flags.size(); ++parameter_idx)
	{
		nbl::hdf5_file material(pos_flags[parameter_idx]);
		data.materials.push_back(cpu_material_t(material));

		std::clog << "  Material " << (parameter_idx - 2) << ":\n"
			"    Name: " << material.get_property_string("name") << "\n"
			"    cstool version: " << material.get_property_string("cstool_version") << "\n";
	}
	timer.stop("Loading materials");

	{
		std::lock_guard<std::mutex> lg(data.cv_m);
		data.status = worker_data::status_t::materials_loaded;
	}
	data.cv.notify_all();


	// Load primaries
	std::clog << "Loading primary electrons..." << std::endl;
	timer.start();
	std::vector<particle> primaries;
	std::tie(primaries, data.pixels) = nbl::load_pri_file(pos_flags[1], data.geometry.AABB_min(), data.geometry.AABB_max(), get_max_energy(data.materials));
	if (sort_primaries)
		nbl::sort_pri_file(primaries, data.pixels);
	nbl::prescan_shuffle(primaries, data.pixels, prescan_size);
	timer.stop("Loading primary electrons");

	if (primaries.empty())
	{
		std::clog << "Error: could not load primary electrons!\n";
		p.print_usage(std::clog);
		return 1;
	}

	// The GPU driver only accepts uint32 tags. So we make a map: GPU tag is
	// the index of the primary particle in the "primaries" / "pixels" array.
	std::vector<uint32_t> gpu_tags(primaries.size());
	std::iota(gpu_tags.begin(), gpu_tags.end(), 0); // Fill with 0, 1, ... tags.size()-1

	data.primaries = work_pool(primaries.data(), gpu_tags.data(), primaries.size());


	{
		std::lock_guard<std::mutex> lg(data.cv_m);
		data.status = worker_data::status_t::primaries_loaded;
	}
	data.cv.notify_all();



	// Print debug data
	std::clog << "\n"
		<< "Loaded " << triangles.size() << " triangles.\n"
		<< "  min = {" << data.geometry.AABB_min().x << ", " << data.geometry.AABB_min().y << ", " << data.geometry.AABB_min().z << "}\n"
		<< "  max = {" << data.geometry.AABB_max().x << ", " << data.geometry.AABB_max().y << ", " << data.geometry.AABB_max().z << "}\n"
		<< "Loaded " << primaries.size() << " primaries.\n"
		<< "Loaded " << data.materials.size() << " materials.\n\n" << std::flush;




	// Simulation
	{
		std::unique_lock<std::mutex> lk(data.cv_m);
		data.cv.wait(lk, [&data]() { return data.status >= worker_data::status_t::prescan_done; });
		lk.unlock();
	}

	timer.start();
	for (;;)
	{
		std::this_thread::sleep_for(std::chrono::seconds(1));
		auto primaries_to_go = data.primaries.get_primaries_to_go();

		std::clog << " \rProgress "
			<< std::fixed << std::setprecision(2) << 100 * (1 - ((double)primaries_to_go / primaries.size()))
			<< "% | Running: ";
		for (int i = 0; i < n_gpus-1; ++i)
			std::clog << data.running_count[i] << " | ";
		std::clog << data.running_count[n_gpus-1];

		if (primaries_to_go == 0)
		{
			if (std::all_of(data.running_count.begin(), data.running_count.end(), [](uint32_t i) { return i==0; }))
				break;
		}
	}
	for (auto& t : threads)
		t.join();
	timer.stop("Simulation");

	cudaError_t err = cudaDeviceSynchronize();
	if (err == 0)
		std::clog << "\nSimulation successful!\n\n";
	else
		std::clog << "\nSimulation ended with CUDA error code " << err << "\n\n";

	timer.print(std::clog);
	return 0;
}
