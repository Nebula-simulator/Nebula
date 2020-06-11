#include "../config/config.h"
#include "hdf5_file.h"
#include <hdf5_hl.h>
#include <array>
#include <vector>
#include <assert.h>

namespace nbl {

namespace
{
	/*
	 * Get size of a dataset's associated dataspace in each dimension.
	 * Verifies that the table has the expected number of dimensions.
	 */
	template<size_t N>
	inline std::array<hsize_t, N> get_size(hid_t dataset)
	{
		hid_t dataspace = H5Dget_space(dataset);

		if (H5Sget_simple_extent_ndims(dataspace) != N)
			throw std::runtime_error("Dataspace has unexpected dimension");

		std::array<hsize_t, N> dim;
		H5Sget_simple_extent_dims(dataspace, dim.data(), nullptr);

		H5Sclose(dataspace);
		return dim;
	}


	/*
	 * Get native string data type for use with the HDF5 library
	 */
	inline hid_t string_type()
	{
		hid_t tid = H5Tcopy(H5T_C_S1);
		H5Tset_size(tid, H5T_VARIABLE);
		return tid;
	}

	/*
	 * Get native data type for use with HDF5 library
	 */
	template<typename T>
	hid_t H5_mem_type();
	template<>
	hid_t H5_mem_type<float>() { return H5T_NATIVE_FLOAT; }
	template<>
	hid_t H5_mem_type<double>() { return H5T_NATIVE_DOUBLE; }


	/*
	 * Read an attribute, expecting a string
	 */
	inline std::string read_attribute_string(hid_t attribute)
	{
		char* buffer;
		hid_t string_t = string_type();

		H5Aread(attribute, string_t, &buffer);
		std::string result(buffer);

		H5Tclose(string_t);
		H5free_memory(buffer);
		return result;
	}

	/*
	 * Read an attribute, expecting a quantity (value + unit)
	 */
	inline units::quantity<double> read_attribute_quantity(
		hid_t attribute, units::unit_parser<double> const & parser)
	{
		struct data_struct
		{
			double value;
			char* unit;
		};

		// Create datatype
		hid_t string_t = string_type();
		hid_t data_type = H5Tcreate(H5T_COMPOUND, sizeof(data_struct));
		H5Tinsert(data_type, "value", HOFFSET(data_struct, value), H5T_NATIVE_DOUBLE);
		H5Tinsert(data_type, "unit", HOFFSET(data_struct, unit), string_t);

		// Read from file
		data_struct h5_result;
		H5Aread(attribute, data_type, &h5_result);

		// Copy to result struct
		const units::quantity<double> result = h5_result.value *
			parser.parse_unit(h5_result.unit);

		// Clean up and return
		H5free_memory(h5_result.unit);
		H5Tclose(data_type);
		H5Tclose(string_t);
		return result;
	}


	/*
	 * Get dimension scale data.
	 * The std::vector contains the actual values, the quantity contains the
	 * associated unit which acts like a multiplier for all values in the vector.
	 */
	std::pair<std::vector<double>, units::quantity<double>> get_dimscale(
		hid_t dataset,
		int dim,
		units::unit_parser<double> const & parser)
	{
		// Find out which dimension scale, if any, is attached.
		// If there are multiple, just get the first one.
		std::vector<double> data;
		units::quantity<double> unit;

		// Assemble data to be sent to the iteration function
		std::tuple<
			std::vector<double>*,
			units::quantity<double>*,
			units::unit_parser<double> const *>
		itdata { &data, &unit, &parser };

		// Iterate through the data scales attached, if any, and retrieve data
		H5DSiterate_scales(dataset, (unsigned int)dim, nullptr,
			[](hid_t /*did*/, unsigned /*dim*/, hid_t dsid, void* data) -> herr_t
			{
				// Recover itdata
				auto itdata = *reinterpret_cast<std::tuple<
					std::vector<double>*,
					units::quantity<double>*,
					units::unit_parser<double> const *>*>(data);

				// Read unit associated to dimension scale
				hid_t at_units = H5Aopen(dsid, "units", H5P_DEFAULT);
				*std::get<1>(itdata) = std::get<2>(itdata)->parse_unit(
					read_attribute_string(at_units));
				H5Aclose(at_units);

				// Read the dataset
				std::array<hsize_t, 1> dimensions = get_size<1>(dsid);

				std::get<0>(itdata)->resize(dimensions[0]);
				H5Dread(dsid, H5T_NATIVE_DOUBLE,
					H5S_ALL, H5S_ALL, H5P_DEFAULT, std::get<0>(itdata)->data());

				// Return 1 to stop the iteration
				return 1;
			}, &itdata);

		return { data, unit };
	}
} // anonymous namespace


hdf5_file::hdf5_file(std::string const & filename)
	: _filename(filename)
{
	H5Eset_auto(H5E_DEFAULT, nullptr, nullptr);
	_file = H5Fopen(filename.c_str(), H5F_ACC_RDONLY, H5P_DEFAULT);
	if (_file < 0)
		throw std::runtime_error("Unable to open file " + filename);
}
hdf5_file::~hdf5_file()
{
	 H5Fclose(_file);
}

std::string const & hdf5_file::get_filename() const
{
	return _filename;
}

bool hdf5_file::exists(std::string const & name) const
{
	return (H5LTpath_valid(_file, name.c_str(), 0) == 1);
}

std::string hdf5_file::get_property_string(std::string const & name) const
{
	hid_t attribute = H5Aopen(_file, name.c_str(), H5P_DEFAULT);
	if (attribute < 0)
		throw std::runtime_error("Could not read string property '" + name + '\'');

	auto result = read_attribute_string(attribute);
	H5Aclose(attribute);
	return result;
}
units::quantity<double> hdf5_file::get_property_quantity(
	std::string const & name,
	units::unit_parser<double> const & parser) const
{
	hid_t attribute = H5Aopen(_file, name.c_str(), H5P_DEFAULT);
	if (attribute < 0)
		throw std::runtime_error("Could not read quantity property '" + name + '\'');

	auto result = read_attribute_quantity(attribute, parser);
	H5Aclose(attribute);
	return result;
}

std::string hdf5_file::get_property_string(
	std::string const & name,
	std::string const & _default) const
{
	hid_t attribute = H5Aopen(_file, name.c_str(), H5P_DEFAULT);
	if (attribute < 0)
		return _default;

	auto result = read_attribute_string(attribute);
	H5Aclose(attribute);
	return result;
}
units::quantity<double> hdf5_file::get_property_quantity(
	std::string const & name,
	units::quantity<double> const & _default,
	units::unit_parser<double> const & parser) const
{
	hid_t attribute = H5Aopen(_file, name.c_str(), H5P_DEFAULT);
	if (attribute < 0)
		return _default;

	auto result = read_attribute_quantity(attribute, parser);
	H5Aclose(attribute);
	return result;
}

template<typename T>
util::table_1D<T, false> hdf5_file::fill_table1D(std::string const & dataset_name) const
{
	hid_t dataset = H5Dopen(_file, dataset_name.c_str(), H5P_DEFAULT);
	if (dataset < 0)
		throw std::runtime_error("Could not read table '" + dataset_name + '\'');

	std::array<hsize_t, 1> dimensions = get_size<1>(dataset);
	assert(dimensions[0] < std::numeric_limits<int>::max());

	auto table = util::table_1D<T, false>::create(0, 1, int(dimensions[0]));
	H5Dread(dataset, H5_mem_type<T>(),
		H5S_ALL, H5S_ALL, H5P_DEFAULT, table.data());

	H5Dclose(dataset);
	return table;
}

template<typename T>
util::table_2D<T, false> hdf5_file::fill_table2D(std::string const & dataset_name) const
{
	hid_t dataset = H5Dopen(_file, dataset_name.c_str(), H5P_DEFAULT);
	if (dataset < 0)
		throw std::runtime_error("Could not read table '" + dataset_name + '\'');

	std::array<hsize_t, 2> dimensions = get_size<2>(dataset);
	assert(dimensions[0] < std::numeric_limits<int>::max());
	assert(dimensions[1] < std::numeric_limits<int>::max());

	auto table = util::table_2D<T, false>::create(
		0, 1, int(dimensions[0]),
		0, 1, int(dimensions[1]));
	H5Dread(dataset, H5_mem_type<T>(),
		H5S_ALL, H5S_ALL, H5P_DEFAULT, table.data());

	H5Dclose(dataset);
	return table;
}

template<typename T>
util::table_3D<T, false> hdf5_file::fill_table3D(std::string const & dataset_name) const
{
	hid_t dataset = H5Dopen(_file, dataset_name.c_str(), H5P_DEFAULT);
	if (dataset < 0)
		throw std::runtime_error("Could not read table '" + dataset_name + '\'');

	std::array<hsize_t, 3> dimensions = get_size<3>(dataset);
	assert(dimensions[0] < std::numeric_limits<int>::max());
	assert(dimensions[1] < std::numeric_limits<int>::max());
	assert(dimensions[2] < std::numeric_limits<int>::max());

	auto table = util::table_3D<T, false>::create(
		0, 1, int(dimensions[0]),
		0, 1, int(dimensions[1]),
		0, 1, int(dimensions[2]));
	H5Dread(dataset, H5_mem_type<T>(),
		H5S_ALL, H5S_ALL, H5P_DEFAULT, table.data());

	H5Dclose(dataset);
	return table;
}

util::linspace<units::quantity<double>> hdf5_file::get_lin_dimscale(
	std::string const & dataset_name,
	int dim,
	int N_expected,
	units::unit_parser<double> const & parser) const
{
	hid_t dataset = H5Dopen(_file, dataset_name.c_str(), H5P_DEFAULT);
	if (dataset < 0)
		throw std::runtime_error("Could not read table '" + dataset_name + '\'');

	auto dimscale_data = get_dimscale(dataset, dim, parser);

	// If no dimension scale was found, fill between 0 and 1
	if (dimscale_data.first.size() == 0)
		return util::linspace<units::quantity<double>>(0*units::dimensionless, 1*units::dimensionless, N_expected);

	// Verify that the size is correct
	if (dimscale_data.first.size() != (size_t)N_expected)
		throw std::runtime_error("Dimension scale has unexpected size.");

	// Verify that the spacing is actually linear
	util::linspace<double> tmp_space(
		dimscale_data.first.front(),
		dimscale_data.first.back(),
		N_expected);
	for (int i = 0; i < N_expected; ++i)
		if (tmp_space[i] - dimscale_data.first[i] > EPSILON)
			throw std::runtime_error("Dataset " + dataset_name + ": Linear dimension scale expected");

	// Return linspace
	H5Dclose(dataset);
	return util::linspace<units::quantity<double>>(
		dimscale_data.first.front() * dimscale_data.second,
		dimscale_data.first.back() * dimscale_data.second,
		N_expected);
}

util::geomspace<units::quantity<double>> hdf5_file::get_log_dimscale(
	std::string const & dataset_name,
	int dim,
	int N_expected,
	units::unit_parser<double> const & parser) const
{
	hid_t dataset = H5Dopen(_file, dataset_name.c_str(), H5P_DEFAULT);
	if (dataset < 0)
		throw std::runtime_error("Could not read table '" + dataset_name + '\'');

	auto dimscale_data = get_dimscale(dataset, dim, parser);

	// Verify that the size is correct
	if (dimscale_data.first.size() != (size_t)N_expected)
		throw std::runtime_error("Dimension scale has unexpected size.");

	// Verify that the spacing is actually logarithmic
	util::geomspace<double> tmp_space(
		dimscale_data.first.front(),
		dimscale_data.first.back(),
		N_expected);
	for (int i = 0; i < N_expected; ++i)
		if (tmp_space[i] - dimscale_data.first[i] > EPSILON)
			throw std::runtime_error("Dataset " + dataset_name + ": Logarithmic dimension scale expected");

	// Return geomspace
	H5Dclose(dataset);
	return util::geomspace<units::quantity<double>>(
		dimscale_data.first.front() * dimscale_data.second,
		dimscale_data.first.back() * dimscale_data.second,
		N_expected);
}


units::quantity<double> hdf5_file::get_unit(
	std::string const & dataset_name,
	units::unit_parser<double> const & parser) const
{
	hid_t dataset = H5Dopen(_file, dataset_name.c_str(), H5P_DEFAULT);
	if (dataset < 0)
		throw std::runtime_error("Could not read table '" + dataset_name + '\'');

	hid_t at_units = H5Aopen(dataset, "units", H5P_DEFAULT);
	auto unit = parser.parse_unit(read_attribute_string(at_units));
	H5Aclose(at_units);
	H5Dclose(dataset);
	return unit;
}

// Explicit instantiations
template util::table_1D<float, false> hdf5_file::fill_table1D(std::string const &) const;
template util::table_2D<float, false> hdf5_file::fill_table2D(std::string const &) const;
template util::table_3D<float, false> hdf5_file::fill_table3D(std::string const &) const;
template util::table_1D<double, false> hdf5_file::fill_table1D(std::string const &) const;
template util::table_2D<double, false> hdf5_file::fill_table2D(std::string const &) const;
template util::table_3D<double, false> hdf5_file::fill_table3D(std::string const &) const;

} // namespace nbl
