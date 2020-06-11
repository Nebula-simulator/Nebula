#ifndef __HDF5_FILE_H_
#define __HDF5_FILE_H_

#include <string>
#include <hdf5.h>
#include "../common/units/unit_system.h"
#include "../common/util/range.h"
#include "../common/util/table_1D.h"
#include "../common/util/table_2D.h"
#include "../common/util/table_3D.h"

namespace nbl {

/**
 * \brief Simple interface to HDF5 files as we use them for storing material data.
 *
 * We can read:
 *   - Properties.
 *     Stored using the HDF5 attribute system, these can either be text or a
 *     number + unit.
 *   - N-dimensional data sets, with axis data stored using the HDF5 dimension
 *     scale API. The data set and axes may have a "unit" attribute. The data
 *     is returned as an ::nbl::util::table_1D, or its 2D or 3D equivalents.
 *     The dimension scales and the unit attribute must be read using separate
 *     functions in this class.
 */
class hdf5_file
{
public:
	/**
	 * \brief Constructor: open a file for reading
	 */
	hdf5_file(std::string const & filename);

	/**
	 * \brief Destructor.
	 */
	~hdf5_file();

	/**
	 * \brief Get the filename.
	 */
	std::string const & get_filename() const;

	/**
	 * \brief Check if a link with a certain name exists.
	 */
	bool exists(std::string const & name) const;

	/**
	 * \brief Get a string property. Throws `std::runtime_error` exception if
	 * not found.
	 *
	 * \param name Name of the property to be found.
	 */
	std::string get_property_string(std::string const & name) const;
	/**
	 * \brief Get a quantity property. Throws `std::runtime_error` exception if
	 * not found.
	 *
	 * \param name   Name of the property to be found.
	 * \param parser Class used to parse the quantity's units
	 */
	units::quantity<double> get_property_quantity(std::string const & name,
		units::unit_parser<double> const & parser = units::default_unit_parser()) const;

	/**
	 * \brief Get a string property. Returns default if not found.
	 *
	 * \param name     Name of the property to be found
	 * \param _default Value returned if property is not found
	 */
	std::string get_property_string(std::string const & name,
		std::string const & _default) const;
	/**
	 * \brief Get a quantity property. Returns default if not found.
	 *
	 * \param name     Name of the property to be found
	 * \param _default Value returned if property is not found
	 * \param parser   Class used to parse the quantity's units
	 */
	units::quantity<double> get_property_quantity(std::string const & name,
		units::quantity<double> const & _default,
		units::unit_parser<double> const & parser = units::default_unit_parser()) const;

	/**
	 * \brief Fill a 1D table with data.
	 *
	 * This function does not take care of the units of the data stored, and it
	 * does not set the min and max values correctly! These things must be done
	 * by the user, using get_unit(), get_lin_dimscale() and get_log_dimscale().
	 *
	 * \tparam T            Desired data type
	 * \param  dataset_name Name of the dataset
	 */
	template<typename T>
	util::table_1D<T, false> fill_table1D(std::string const & dataset_name) const;

	/**
	 * \brief Fill a 2D table with data.
	 *
	 * This function does not take care of the units of the data stored, and it
	 * does not set the min and max values correctly! These things must be done
	 * by the user, using get_unit(), get_lin_dimscale() and get_log_dimscale().
	 *
	 * \tparam T            Desired data type
	 * \param  dataset_name Name of the dataset
	 */
	template<typename T>
	util::table_2D<T, false> fill_table2D(std::string const & dataset_name) const;

	/**
	 * \brief Fill a 3D table with data.
	 *
	 * This function does not take care of the units of the data stored, and it
	 * does not set the min and max values correctly! These things must be done
	 * by the user, using get_unit(), get_lin_dimscale() and get_log_dimscale().
	 *
	 * \tparam T            Desired data type
	 * \param  dataset_name Name of the dataset
	 */
	template<typename T>
	util::table_3D<T, false> fill_table3D(std::string const & dataset_name) const;

	/**
	 * \brief Get a linearly spaced dimension scale.
	 *
	 * This function should be used when you expect the file to contain a linearly
	 * spaced dimension scale. This function reads it, checks that it is linear,
	 * and returns a ::nbl::util::linspace.
	 *
	 * \param dataset_name Name of the dataset
	 * \param dim          Index of the dimension of interest, starting at 0
	 * \param N_expected   Expected number of entries in the dimension scale.
	 *                     This should be equal to the width/height/depth of the
	 *                     associated table_xD.
	 * \param parser       Class used to parse the dimension scale's units
	 */
	util::linspace<units::quantity<double>> get_lin_dimscale(
		std::string const & dataset_name,
		int dim,
		int N_expected,
		units::unit_parser<double> const & parser = units::default_unit_parser()) const;

	/**
	 * \brief Get a logarithmically spaced dimension scale.
	 *
	 * This function should be used when you expect the file to contain a
	 * logarithmically spaced dimension scale. This function reads it, checks
	 * that it is logarithmic, and returns a ::nbl::util::geomspace.
	 *
	 * \param dataset_name Name of the dataset
	 * \param dim          Index of the dimension of interest, starting at 0
	 * \param N_expected   Expected number of entries in the dimension scale.
	 *                     This should be equal to the width/height/depth of the
	 *                     associated table_xD.
	 * \param parser       Class used to parse the dimension scale's units
	 */
	util::geomspace<units::quantity<double>> get_log_dimscale(
		std::string const & dataset_name,
		int dim,
		int N_expected,
		units::unit_parser<double> const & parser = units::default_unit_parser()) const;

	/**
	 * \brief Get the unit associated with a dataset.
	 *
	 * \param dataset_name Name of the dataset
	 * \param parser       Class used to parse the unit
	 */
	units::quantity<double> get_unit(std::string const & dataset_name,
		units::unit_parser<double> const & parser = units::default_unit_parser()) const;

private:
	std::string _filename;
	hid_t _file;
};

} // namespace nbl

#endif // __HDF5_FILE_H_
