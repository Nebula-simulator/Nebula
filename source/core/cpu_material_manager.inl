namespace nbl {

template<typename material_t>
cpu_material_manager<material_t>::~cpu_material_manager()
{
	for (material_t& mat : _materials)
		material_t::destroy(mat);
}

template<typename material_t>
void cpu_material_manager<material_t>::add(hdf5_file const & mat)
{
	if (_materials.size() == std::numeric_limits<material_index_t>::max())
		throw std::runtime_error("Too many materials submitted");
	_materials.push_back(material_t::create(mat));
}

template<typename material_t>
auto cpu_material_manager<material_t>::size() const
	-> material_index_t
{
	return static_cast<material_index_t>(_materials.size());
}

template<typename material_t>
material_t & cpu_material_manager<material_t>::operator[](material_index_t i)
{
	return _materials[i];
}
template<typename material_t>
material_t const & cpu_material_manager<material_t>::operator[](material_index_t i) const
{
	return _materials[i];
}

template<typename material_t>
bool cpu_material_manager<material_t>::is_physical(material_index_t material_idx) const
{
	return material_idx >= 0;
}

template<typename material_t>
real cpu_material_manager<material_t>::get_max_energy() const
{
	real max_energy = std::numeric_limits<real>::infinity();

	for (material_t const & mat : _materials)
	{
		const real e = mat.get_max_energy();
		if (e < max_energy)
			max_energy = e;
	}

	return max_energy;
}

} // namespace nbl
