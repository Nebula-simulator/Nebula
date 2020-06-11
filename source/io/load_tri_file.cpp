#include "../config/config.h"
#include "load_tri_file.h"
#include <fstream>
#include <sstream>

namespace nbl {

std::vector<triangle> load_tri_file(std::string const & filename)
{
	std::vector<triangle> triangle_vec;

	std::ifstream ifs(filename);
	if(!ifs.is_open())
		return triangle_vec;

	size_t line_num = 1;
	while(!ifs.eof())
	{
		std::string line_str;
		std::getline(ifs, line_str);
		line_num++;
		if(line_str.empty())
			continue;

		const size_t i = line_str.find_first_not_of(" \t");
		if(i < line_str.size() && line_str[i] == '#')
			continue;

		// Split columns in line into a vector
		std::stringstream ss;
		ss << line_str;
		std::vector<std::string> tag_vec;
		while(!ss.eof())
		{
			std::string tag_str;
			ss >> tag_str;
			if(!tag_str.empty())
				tag_vec.push_back(tag_str);
		}
		if(tag_vec.size() != 11)
		{
			std::ostringstream oss;
			oss << "invalid number of columns in line '" << line_num << "'";
			throw std::runtime_error(oss.str());
		}

		// Interpret line
		int in, out;
		in = std::stoi(tag_vec[0]);
		out = std::stoi(tag_vec[1]);
		vec3 A, B, C;
		A.x = (real)std::stod(tag_vec[2]); A.y = (real)std::stod(tag_vec[3]); A.z = (real)std::stod(tag_vec[4]);
		B.x = (real)std::stod(tag_vec[5]); B.y = (real)std::stod(tag_vec[6]); B.z = (real)std::stod(tag_vec[7]);
		C.x = (real)std::stod(tag_vec[8]); C.y = (real)std::stod(tag_vec[9]); C.z = (real)std::stod(tag_vec[10]);
		triangle_vec.push_back(triangle(A, B, C, in, out));
	}
	ifs.close();

	return triangle_vec;
}

} // namespace nbl
