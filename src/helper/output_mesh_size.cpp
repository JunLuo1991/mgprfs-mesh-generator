// __START_OF_LICENSE__
// 
// Copyright (c) 2019 Jun Luo
// All rights reserved.
//
// This program is free software; you can redistribute it and/or
// modify it under the terms of the GNU General Public License as
// published by the Free Software Foundation; either version 3,
// or (at your option) any later version.
// 
// This program is distributed in the hope that it will be useful,
// but WITHOUT ANY WARRANTY; without even the implied warranty of
// MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
// GNU General Public License for more details.
// 
// You should have received a copy of the GNU General Public
// License along with this program; see the file LICENSE.  If not,
// see <http://www.gnu.org/licenses/>.
// 
// __END_OF_LICENSE__


/*!
@brief
This program reads a mesh model from stardard input and
output the mesh size.
*/

#include <CGAL/Cartesian.h>
#include <SPL/Array2.hpp>
#include <SPL/cgal/Delaunay_triangulation_2.hpp>
#include <SPL/cgal/Triangulation_hierarchy_2.hpp>
#include <SPL/rasterize.hpp>
#include <array>
#include <boost/program_options.hpp>
#include <iostream>
#include <string>

#include "pixel.hpp"

/*
@brief A structure for making the trinagulation
*/
struct Make_tri
{
  class Tri_vertex;
  struct Tri_traits
  {
    using Geom_kernel = CGAL::Filtered_kernel<CGAL::Cartesian<double>>;
    using Vertex0 = SPL::Triangulation_hierarchy_vertex_base_2<Tri_traits>;
    using Vertex = Tri_vertex;
    using Edge = SPL::Triangulation_edge_base_2<Tri_traits>;
    using Face = SPL::Triangulation_face_base_2<Tri_traits>;
    using Halfedge = SPL::Triangulation_halfedge_base_2<Tri_traits>;
  };

  class Tri_vertex : public Tri_traits::Vertex0
  {
  public:
    int z;
    pixel<int, 3> color;
  };
  using Triangulation0 = SPL::Delaunay_triangulation_2<Tri_traits>;
  using Triangulation = SPL::Triangulation_hierarchy_2<Triangulation0>;
};

using Tri = Make_tri::Triangulation;

// Read the off data into a triangulation from specified input stream
int input_off_model(std::istream& in, Tri& tri)
{

  auto set_vertex_z = [](Tri::Vertex_handle v, double z) { v->z = z; };
  auto set_vertex_color = [](Tri::Vertex_handle v,
                             const std::vector<double>& color, bool as_int) {
    assert(color.size() == 4);
    std::copy(color.begin(), color.begin() + 3, v->color.begin());
  };

  if (tri.input_off(in, set_vertex_z, nullptr, set_vertex_color, nullptr,
                    nullptr)) {
    std::cerr << "read OFF data failed!\n";
    return 1;
  }
  return 0;
}


int main(int argc, char** argv)
{
  std::string out_file_name;
  std::string in_file_name;

  namespace po = boost::program_options;
  po::options_description desc{ "Options" };
  desc.add_options()
    ("help,h", "print help information only")
    ("input-file,i", po::value<std::string>(&in_file_name),
      "specify input file name (in OFF format)");
  po::variables_map vm;

  try {
    po::store(po::parse_command_line(argc, argv, desc), vm);
    po::notify(vm);
  } catch (const po::error& e) {
    std::cerr << e.what() << '\n';
    return -1;
  }

  if (vm.count("help")) {
    std::cout << desc << '\n';
    return 1;
  }

  int width = 0;
  int height = 0;
  int prec = 0;
  int num_comps = 0;
  Tri tri;
  if (!in_file_name.empty()) {
    std::ifstream in(in_file_name);
    if (in.is_open()) {
      // Read header
      if (!(in >> width >> height >> num_comps >> prec)) {
        std::cerr << "cannot read mesh header.\n";
        return 1;
      }
      // Read OFF data
      if (input_off_model(in, tri)) {
        std::cerr << "cannot read OFF data.\n";
        return 1;
      }
      in.close();
    }
  } else {  // read file from standard input
    // Read header
    if (!(std::cin >> width >> height >> num_comps >> prec)) {
      std::cerr << "cannot read mesh header.\n";
      return 1;
    }
    // Read OFF data
    if (input_off_model(std::cin, tri)) {
      std::cerr << "cannot read OFF data.\n";
      return 1;
    }
  }

  std::cout << tri.number_of_vertices() << "\n";
  return 0;
}
