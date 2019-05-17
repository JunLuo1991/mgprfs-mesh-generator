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
@file off_to_triangulation.cpp
@brief
This program reads a triangle mesh model in MODEL (i.e., header + OFF) format
and output the triangulation file in TRI fomrat.
*/

#include <CGAL/Cartesian.h>
#include <CGAL/bounding_box.h>
#include <SPL/Array2.hpp>
#include <SPL/cgal/Delaunay_triangulation_2.hpp>
#include <SPL/cgal/Triangulation_hierarchy_2.hpp>
#include <SPL/rasterize.hpp>
#include <array>
#include <boost/program_options.hpp>
#include <iostream>
#include <string>

#include "pixel.hpp"

class Tri_vertex;
struct Tri_traits
{
    using Geom_kernel = CGAL::Filtered_kernel<CGAL::Cartesian<double>>;
    using Vertex = Tri_vertex;
    using Edge = SPL::Triangulation_edge_base_2<Tri_traits>;
    using Face = SPL::Triangulation_face_base_2<Tri_traits>;
    using Halfedge = SPL::Triangulation_halfedge_base_2<Tri_traits>;
};

class Tri_vertex : public SPL::Triangulation_hierarchy_vertex_base_2<Tri_traits>
{
public:
    Tri_vertex() : z(0){};
    int z;
    pixel<int, 3> color;
};

#if defined(USE_TRIANGULATION_HIERARCHY)
  using Tri =
    SPL::Triangulation_hierarchy_2<SPL::Delaunay_triangulation_2<Tri_traits>>;
#else
  using Tri = SPL::Delaunay_triangulation_2<Tri_traits>;
#endif


int input_off_model(std::istream& in, Tri& tri)
{

  auto set_vertex_z = [](Tri::Vertex_handle v, double z) { v->z = z; };
  auto set_vertex_color = [](Tri::Vertex_handle v,
                             const std::vector<double>& color, bool as_int) {
    assert(color.size() == 4);
    // we only need the r,g,b components, ignore the alpha component here
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

  using Kernel = CGAL::Filtered_kernel<CGAL::Cartesian<double>>;
  using Point = Tri::Point;

  Tri tri;
  std::vector<Point> points;
  int width = 0;
  int height = 0;
  int prec = 0;
  int num_comps = 0;

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

  tri.output(std::cout);

  return 0;
}


