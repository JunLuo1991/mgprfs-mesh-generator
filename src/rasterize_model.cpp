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
@file rasterize_model.cpp
@brief
This file implements the rasterize_model program.
It is used to read an OFF mesh model, do linear interpolation to each triangle
face and output the rasterized result as an image file in PNM/PPM/PGM format.
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

#include "image.hpp"
#include "rasterize_helper.hpp"
#include "scan_triangle.hpp"

// The make triangulation data structure
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
int
input_off_model(std::istream& in, Tri& tri)
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

/*!
@brief Scan a triangulation face and write the scan-converted
data into recon_image.
@tparam N number of components in the image.
@param face A Face_const_handle that refers to the scanned face.
@param recon_image The re-constructed image,
*/
template <std::size_t N>
void scan_face(Tri::Face_const_handle face, Image& recon_image)
{
  Tri::Halfedge_const_handle h = face->halfedge();
  std::array<Tri::Vertex_const_handle, 3> vs{
    { h->vertex(), h->next()->vertex(), h->prev()->vertex() }
  };

  std::array<Tri::Point, 3> ps{ { vs[0]->point(), vs[1]->point(),
                                  vs[2]->point() } };

  Scan_stats stats;
  int flags = Scan_triangle_options::write_scan_buf;
  Image ori_image;
  // The error metric doesn't matter which one to select here since
  // we are only trying to get the recon image, doesn't care about
  // the error statistics.
  Error_metric error_metric = Error_metric::mean_comp_se;

  if (N == 1) {
    scan_triangle_grayscale(ps[0].x(), ps[0].y(), ps[1].x(), ps[1].y(),
                            ps[2].x(), ps[2].y(), vs[0]->z, vs[1]->z, vs[2]->z,
                            ori_image, error_metric, recon_image, stats, flags);

  } else if (N == 3) {
    scan_triangle_color(ps[0].x(), ps[0].y(), ps[1].x(), ps[1].y(), ps[2].x(),
                        ps[2].y(), vs[0]->color, vs[1]->color, vs[2]->color,
                        ori_image, error_metric, recon_image, stats, flags);
  } else {
    std::cerr << "invalid number of components\n";
    std::exit(1);
  }
}

int main(int argc, char** argv)
{
  std::string out_file_name;
  std::string in_file_name;

  namespace po = boost::program_options;
  po::options_description desc{ "Options" };
  desc.add_options()
  ("help,h", "print help information only")
  ("output-file,o", po::value<std::string>(&out_file_name),
   "specify output file name (in PNM/PPM/PGM format)")
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
  } else {
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

  Image recon_image(width, height, prec, num_comps);

  // Scan each face and do rasterization...
  for (Tri::Face_const_iterator faceIter = tri.faces_begin();
       faceIter != tri.faces_end(); ++faceIter) {
    if (num_comps == 3) {
      scan_face<3>(&*faceIter, recon_image);
    } else {
      scan_face<1>(&*faceIter, recon_image);
    }
  }

  if (!out_file_name.empty()) {
    std::ofstream out(out_file_name);
    if (out.is_open()) {
      if (recon_image.output(out)) {
        std::cerr << "encode Pnm file failed\n";
        return 1;
      }
      out.close();
    }
  } else {
    if (recon_image.output(std::cout)) {
      std::cerr << "encode Pnm file failed\n";
      return 1;
    }
  }
  return 0;
}
