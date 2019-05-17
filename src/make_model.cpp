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

#include "mesh_generator.hpp"
#include <iostream>

int main(int argc, char** argv)
{
  try {
    Mesh_generator mesh_generator;
    if (auto status = mesh_generator(std::cin, std::cout, argc, argv); status) {
      std::cerr << "mesh generation failed with status " << status << '\n';
      return 1;
    }
  } catch(const std::exception& e) {
    std::cerr << "mesh generation failed due to exception: " << e.what() << '\n';
    return -1;
  }
  return 0;
}
