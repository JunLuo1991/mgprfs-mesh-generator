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



// This program reads an image in PNM format from standard input
// stream and output its mirrored image (i.e., vertical flipped image)
// to standard output stream.

#include <iostream>
#include "image.hpp"
#include "util.hpp"

int main(int argc, char ** argv )
{
  Image input_image;
  // Read the image from standard input stream
  if (auto status = input_image.input(std::cin); status) {
    return status;
  }

  int prec = input_image.precision();
  int width = input_image.width();
  int height = input_image.height();
  int num_comps = input_image.num_components();

  // Construct the mirrored image
  Image mirror(width, height, prec, num_comps);
  for(int y = height - 1; y >= 0; --y) {
    for(int x = 0; x < width; ++x) {
       for(int k = 0; k < num_comps; ++k) {
          mirror[k](x, height-1-y) = input_image[k](x, y);
       }
    }
  }

  // Output the mirrored image to standard output stream
  mirror.output(std::cout);
  std::cerr << "output mirror image complete!\n";
  return 0;
}
