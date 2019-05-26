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



// This program reads an image in PNM format from standard input stream
// and output its basic information (i.e., width, height, number of components,
// precision and max value) to standard output stream.

#include <iostream>
#include "image.hpp"
#include "util.hpp"

int main(int argc, char ** argv )
{
  Image input_image;
  // read the image from standard input stream
  if (auto status = input_image.input(std::cin); status) {
    return status;
  }

  int prec = input_image.precision();
  int max_val = precision_to_value(prec);

  // Output the image information
  std::cout << input_image.width() << "\n"
            << input_image.height() << "\n"
            << input_image.num_components() << "\n"
            << input_image.precision() << "\n"
            << max_val << "\n";
  return 0;

}
