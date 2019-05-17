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

#ifndef image_hpp
#define image_hpp

#include <SPL/Array2.hpp>
#include <algorithm>
#include <iostream>
#include <random>
#include "pixel.hpp"

/*!
@breif
A class representing an image.
@detail
The image can contain one or three components (for grayscale
or rgb-color images respectively).
*/
class Image
{
public:
  //! The image component type.
  using Component = SPL::Array2<int>;

  //! The image component iterator type
  using Component_iterator = std::vector<Component>::iterator;

  //! The image component const iterator type
  using Component_const_iterator = std::vector<Component>::const_iterator;

  /*
   @brief  Create an empty image
   @detail  The image created has no components.
   */
  Image();

  /*
  @brief
  Create an empty image with given size, precision, and number of components.
  @param width  The width of the image.
  @param height  The height of the image.
  @param prec  The precision of the image.
  @param num_comps  The number of components of the image.
  */
  Image(int width, int height, int prec, int num_comps);

  /*
  @brief
  Create an image by reading the image data in PNM format from a stream.
  @param in  The stream to read from.
  */
  Image(std::istream& in);

  /*
  @brief  The move constructor.
  */
  Image(Image&&) = default;

  /*
  @brief  The copy constructor.
  */
  Image(const Image&) = default;

  /*
  @brief  The move assignment operator.
  */
  Image& operator=(Image&&) = default;

  /*
  @brief  The copy assignment operator.
  */
  Image& operator=(const Image&) = default;

  /*
  @brief  Destroy an image.
  */
  ~Image() = default;

  /*
  @brief  Get the width of the image components.
  */
  int width() const;

  /*
  @brief  Get the height of the image components.
  */
  int height() const;

  /*
  @brief  Get the sample precision of the image components.
  */
  int precision() const;

  /*
  @brief  Get the number of image components.
  */
  int num_components() const;

  /*
  @brief
  Get the mutable iterator for the first component in the image.
  */
  Component_iterator components_begin();

  /*
  @brief
  Get the const iterator for the first component in the image.
  */
  Component_const_iterator components_begin() const;

  /*
  @brief
  Get the mutable iterator for one past the last component in the image.
  */
  Component_iterator components_end();

  /*
  @breif
  Get the const iterator for the one past the last component in the image.
  */
  Component_const_iterator components_end() const;

  /*
  @breif  Get a const reference to a particular image component.
  @param i  The index of the component
  */
  const Component& operator[](int i) const;

  /*
   @brief  Get a mutable reference to a particular image component
   @param i  The index of the component
   */
  Component& operator[](int i);

  /*
   @breif  Read an image from a stream in PNM format.
   @param in  The stream to read from.
   @return
   Upon success, zero is returned; otherwise, a non-zero value is returned.
   */
  int input(std::istream& in);

  /*
   @brief  Write an image to a stream in PNM format.
   @param out  The stream to write to.
   @return
   Upon success, zero is returned; otherwise, a non-zero value is returned.
   */
  int output(std::ostream& out) const;

private:
  //! The width of the image components.
  int width_;

  //! The height of the image components.
  int height_;

  //! The sample precision (i.e., bits/sample) of the image components.
  int prec_;

  //! The image component data.
  std::vector<Component> components_;
};

/*
@brief  Convert a RGB image to a grayscale image.
@param rgb_image  The input RGB image.
@param gray_image  The converted grayscale image.
*/
void rgb_to_gray(const Image& rgb_image, Image& gray_image);

#endif
