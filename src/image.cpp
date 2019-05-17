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

#include "image.hpp"
#include "util.hpp"
#include <cmath>
#include <iostream>

Image::Image()
  : width_(0)
  , height_(0)
  , prec_(0)
{
}

Image::Image(int width, int height, int prec, int num_comps)
  : width_(width)
  , height_(height)
  , prec_(prec)
{
  assert(num_comps == 1 || num_comps == 3);
  for (int i = 0; i < num_comps; ++i) {
    components_.push_back(SPL::Array2<int>(width, height));
  }
}

Image::Image(std::istream& in)
{
  if (input(in)) {
    throw std::runtime_error("failed to construct an image.");
  }
}

int Image::width() const {
  return width_;
}

int Image::height() const {
  return height_;
}

int Image::precision() const {
  return prec_;
}

int Image::num_components() const{
  return components_.size();
}

Image::Component_iterator Image::components_begin() {
  return components_.begin();
}

Image::Component_const_iterator Image::components_begin() const {
  return components_.begin();
}

Image::Component_iterator Image::components_end() {
  return components_.end();
}

Image::Component_const_iterator Image::components_end() const {
  return components_.end();
}

const Image::Component& Image::operator[](int i) const {
  return components_[i];
}

Image::Component& Image::operator[](int i) {
  return components_[i];
}

int Image::input(std::istream& in)
{
  int max_val = 0;
  bool sgnd = false;
  if (SPL::decodePnm(in, components_, max_val, sgnd)) {
    return -1;
  }

  assert(components_.size() == 1 || components_.size() == 3);

  width_ = components_begin()->getWidth();
  height_ = components_begin()->getHeight();
  prec_ = value_to_precision(max_val);

  return 0;
}

int Image::output(std::ostream& out) const
{
  int max_value = precision_to_value(prec_);
  return SPL::encodePnm(out, components_, max_value, false);
}

void rgb_to_gray(const Image& rgb_image, Image& gray_image)
{
  assert(rgb_image.num_components() == 3);

  const int width = rgb_image.width();
  const int height = rgb_image.height();
  const int prec = rgb_image.precision();
  const int z_max = precision_to_value(prec);
  gray_image = Image(width, height, prec, 1);

  // for each pixel in the image...
  for (int y = 0; y < height; ++y) {
    for (int x = 0; x < width; ++x) {
      pixel<double, 3> val{ static_cast<double>(rgb_image[0](x, y)),
                            static_cast<double>(rgb_image[1](x, y)),
                            static_cast<double>(rgb_image[2](x, y)) };
      // convert from rgb to luma
      int gray_val = static_cast<int>(rgb_to_luma<double>(val));
      gray_image[0](x, y) = SPL::clip<int>(gray_val, 0, z_max);
    }
  }
}
