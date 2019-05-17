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

#ifndef pixel_hpp
#define pixel_hpp

#include <SPL/rasterize.hpp>
#include <algorithm>
#include <boost/rational.hpp>
#include <cassert>
#include <cstddef>
#include <initializer_list>
#include <iostream>

/*
@brief A variable template (since C++14) to make a ratio.
@tparam T  The type of the variable.
@tparam numerator  The value of numerator.
@tparam denominator  The value of denominator.
@return  The ratio result.
*/
template <class T, int numerator, int denominator>
constexpr T make_ratio = T(numerator) / T(denominator);

/*!
@brief  A class template representing an image pixel.
@tparam S  The pixel value type.
@tparam N  The number of components in the pixel.
@details
This is a helper class for the Image class.
It can be used to represent a single pixel value.
*/
template <class S, std::size_t N>
class pixel
{
public:
  using value_type = S;
  using iterator = value_type*;
  using const_iterator = const value_type*;

  /*
  @brief  Create an empty pixel.
  */
  pixel();

  /*
  @brief  Create a pixel using initializer list
  */
  pixel(std::initializer_list<value_type>);

  /*
  @brief  Get the number of components
  */
  static constexpr std::size_t num_components();

  /*
  @brief
  Get the const reference to the component at position i in the pixel.
  @param i  Position of a component in the pixel.
  @return  The const reference to the component at the specified position.
  in the pixel.
  */
  const value_type& operator[](int i) const;

  /*
  @brief
  Get the mutable reference to the component at position i in the pixel.
  @param i  Position of a component in the pixel.
  @return  The mutable reference to the component at the specified position
  in the pixel.
  */
  value_type& operator[](int i);

  /*
  @brief
  Get the mutable iterator for the first component in the pixel.
  */
  iterator begin();

  /*
  @brief
  Get the mutable iterator for one past the last component in the pixel.
  */
  iterator end();

  /*
  @brief
  Get the const iterator for the first component in the pixel.
  */
  const_iterator begin() const;

  /*
  @brief
  Get the const iterator for the one past the last component in the pixel. 
  */
  const_iterator end() const;

private:

  //! The pixel data
  value_type values_[num_components()];
};

/*
@brief  Output a pixel to the specified stream.
@detail
An extra empty space will be output between components in the pixel.
*/
template <class S, std::size_t N>
std::ostream& operator<<(std::ostream& out, const pixel<S, N>& p);

/*
@brief  Input a pixel from the specified stream.
*/
template <class S, std::size_t N>
std::istream& operator>>(std::istream& in, pixel<S, N>& p);

/*!
@brief A function template to convert a RGB pixel to luma space value.
@tparam T  The pixel value type.
@param rgb  The rgb pixel value.
@return  The converted luma color space value.

@detail
Convert the rgb value to luma color space value.
Note: The returned luma value is not rounded in this function,
you may need to round it outside the function body.
*/
template <class T>
inline constexpr T rgb_to_luma(const pixel<T, 3>& rgb);

/*!
@brief
A function template to convert a pixel from RGB to YCbCr color space.
@tparam T  The pixel value type.
@param rgb  The input RGB pixel value.
@param ycrcb  The output YCbCr pixel value.

@detail
Convert from RGB to full-range YCrCb color space.
Refered from http://www.equasys.de/colorconversion.html
Note: The value ycbcr is not rounded in this function,
you may need to round it outside the function body.
*/
template <typename T>
inline constexpr void 
rgb_to_ycbcr(const pixel<T, 3>& rgb, pixel<T, 3>& ycbcr);



/************************************************************

The following shows the definition of above function templates
and member functions of class template

************************************************************/

template <class S, std::size_t N>
inline pixel<S, N>::pixel()
{
}

template <class S, std::size_t N>
inline pixel<S, N>::pixel(typename std::initializer_list<S> list)
{
  assert(list.size() == num_components());
  std::copy_n(list.begin(), num_components(), values_);
}

template <class S, std::size_t N>
inline constexpr std::size_t pixel<S, N>::num_components()
{
  return N;
}

template <class S, std::size_t N>
inline const typename pixel<S, N>::value_type&
pixel<S, N>::operator[](int i) const
{
  return values_[i];
}

template <class S, std::size_t N>
inline typename pixel<S, N>::value_type& pixel<S, N>::operator[](int i)
{
  return values_[i];
}

template <class S, std::size_t N>
inline typename pixel<S, N>::iterator pixel<S, N>::begin()
{
  return &values_[0];
}

template <class S, std::size_t N>
inline typename pixel<S, N>::iterator pixel<S, N>::end()
{
  return &values_[num_components()];
}

template <class S, std::size_t N>
inline typename pixel<S, N>::const_iterator pixel<S, N>::begin() const
{
  return &values_[0];
}

template <class S, std::size_t N>
inline typename pixel<S, N>::const_iterator pixel<S, N>::end() const
{
  return &values_[num_components()];
}

template <class S, std::size_t N>
std::ostream& operator<<(std::ostream& out, const pixel<S, N>& p)
{
  for (auto i = p.begin(); i != p.end(); ++i) {
    if (i != p.begin()) {
      out << ' ';
    }
    out << *i;
  }
  return out;
}

template <class S, std::size_t N>
std::istream& operator>>(std::istream& in, pixel<S, N>& p)
{
  for (auto i = p.begin(); i != p.end(); ++i) {
    in >> *i;
  }
  return in;
}

template <class T>
inline constexpr T rgb_to_luma(const pixel<T, 3>& rgb)
{
  return make_ratio<T, 299, 1000> * rgb[0] + make_ratio<T, 587, 1000> * rgb[1] +
         make_ratio<T, 114, 1000> * rgb[2];
}

template <typename T>
inline constexpr void rgb_to_ycbcr(const pixel<T, 3>& rgb, pixel<T, 3>& ycbcr)
{
  ycbcr[0] = make_ratio<T, 299, 1000> * rgb[0] +
             make_ratio<T, 587, 1000> * rgb[1] +
             make_ratio<T, 114, 1000> * rgb[2];
  ycbcr[1] = make_ratio<T, 128, 1> + make_ratio<T, -169, 1000> * rgb[0] +
             make_ratio<T, -331, 1000> * rgb[1] +
             make_ratio<T, 500, 1000> * rgb[2];
  ycbcr[2] = make_ratio<T, 128, 1> + make_ratio<T, 500, 1000> * rgb[0] +
             make_ratio<T, -419, 1000> * rgb[1] +
             make_ratio<T, -81, 1000> * rgb[2];
}

#endif
