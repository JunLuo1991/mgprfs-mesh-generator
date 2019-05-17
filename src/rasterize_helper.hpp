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


#ifndef rasterize_helper_hpp
#define rasterize_helper_hpp

#include "pixel.hpp"
#include <SPL/rasterize.hpp>
#include <array>
#include <boost/container/static_vector.hpp>

namespace rasterize {

/*!
@brief
A helper wrapper class for triangle linear interpolation with
multi-component function values.

@tparam I (\c Int)  The interpolated function value (scan-converted value) type.
@tparam R (\c Real)  The function value type, which must be a real type.
@tparam F (\c Round_operator)  A functor/function to do rounding operation.
@tparam N  The number of components.

@details
This class is used to compute a linear interpolant during the scan
conversion of a triangle with multi-component function values.
The class is intended to be used to obtain an input iterator for
the Triangle_scan_line class in the SPLEL library.

The type Int is the interpolated function value type, it represents the
type of the value that passed to scan_buffer and statitics_collector.
This type can be any type that behaves like an integer or a real number.
The Rational class in SPLEL library is also supported.

The type Real is the type used internally by the class for performing
the numerical calculation for interpolation.
This type can be any built-in real type (e.g. float, double, long double)
The Rational class in SPLEL library is also supported.

NOTE: Type Real should behave like a real number. It should not behave 
like an integer type since this would cause extremely roundoff error.

NOTE: It must be possible to convert from type the Int type to Real type.
A conversion from float/double to Rational is not allowed.

The \c Round_operator has the following signature:
<ul>
       <li>Int round_operator(Real x)
<\ul>
The Round_operator is used to convert the real value produced by
the internal approximating function from type Real to Int.
User can choose not to do rounding here in the Round_operator if
needed, by specifying Int and Real with the same type and directly
return the value x in round_operator.

*/
template <class I, class R, class F, std::size_t N>
class Pixel_triangle_linear_interpolator
{
public:
  using Int = I;
  using Real = R;
  using Round_operator = F;
  using Pixel_value = pixel<Int, N>;
  using Interpolator = SPL::rasterize::Triangle_linear_interpolator<int,
                       Int, Real, Round_operator>;

  /*!
  @brief  Create a pixel triangle linear interpolant.
  @param a_x,a_y  The x and y coordinates of the first vertex.
  @param a_zs  The mult-component function values of the first vertex.
  @param b_x,b_y  The x and y coordinates of the second vertex.
  @param b_zs  The mult-component function values of the second vertex.
  @param c_x,c_y  The x and y coordinates of the third vertex.
  @param c_zs  The mult-component function values of the third vertex.
  @param round_operator  Functor/function to convert from Real to Int.
  */
  Pixel_triangle_linear_interpolator(int a_x, int a_y, const Pixel_value& a_zs,
                                     int b_x, int b_y, const Pixel_value& b_zs,
                                     int c_x, int c_y, const Pixel_value& c_zs,
                                     Round_operator round_operator);

  class Iterator
  {
  public:
    using Interpolator_element_iterator = typename Interpolator::X_iterator;

    Iterator(const boost::container::static_vector<Interpolator, N>& m, int x,
             int y, unsigned mask)
    {
      for (int i = 0; i < N; ++i) {
        per_component_iterator_.push_back(m[i].row_begin_border(x, y, mask));
      }
    }

    Iterator(const boost::container::static_vector<Interpolator, N>& m, int x,
             int y)
    {
      for (int i = 0; i < N; ++i) {
        per_component_iterator_.push_back(m[i].row_begin(x, y));
      }
    }

    Iterator(const Iterator&) = default;
    Iterator& operator=(const Iterator&) = default;

    Iterator& operator++() {
      for (auto&& i : per_component_iterator_) {
        ++i;
      }
      return *this;
    }

    Pixel_value operator*() const {
      Pixel_value values;
      for (int i = 0; i < N; ++i) {
        values[i] = *per_component_iterator_[i];
      }
      return values;
    }

  private:
    boost::container::static_vector<Interpolator_element_iterator, N>
      per_component_iterator_;
  };

  Iterator row_begin_border(int x, int y, unsigned mask) const {
    return Iterator(per_component_interpolator_, x, y, mask);
  }

  Iterator row_begin(int x, int y) const {
    return Iterator(per_component_interpolator_, x, y);
  }

private:
  boost::container::static_vector<Interpolator, N> per_component_interpolator_;
};

template <class I, class R, class F, std::size_t N>
Pixel_triangle_linear_interpolator<I, R, F, N>::Pixel_triangle_linear_interpolator(
                                         int a_x, int a_y, const Pixel_value& a_zs,
                                         int b_x, int b_y, const Pixel_value& b_zs,
                                         int c_x, int c_y, const Pixel_value& c_zs,
                                         Round_operator round_operator)
{
  for (int i = 0; i < N; ++i) {
    per_component_interpolator_.push_back(Interpolator(
      a_x, a_y, a_zs[i], b_x, b_y, b_zs[i], c_x, c_y, c_zs[i], round_operator));
  }
}

/*!
@brief
A helper class for making statistics collecting functor for triangle
scan conversion.

@tparam T1 (\c Update_stat)
A functor/function to update the statistics for each scan-converted value.
@tparam T2 (\c Get_func_iter)
A functor/function returning an input iterator that can be used to
access the original/current scan buffer contents for a horizontal
scan line.
@tparam I  (\c Int)
The interpolated funciton value (scan-converted value) type.
NOTE: This type should be consistent with the type Int used in class
Pixel_triangle_linear_interpolator.
@tparam N  Number of components.

@details
This functor class is intended to be used to obtain a get_stat_iter
functor suitable for use with the \ref Triangle_scan_line class in
SPLEL library.

The Update_stat functor/function has the signature:
<ul>
        <li>void update_stat(const pixel<Int, N>& value,
                             const pixel<int, N>& old_value)
</ul>
The Update_stat functor/function is used to update the statistics
being collected on the scan-converted data. It is called for each
value that is scan converted.
The parameter value is the new scan-converted value on each pixel.
Its pixel value type is Int.
The parameter old_value is the corresponding old value read using
the iterator obtained via Get_func_iter. Its pixel value type is integer.

The Get_func_iter has the signature:
<ul>
        <li>Stat_input_iterator get_func_iter(int x, int y)
</ul>
It returns an input iterator that can be used to read Pixel value of the
image at position (x, y)
*/

template <class T1, class T2, class I, std::size_t N>
class Stat_collector
{
public:
  using Update_stat = T1;
  using Get_func_iter = T2;
  using Int = I;
  using Func_iter = decltype((*((Get_func_iter*)nullptr))(0, 0));
  using Pixel_value = pixel<Int, N>;

  Stat_collector(Update_stat update_stat, Get_func_iter get_func_iter)
    : update_stat_(update_stat)
    , get_func_iter_(get_func_iter)
  {
  }

  class Stat_iter
  {
  public:
    Stat_iter(Update_stat update_stat, Func_iter func_iter)
      : update_stat_(update_stat)
      , func_iter_(func_iter)
    {
    }
    Stat_iter(const Stat_iter&) = default;
    Stat_iter& operator=(const Stat_iter&) = default;
    Stat_iter& operator++() {
      ++func_iter_;
      return *this;
    }
    Stat_iter& operator*() {
      return *this;
    }

    void operator=(const Pixel_value& per_comp_value) {
      update_stat_(per_comp_value, *func_iter_);
    }

  private:
    Update_stat update_stat_;
    Func_iter func_iter_;
  };

  Stat_iter get_iter(int x, int y) const {
    return Stat_iter(update_stat_, get_func_iter_(x, y));
  }

  Stat_iter operator()(int x, int y) const {
    return get_iter(x, y);
  }

private:
  Update_stat update_stat_;
  Get_func_iter get_func_iter_;
};

/*!
@brief A template input iterator class used to read the Pixel value
from an Image.
@tparam N Number of components.

@details
This class is used to provide an input iterator for the Get_func_iter
in Stat_collector. It can be used to read Pixel value when dereferenced.
*/
template <std::size_t N>
class Stat_input_iterator
{
public:
  using Component_element_const_iterator = SPL::Array2<int>::ConstXIterator;
  using Pixel_value = pixel<int, N>;

  /*!
  @brief Constructor
  @param m  The input Image.
  @param x  The x value.
  @param y  The y value.
  @pre m.num_components() = N
  */
  Stat_input_iterator(const Image& m, int x, int y)
  {
    assert(m.num_components() == N);
    for (int i = 0; i < N; ++i) {
      per_component_const_iterator_[i] = m[i].rowBegin(y) + x;
    }
  }

  Stat_input_iterator& operator++()
  {
    for (auto&& i : per_component_const_iterator_) {
      ++i;
    }
    return *this;
  }
  Pixel_value operator*() const
  {
    Pixel_value values;
    for (int i = 0; i < N; ++i) {
      values[i] = *per_component_const_iterator_[i];
    }
    return values;
  }

private:
  std::array<Component_element_const_iterator, N> per_component_const_iterator_;
};

/*!
@brief
An output iterator class that can be used to assign Pixel value
into an image when dereferenced.

@tparam I  (\c Int)
The interpolated function value (scan-converted value) type.
NOTE: This type should be consistent with the type Int used in class
Pixel_triangle_linear_interpolator.
@tparam F (\c Real_to_int)
A functor/functin converting from type Int to integer.
@tparam N Number of components.

@details
This class is used to provide an output iterator for writing data into
the scan buffer. It can be assigned Pixel value when derefernced.

The Real_to_int functor/function has the signature:
<ul>
       <li>int real_to_int(Int x)
<\ul>
The Real_to_int functor/function is used to convert the scan-converted
value from Int to integer. The converted value will be write into the
scan buffer.

*/
template <class I, class F, std::size_t N>
class Scan_buf_output_iterator
{
public:
  using Component_element_iterator = SPL::Array2<int>::XIterator;
  using Int = I;
  using Real_to_int = F;
  using Pixel_value = pixel<Int, N>;

  /*!
  @brief Constructor.
  @param m  The target Image to be written.
  @param x  The x value.
  @param y  The y value.
  @param real_to_int  The Real_to_int functor/function.
  @pre  m.num_components() = N
  */
  Scan_buf_output_iterator(Image& m, int x, int y, Real_to_int real_to_int)
    : real_to_int_(real_to_int)
  {
    assert(m.num_components() == N);
    for (int i = 0; i < N; ++i) {
      per_component_iterator_[i] = m[i].rowBegin(y) + x;
    }
  }

  Scan_buf_output_iterator& operator++()
  {
    for (auto&& i : per_component_iterator_) {
      ++i;
    }
    return *this;
  }

  Scan_buf_output_iterator& operator*() { return *this; }

  void operator=(const Pixel_value& value)
  {
    for (int i = 0; i < N; ++i) {
      *per_component_iterator_[i] = real_to_int_(value[i]);
    }
  }

private:
  std::array<Component_element_iterator, N> per_component_iterator_;
  Real_to_int real_to_int_;
};
};

#endif
