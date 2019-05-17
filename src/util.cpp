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

#include "util.hpp"
#include <SPL/Array2.hpp>
#include <SPL/Sequence2.hpp>
#include <SPL/math.hpp>
#include <stdexcept>

int precision_to_value(int prec)
{
  return (1 << prec) - 1;
}

int value_to_precision(int value)
{
  int prec = 1;
  while (value >> 1) {
    ++prec;
    value = value >> 1;
  }
  return prec;
}


void get_binomial_filter(int smooth_order, SPL::RealSequence1& result_filter)
{
  assert(smooth_order % 2 == 1);

  // based on smooth_order get the mask of the binomial filter.
  // The n-order binomial filter correspond to the nth row of Pascal's triangle.
  SPL::RealSequence1::ElemType mask[smooth_order];
  double x = 1;
  int n = smooth_order - 1;
  for (int i = 0; i <= n; ++i) {
    mask[i] = x;
    x = x * (n - i) / (i + 1);
  }

  // get the start index
  int start_index = 0 - (smooth_order - 1) / 2;
  // get the filter
  result_filter = SPL::RealSequence1(start_index, smooth_order, mask);
  result_filter /= result_filter.sum();
}

void get_mean_filter(int smooth_order, SPL::RealSequence1& result_filter)
{
  assert(smooth_order % 2 == 1);
  int start_index = 0 - (smooth_order - 1) / 2;
  result_filter = SPL::RealSequence1(start_index, smooth_order, 1.0);
  result_filter /= result_filter.sum();
}

void get_smooth_filter(int smooth_order, Smooth_operator smooth_operator,
                       SPL::RealSequence1& result_filter)
{
  switch (smooth_operator) {
    case Smooth_operator::binomial:
      get_binomial_filter(smooth_order, result_filter);
      break;
    case Smooth_operator::mean:
      get_mean_filter(smooth_order, result_filter);
      break;
    default:
      throw std::invalid_argument("unsupported smooth operator");
  }
}

void get_derivative_filter(int x_order, int y_order,
  SPL::RealSequence1& horz_filter, SPL::RealSequence1& vert_filter)
{
  SPL::RealSequence1::ElemType mask1[] = { -0.5, 0, 0.5 };
  SPL::RealSequence1::ElemType mask2[] = { 1, -2, 1 };
  SPL::RealSequence1 first_order_filter = SPL::RealSequence1(-1, 3, mask1);
  SPL::RealSequence1 second_order_filter = SPL::RealSequence1(-1, 3, mask2);
  SPL::RealSequence1 delta(0, 1, 1.0);

  if ((x_order == 1 && y_order == 0)) {
    horz_filter = first_order_filter;
    vert_filter = delta;
  } else if (x_order == 0 && y_order == 1) {
    horz_filter = delta;
    vert_filter = first_order_filter;
  } else if (x_order == 1 && y_order == 1) {
    horz_filter = first_order_filter;
    vert_filter = first_order_filter;
  } else if (x_order == 2 && y_order == 0) {
    horz_filter = second_order_filter;
    vert_filter = delta;
  } else if (x_order == 0 && y_order == 2) {
    horz_filter = delta;
    vert_filter = second_order_filter;
  } else {
    std::string err_str = "invalid derivative orders: x_order = " +
      std::to_string(x_order) + ", y_order = " + std::to_string(y_order);
    throw std::invalid_argument(err_str);
  }
}

void get_smoothing_and_derivative_filters(int x_order, int y_order,
                                     const Smooth_options& smooth_options,
                                     SPL::RealSequence1& horz_smooth_filter,
                                     SPL::RealSequence1& vert_smooth_filter,
                                     SPL::RealSequence1& horz_deriv_filter,
                                     SPL::RealSequence1& vert_deriv_filter)
{
  SPL::RealSequence1 smooth_filter;
  SPL::RealSequence1 derivative_filter;
  SPL::RealSequence1 delta(0, 1, 1.0);
  Smooth_direction smooth_direction = smooth_options.smooth_direction_;

  get_derivative_filter(x_order, y_order, horz_deriv_filter, vert_deriv_filter);

  get_smooth_filter(smooth_options.smooth_order_,
                    smooth_options.smooth_operator_, smooth_filter);

  if (smooth_direction == Smooth_direction::both) {
    horz_smooth_filter = smooth_filter;
    vert_smooth_filter = smooth_filter;
  } else if (smooth_direction == Smooth_direction::orthogonal) {
    horz_smooth_filter = (y_order > 0 ? smooth_filter : delta);
    vert_smooth_filter = (x_order > 0 ? smooth_filter : delta);
  }
}

void get_smoothed_image(const SPL::Array2<double>& input,
                        const Smooth_options& smooth_options,
                        SPL::Array2<double>& result)
{
  SPL::RealSequence2 f(input);
  SPL::RealSequence2 result_seq;
  SPL::RealSequence1 smooth_filter;
  get_smooth_filter(smooth_options.smooth_order_,
                    smooth_options.smooth_operator_, smooth_filter);
  result_seq = SPL::convolveSeparable(f, smooth_filter, smooth_filter,
                                      smooth_options.bound_policy_);
  result = result_seq.getArray();
}

void get_smoothed_derivative(const SPL::Array2<double>& input, int x_order,
                             int y_order, int derivative_bound_policy,
                             const Smooth_options& smooth_options,
                             bool use_double_convolution,
                             SPL::Array2<double>& result)
{

  SPL::RealSequence2 f(input);
  SPL::RealSequence2 smoothed_seq;
  SPL::RealSequence2 result_seq;
  SPL::RealSequence1 horz_smooth_filter;
  SPL::RealSequence1 vert_smooth_filter;
  SPL::RealSequence1 horz_deriv_filter;
  SPL::RealSequence1 vert_deriv_filter;

  get_smoothing_and_derivative_filters(x_order, y_order, smooth_options,
     horz_smooth_filter, vert_smooth_filter, horz_deriv_filter, vert_deriv_filter);

  // if chosing smoothing image first and then convolve with derivative filter
  if (use_double_convolution) {
    smoothed_seq = SPL::convolveSeparable(
      f, horz_smooth_filter, vert_smooth_filter, smooth_options.bound_policy_);

    result_seq = SPL::convolveSeparable(
      smoothed_seq, horz_deriv_filter, vert_deriv_filter, derivative_bound_policy);
  } else {
    auto horz_filter = SPL::convolve(horz_smooth_filter, horz_deriv_filter);
    auto vert_filter = SPL::convolve(vert_smooth_filter, vert_deriv_filter);
    result_seq = SPL::convolveSeparable(f, horz_filter, vert_filter,
                                        derivative_bound_policy);
  }

  result = result_seq.getArray();
}

void get_mag_laplacian_grayscale(const SPL::Array2<double>& input,
                             int derivative_bound_policy,
                             const Smooth_options& smooth_options,
                             bool use_double_convolution,
                             SPL::Array2<double>& result)
{

  const int width = input.getWidth();
  const int height = input.getHeight();

  SPL::Array2<double> fxx(width, height);
  SPL::Array2<double> fyy(width, height);
  get_smoothed_derivative(input, 2, 0, derivative_bound_policy, smooth_options,
                          use_double_convolution, fxx);
  get_smoothed_derivative(input, 0, 2, derivative_bound_policy, smooth_options,
                          use_double_convolution, fyy);

  result = SPL::Array2<double>(width, height);
  for (int y = 0; y < height; ++y) {
    for (int x = 0; x < width; ++x) {
      result(x, y) = SPL::absVal(fxx(x, y) + fyy(x, y));
    }
  }
}

void get_mmsodd_grayscale(const SPL::Array2<double>& input,
                          int derivative_bound_policy,
                          const Smooth_options& smooth_options,
                          bool use_double_convolution,
                          SPL::Array2<double>& result)
{
  const int width = input.getWidth();
  const int height = input.getHeight();
  result = SPL::Array2<int>(width, height);

  SPL::Array2<double> fxx(width, height);
  SPL::Array2<double> fxy(width, height);
  SPL::Array2<double> fyy(width, height);
  get_smoothed_derivative(input, 2, 0, derivative_bound_policy, smooth_options,
                          use_double_convolution, fxx);
  get_smoothed_derivative(input, 1, 1, derivative_bound_policy, smooth_options,
                          use_double_convolution, fxy);
  get_smoothed_derivative(input, 0, 2, derivative_bound_policy, smooth_options,
                          use_double_convolution, fyy);

  for (int y = 0; y < height; ++y) {
    for (int x = 0; x < width; ++x) {
      double alpha = (fxx(x, y) + fyy(x, y)) / 2.0;
      double beta =
        std::sqrt(SPL::sqr(fxx(x, y) - fyy(x, y)) / 4.0 + SPL::sqr(fxy(x, y)));
      result(x, y) =
        std::max(SPL::absVal(alpha + beta), SPL::absVal(alpha - beta));
    }
  }

}

void get_mdg_color(const std::array<SPL::Array2<double>, 3>& input,
                   int derivative_bound_policy,
                   const Smooth_options& smooth_options,
                   bool use_double_convolution,
                   SPL::Array2<double>& result)
{
  const int width = input[0].getWidth();
  const int height = input[0].getHeight();
  result = SPL::Array2<int>(width, height);

  SPL::Array2<double> fRx(width, height);
  SPL::Array2<double> fRy(width, height);
  SPL::Array2<double> fGx(width, height);
  SPL::Array2<double> fGy(width, height);
  SPL::Array2<double> fBx(width, height);
  SPL::Array2<double> fBy(width, height);

  get_smoothed_derivative(input[0], 1, 0, derivative_bound_policy,
                          smooth_options, use_double_convolution, fRx);
  get_smoothed_derivative(input[0], 0, 1, derivative_bound_policy,
                          smooth_options, use_double_convolution, fRy);
  get_smoothed_derivative(input[1], 1, 0, derivative_bound_policy,
                          smooth_options, use_double_convolution, fGx);
  get_smoothed_derivative(input[1], 0, 1, derivative_bound_policy,
                          smooth_options, use_double_convolution, fGy);
  get_smoothed_derivative(input[2], 1, 0, derivative_bound_policy,
                          smooth_options, use_double_convolution, fBx);
  get_smoothed_derivative(input[2], 0, 1, derivative_bound_policy,
                          smooth_options, use_double_convolution, fBy);

  for (int y = 0; y < height; ++y) {
    for (int x = 0; x < width; ++x) {
      double gxx =
        SPL::sqr(fRx(x, y)) + SPL::sqr(fGx(x, y)) + SPL::sqr(fBx(x, y));
      double gyy =
        SPL::sqr(fRy(x, y)) + SPL::sqr(fGy(x, y)) + SPL::sqr(fBy(x, y));
      double gxy =
        fRx(x, y) * fRy(x, y) + fGx(x, y) * fGy(x, y) + fBx(x, y) * fBy(x, y);

      double theta = 0.0;
      if (gxx == gyy) {
        theta = pi / 2.0;
      } else {
        theta = 0.5 * std::atan2(2.0 * gxy, gxx - gyy);
      }
      double theta_pi_half = theta + pi / 2.0;
      double f1 = 0.5 * (gxx + gyy + (gxx - gyy) * std::cos(2.0 * theta) +
                  2.0 * gxy * std::sin(2.0 * theta));
      double f2 = 0.5 * (gxx + gyy + (gxx - gyy) * std::cos(2.0 * theta_pi_half)
                  + 2.0 * gxy * std::sin(2.0 * theta_pi_half));
      result(x, y) =
        std::max(std::sqrt(SPL::absVal(f1)), std::sqrt(SPL::absVal(f2)));
    }
  }
}

