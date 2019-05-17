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

#include "error_diffusion.hpp"
#include "util.hpp"
#include <SPL/Sequence1.hpp>
#include <SPL/Sequence2.hpp>
#include <SPL/math.hpp>
#include <cmath>

void perform_fs_error_diffusion_by_threshold(
  const SPL::Array2<double>& density, double threshold,
  const Fs_error_diffusion_options& options,
  const std::vector<double>& initial_errors,
  std::vector<double>& last_row_errors, SPL::Array2<int>& result)
{

  const int width = density.getWidth();
  const int height = density.getHeight();
  assert(width > 0 && height > 0);

  // copy the initial density function
  SPL::Array2<double> fmap = density;
  // initialize the quantization errors
  SPL::Array2<double> e(width, height, 0.0);

  last_row_errors = std::vector<double>(width, 0.0);
  result = SPL::Array2<int>(width, height);

  // diffuse the initial errors
  if (!initial_errors.empty()) {
    assert(initial_errors.size() == width);
    for (int x = 0; x < width; ++x) {
      fmap(x, 0) += initial_errors[x];
    }
  }

  // for each line of the image from bottom to top...
  for (int y = 0; y < height; ++y) {
    int step_factor = 0;
    int start_pos = 0;
    // if use raster scan order
    if (options.scan_order_ == Scan_order::raster) {
      start_pos = 0;
      step_factor = 1;

      // if use serpentine scan order
    } else if (options.scan_order_ == Scan_order::serpentine) {
      if (y % 2) {
        start_pos = width - 1;
        step_factor = -1;
      } else {
        start_pos = 0;
        step_factor = 1;
      }
    }
    assert((step_factor == 1) || (step_factor == -1));

    // for each point in the line from position start_pos to the other border...
    for (int x = start_pos; (x >= 0) && (x < width); x = x + step_factor) {
      // set result(x, y)
      result(x, y) = (fmap(x, y) >= threshold ? 1 : 0);
      // compute quantization error e(x, y)
      e(x, y) = fmap(x, y) - 2.0 * threshold * result(x, y);

      // the error diffusion weights
      std::vector<double> w = { 7.0 / 16.0, 3.0 / 16.0, 5.0 / 16.0,
                                1.0 / 16.0 };
      // if no_leaky mode, adapt the weights on boundary
      if (!options.is_leaky_) {
        if (y + 1 >= height) {
          w[1] = 0;
          w[2] = 0;
          w[3] = 0;
        }
        if ((x + step_factor < 0) || (x + step_factor >= width)) {
          w[0] = 0;
          w[3] = 0;
        }
        if ((x - step_factor < 0) || (x - step_factor >= width)) {
          w[1] = 0;
        }
        double sum = std::accumulate(w.begin(), w.end(), 0.0);
        if (sum != 0) {
          for (auto& a : w) {
            a /= sum;
          }
        }
      }

      // diffuse errors to position (x + step_factor, y)
      if ((x + step_factor >= 0) && (x + step_factor < width)) {
        fmap(x + step_factor, y) += w[0] * e(x, y);
      }
      if (y + 1 < height) {
        // diffuse errors to position (x - step_factor, y + 1)
        if ((x - step_factor >= 0) && (x - step_factor < width)) {
          fmap(x - step_factor, y + 1) += w[1] * e(x, y);
        }

        // diffuse errors to position (x, y + 1)
        fmap(x, y + 1) += w[2] * e(x, y);

        // diffuse errors to position (x + step_factor, y + 1)
        if ((x + step_factor >= 0) && (x + step_factor < width)) {
          fmap(x + step_factor, y + 1) += w[3] * e(x, y);
        }
      }

      // store the errors that diffused to the last row
      if (y == height - 2) {
        if ((x - step_factor >= 0) && (x - step_factor < width)) {
          last_row_errors[x - step_factor] += w[1] * e(x, y);
        }

        last_row_errors[x] += w[2] * e(x, y);

        if ((x + step_factor >= 0) && (x + step_factor < width)) {
          last_row_errors[x + step_factor] += w[3] * e(x, y);
        }
      }
    }
  }
}

int perform_fs_error_diffusion_by_size(const SPL::Array2<double>& density,
                                   int num_points, int tolerance,
                                   const Fs_error_diffusion_options& options,
                                   SPL::Array2<int>& result)
{

  const int width = density.getWidth();
  const int height = density.getHeight();
  assert(width > 0 && height > 0);
  assert(num_points >= 0);

  if (num_points == 0) {
    result = SPL::Array2<int>(width, height, 0);
    return 0;
  }

  // compute the acceptable range of resulting number of points in the bitmap
  tolerance = std::abs(tolerance);
  const int lower_bound = num_points - tolerance;
  const int upper_bound = num_points + tolerance;
  std::vector<double> startup_errors;  //the startup errors for error diffusion
  std::vector<double> output_errors;

  const SPL::Array2<double>& d = density;
  double d_sum = std::accumulate(d.begin(), d.end(), 0.0);
  // WARNING: if d_sum is zero (in other words, all pixel value are zero)
  // the result mesh will have no point been selected. We do it here to
  // avoid dead loop when trying to find a good start t_low
  if (d_sum == 0) {
    std::cerr << "Can not run error diffusion by size due to all density "
                 "values are zero\n";
    return -1;
  }

  // set the threshold
  double t = d_sum / (2.0 * num_points);


  // get the startup errors
  if(get_fs_error_diffusion_startup_errors(d, num_points,
       options, startup_errors)) {
    return -1;
  }

  // perform Floyd-steinberg error diffusion with the initial threshold
  perform_fs_error_diffusion_by_threshold(d, t, options, startup_errors,
                                          output_errors, result);

  // count the number of points selected in the result bitmap
  int size = std::count(result.begin(), result.end(), 1);

  double t_high = 0;
  double t_low = 0;
  // find the range of threshold t_high and t_low to be used for binary search
  if (size > upper_bound) {
    while (size > upper_bound) {
      t_low = t;
      t *= 2;
      perform_fs_error_diffusion_by_threshold(d, t, options, startup_errors,
                                              output_errors, result);
      size = std::count(result.begin(), result.end(), 1);
    }
    t_high = t;
  } else if (size < lower_bound) {
    while (size < lower_bound) {
      t_high = t;
      t /= 2;
      perform_fs_error_diffusion_by_threshold(d, t, options, startup_errors,
                                              output_errors, result);
      size = std::count(result.begin(), result.end(), 1);
    }
    t_low = t;
  } else {
    t_high = t;
    t_low = t;
  }

  // perform binary search to find the optimal threshold t
  // and get the result bitmap
  while (t_low <= t_high) {

    t = (t_high + t_low) / 2.0;
    perform_fs_error_diffusion_by_threshold(d, t, options, startup_errors,
                                            output_errors, result);
    size = std::count(result.begin(), result.end(), 1);

    if (size > upper_bound) {
      t_low = t;
    } else if (size < lower_bound) {
      t_high = t;
    } else {
      break;
    }

    // if same threshold will be used in the next iteration
    // this indicates that the target size is not achievable.
    if ((t_low + t_high) / 2.0 == t) {
      std::cerr << "The target size is not achievable by performing error diffusion\n";
      return -1;
    }
  }

  return 0;
}


int get_fs_error_diffusion_startup_errors(const SPL::Array2<double>& density,
                     int num_points, const Fs_error_diffusion_options& options,
                     std::vector<double>& startup_errors)
{
  
  const int width = density.getWidth();
  const int height = density.getHeight();
  assert(width > 0 && height > 0);
  assert(num_points >= 0);

  double d_sum = std::accumulate(density.begin(), density.end(), 0.0);
  double t = d_sum / (2.0 * num_points);  // threshold
  startup_errors = std::vector<double>(width);

  switch (options.initial_error_mode_) {
    case Initial_error_mode::zero:
      startup_errors = std::vector<double>(width, 0.0);
      break;
    case Initial_error_mode::random:
      // random generate some errors in the range -2t to 1.0 - 2t
      for (int x = 0; x < width; ++x) {
        double r = generate_random_number<double>(-2.0 * t, 1.0 - 2.0 * t, 0);
        startup_errors[x] = r;
      }
      break;
    case Initial_error_mode::mirror: {
        std::vector<double> output_errors;
  	SPL::Array2<int> bitmask;
        SPL::Array2<double> mirror_d = density;

        // flip the density function to get a mirrored version
        mirror_d.flipud();

        // perform error diffusion using mirror_d to get the errors diffused to the last row
        perform_fs_error_diffusion_by_threshold(mirror_d, t, options, startup_errors,
                                                output_errors, bitmask);

        startup_errors = output_errors;
      }
      break;
    default:
      std::cerr << "Invalid initial error mode argument\n";
      return -1;
  }
  return 0;
}
