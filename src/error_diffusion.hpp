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
@file error_diffusion.hpp
@brief This file provides the interface of some enumerations and functions
related to Floyd-steinberg error diffusion method.
*/
#ifndef error_diffusion_hpp
#define error_diffusion_hpp

#include "image.hpp"
#include <SPL/Array2.hpp>
#include <vector>

/*!
 @brief The error diffusion scan order.
 */
enum class Scan_order
{
  raster,
  serpentine,
};

/*!
@brief
The error diffusion startup error mode.
@detail
This specifies the initial errors that will be distributed
to the startup row when doing error diffusion.
*/
enum class Initial_error_mode
{
  zero,
  random,
  mirror
};

/*!
@brief The Floyd-steinberg error diffusion options.
*/
struct Fs_error_diffusion_options
{
  bool is_leaky_;                         // if use leaky mode or not
  Scan_order scan_order_;
  Initial_error_mode initial_error_mode_;
};

/*!
@brief
Perform the Floyd-steinberg error diffusion method by a given threshold.
@param density  The sample point density function.
@param threshold  The threshold used for the error diffusion method.
@param options  The options used for the error diffusion method.
@param initial_errors  An array containing the errors to be diffused to the startup row.
@param last_row_errors  A vector to store the diffused errors on the last row.
@param result  A 2d array to store the bitmap result.

@details
This fuction performs the Floyd-steinberg error diffusion method.
It takes the image density function, a threshold, the specified error diffusion 
options and initial errors as input, then performs the Floyd-steinberg 
error diffusion method. It will output the result bitmap as well as the errors
that diffused to the last row.

Note: if initial_errors is an empty vector, the program will take inital errors
as zero by default. Otherwise, if the initial errors are specified, the number 
of elements in initial_errors must be the same as the width of the input 2-d
density function values.
*/
void perform_fs_error_diffusion_by_threshold(
  const SPL::Array2<double>& density, double threshold,
  const Fs_error_diffusion_options& options,
  const std::vector<double>& initial_errors,
  std::vector<double>& last_row_errors, SPL::Array2<int>& result);

/*!
@brief
Perform the Floyd-steinberg error diffusion method by a given target size.
@param density  The sample point density function.
@param num_points  The desired number of sampling points in the result bitmap.
@param tolerance  The tolerance of the result size.
@param options  The options used for the error diffusion method.
@param result  A 2d array to store the bitmap result.
@return  
Upon success, zero is returned; otherwise a non-zero integer is returned.

@details
This function will internally will call function 
get_fs_error_diffusion_startup_errors first to get the startup errors
for Floy-steinberg error diffusion.

This function performs the Floy-steinberg error diffusion method and generates
the desired number of sampling points in the result bitmap. It will call
perform_fs_error_diffusion_by_threshold iteratively and do binary search
in order to generate the desired number of sampling points in the result bitmap.

The result size will be in the range:
  (num_points - tolerance, num_points + tolerance).

The function can fail in the case that the desired size can not be reached
while doing binary search, a non-zero integer will be returned in this case.

*/
int perform_fs_error_diffusion_by_size(
  const SPL::Array2<double>& density, int num_points, int tolerance,
  const Fs_error_diffusion_options& options, SPL::Array2<int>& result);


/*!
@brief
Get the startup error for the Floyd-steinberg error diffusion
based on specified options .
@param density  The sample point density function.
@param num_points  The desired number of sampling points in the result bitmap.
@param options  The options used for the error diffusion method.
@param startup_errors  A vector to store the result startup errors.
@return  
Upon success, zero is returned; otherwise a non-zero integer is returned.
*/
int get_fs_error_diffusion_startup_errors(const SPL::Array2<double>& density,
  int num_points, const Fs_error_diffusion_options& options, 
  std::vector<double>& startup_errors);


#endif
