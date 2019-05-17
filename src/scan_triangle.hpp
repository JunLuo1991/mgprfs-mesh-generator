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

#ifndef scan_triangle_hpp
#define scan_triangle_hpp

#include "image.hpp"
#include "pixel.hpp"
#include <SPL/rasterize.hpp>
#include <vector>

/*!
@brief  A struct representing the scan statistics.
*/
struct Scan_stats
{
  Scan_stats() : error_(0) {}
  double error_;
};

/*
@brief  The avaliable error metrics
*/
enum class Error_metric
{
  vec_se,
  mean_comp_se,
  mean_comp_ae,
  luma_se,
  luma_ae,
  ycbcr_vec_se,
  ycbcr_vec_ae
};

/*!
@brief  The avaliable options while scanning triangles.
@details
This is used to help set the "flags" (bitmask) argument in the
scan_triangle_color and scan_triangle_grayscale functions.
*/
enum Scan_triangle_options
{
  write_scan_buf = 0x01,
  write_statistics = 0x02
};

/*!
@brief  Scan a triangle with rgb values on each point.
@param a_x,a_y  The x and y coordinates of the first vertex.
@param b_x,b_y  The x and y coordinates of the second vertex.
@param c_x,c_y  The x and y coordinates of the third vertex.
@param a_fs  The function values (r,g,b) of the first vertex.
@param b_fs  The function values (r,g,b) of the second vertex.
@param c_fs  The function values (r,g,b) of the third vertex.
@param ori_image  The original color image.
@param error_metric  The error metric to be used.
@param recon_image  The re-constructed image.
@param stats  A Scan_stats object to store the scan-converted statistics.
@param flags  A bitmask that selects which operation to be used.

@details
This function can be used to scan a triangle whose vertices have
rgb values. Whether to write the scan-converted data into the
recon_image, or compute the face error based on the error metrix
can be controlled by flags.

If the ori_image is not empty (i.e, ori_image.num_components() == 3),
the function will set the ori_image as the reference image, and
compute the metric error by comparing with the scan-converted value.
If the ori_image is empty (i.e., ori_image.num_components() == 0),
the funciton will set the the recon_image as the reference image.

The flags is a bitmask used to select what operation to be used.
This should be set based on the enum Scan_triangle_options. User can use
the or operator (i.g '|') to combine different options.
It has the following options:
  write_scan_buf: write scan-converted value to recon_image.
  write_statistics: compute the error based on error_metric and
                            write statistics to stats.
*/
void scan_triangle_color(int a_x, int a_y, int b_x, int b_y, int c_x, int c_y,
                         const pixel<int, 3>& a_fs, const pixel<int, 3>& b_fs,
                         const pixel<int, 3>& c_fs, const Image& ori_image,
                         Error_metric error_metric, Image& recon_image,
                         Scan_stats& stats, unsigned int flags);

/*!
@brief  Scan a triangle with grayscale values on each point.
@param a_x,a_y  The x and y coordinates of the first vertex.
@param b_x,b_y  The x and y coordinates of the second vertex.
@param c_x,c_y  The x and y coordinates of the third vertex.
@param a_fs  The function value (gray value) of the first vertex.
@param b_fs  The function value (gray value) of the second vertex.
@param c_fs  The function value (gray value) of the third vertex.
@param ori_image  The original grayscale image.
@param error_metric  The error metric to be used.
@param recon_image  The re-constructed image.
@param stats  A Scan_stats object to store the scan-converted statistics.
@param flags  A bitmask that selects which operation to be used.

@details
This function can be used to scan a triangle whose vertices have
grayscale values. Whether to write the scan-converted data into
recon_image, or compute the face error based on the error metric
can be controlled by flags.

If the ori_image is provided (e.g ori_image.num_components() == 1),
the function will set the ori_image as the reference image, which is used
to compute the metric error by comparing with the scan-converted value.
If the ori_image is not provided (e.g ori_image.num_components() == 0),
the funciton will set the the recon_image as the reference image.

The flags is a bitmask used to select what operation to be used.
This should be set based on the enum Scan_triangle_options. User can use
the or operator (i.g '|') to combine different options.
It has the following options:
  write_scan_buf: write scan-converted value to recon_image.
  write_statistics: compute the error based on error_metric and
                    write statistics to stats.
*/
void scan_triangle_grayscale(int a_x, int a_y, int b_x, int b_y, int c_x,
                             int c_y, int a_f, int b_f, int c_f,
                             const Image& ori_image, Error_metric error_metric,
                             Image& recon_image, Scan_stats& stats,
                             unsigned int flags);

#endif
