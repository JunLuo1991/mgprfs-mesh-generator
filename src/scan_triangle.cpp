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

#include "scan_triangle.hpp"
#include "pixel.hpp"
#include "rasterize_helper.hpp"
#include "util.hpp"

#include <SPL/iterator.hpp>
#include <SPL/math.hpp>
#include <SPL/rasterize.hpp>
#include <boost/rational.hpp>
#include <cmath>

void
scan_triangle_color(int a_x, int a_y, int b_x, int b_y, int c_x, int c_y,
                    const pixel<int, 3>& a_fs, const pixel<int, 3>& b_fs,
                    const pixel<int, 3>& c_fs, const Image& ori_image,
                    Error_metric error_metric, Image& recon_image,
                    Scan_stats& stats, unsigned int flags)
{

#if defined(SCAN_CONVERSION_REAL_TYPE_USE_RATIONAL)
  using Real = SPL::rasterize::Rational<long long>;
#else
  using Real = double;
#endif

  using Int_pixel_value = pixel<int, 3>;
  using Real_pixel_value = pixel<Real, 3>;

  // If write data to scan buffer (i.g. recon_image)
  bool write_data_to_buffer = flags & (Scan_triangle_options::write_scan_buf);
  // If compute the error and write to stats
  bool write_stats = flags & (Scan_triangle_options::write_statistics);

  // If original image is provided, set ref_image as a reference to ori_image,
  // otherwise set ref_image as a reference to recon_image.
  const Image& ref_image =
    (ori_image.num_components() > 0) ? ori_image : recon_image;
  assert(ref_image.num_components() == 3);

  const int xmin = 0;
  const int ymin = 0;
  const int xmax = ref_image.width() - 1;
  const int ymax = ref_image.height() - 1;
  const int prec = ref_image.precision();
  const int zmin = 0;
  const int zmax = precision_to_value(prec);

  Real_pixel_value a_zs{ static_cast<Real>(a_fs[0]), static_cast<Real>(a_fs[1]),
                         static_cast<Real>(a_fs[2]) };
  Real_pixel_value b_zs{ static_cast<Real>(b_fs[0]), static_cast<Real>(b_fs[1]),
                         static_cast<Real>(b_fs[2]) };
  Real_pixel_value c_zs{ static_cast<Real>(c_fs[0]), static_cast<Real>(c_fs[1]),
                         static_cast<Real>(c_fs[2]) };

  auto round_operator = [](const Real& x) { return x; };

  rasterize::Pixel_triangle_linear_interpolator<
    Real, Real, decltype(round_operator), 3>
    func(a_x, a_y, a_zs, b_x, b_y, b_zs, c_x, c_y, c_zs, round_operator);

  auto get_func_iter_border = [&func](int x, int y, unsigned mask) {
    return func.row_begin_border(x, y, mask);
  };
  auto get_func_iter = [&func](int x, int y) {
    return func.row_begin(x, y);
  };

  auto real_to_int = [zmin, zmax](Real x) {
    return SPL::clip<int>(round(x), zmin, zmax);
  };

  double sum_err = 0;
  // update_stat type 1: this will compute the error based on error_metric.
  auto update_stat_enabled = [&sum_err, error_metric, &real_to_int](
    const Real_pixel_value& values, const Int_pixel_value& ref_values) {

    Int_pixel_value cliped_values = { real_to_int(values[0]),
                                      real_to_int(values[1]),
                                      real_to_int(values[2]) };
    double err = 0;
    switch (error_metric) {
      case Error_metric::vec_se: {
        std::array<double, 3> se {{
            static_cast<double>(SPL::sqr(cliped_values[0] - ref_values[0])),
            static_cast<double>(SPL::sqr(cliped_values[1] - ref_values[1])),
            static_cast<double>(SPL::sqr(cliped_values[2] - ref_values[2]))
        }};
        err = se[0] + se[1] + se[2];
      } break;
      case Error_metric::mean_comp_se: {
        std::array<double, 3> se {{
            static_cast<double>(SPL::sqr(cliped_values[0] - ref_values[0])),
            static_cast<double>(SPL::sqr(cliped_values[1] - ref_values[1])),
            static_cast<double>(SPL::sqr(cliped_values[2] - ref_values[2]))
        }};
        err = (se[0] + se[1] + se[2]);
      } break;
      case Error_metric::mean_comp_ae: {
        std::array<double, 3> ae {{
            static_cast<double>(SPL::absVal(cliped_values[0] - ref_values[0])),
            static_cast<double>(SPL::absVal(cliped_values[1] - ref_values[1])),
            static_cast<double>(SPL::absVal(cliped_values[2] - ref_values[2]))
        }};
        err = (ae[0] + ae[1] + ae[2]) / 3.0;
      } break;
      case Error_metric::luma_se: {
        Real_pixel_value ref = { static_cast<Real>(ref_values[0]),
                                 static_cast<Real>(ref_values[1]),
                                 static_cast<Real>(ref_values[2]) };
        err = SPL::sqr(round(rgb_to_luma(values)) - round(rgb_to_luma(ref)));
      } break;
      case Error_metric::luma_ae: {
        Real_pixel_value ref = { static_cast<Real>(ref_values[0]),
                                 static_cast<Real>(ref_values[1]),
                                 static_cast<Real>(ref_values[2]) };
        err = SPL::absVal(round(rgb_to_luma(values)) - round(rgb_to_luma(ref)));
      } break;
      case Error_metric::ycbcr_vec_se: {
        Real_pixel_value ref = { static_cast<Real>(ref_values[0]),
                                 static_cast<Real>(ref_values[1]),
                                 static_cast<Real>(ref_values[2]) };
        Real_pixel_value ycbcr_value;
        Real_pixel_value ycbcr_ref;
        rgb_to_ycbcr(values, ycbcr_value);
        rgb_to_ycbcr(ref, ycbcr_ref);
        err = SPL::sqr(round(ycbcr_value[0]) - round(ycbcr_ref[0])) +
              SPL::sqr(round(ycbcr_value[1]) - round(ycbcr_ref[1])) +
              SPL::sqr(round(ycbcr_value[2]) - round(ycbcr_ref[2]));
      } break;
      case Error_metric::ycbcr_vec_ae: {
        Real_pixel_value ref = { static_cast<Real>(ref_values[0]),
                                 static_cast<Real>(ref_values[1]),
                                 static_cast<Real>(ref_values[2]) };
        Real_pixel_value ycbcr_value;
        Real_pixel_value ycbcr_ref;
        rgb_to_ycbcr(values, ycbcr_value);
        rgb_to_ycbcr(ref, ycbcr_ref);
        double se = SPL::sqr(round(ycbcr_value[0]) - round(ycbcr_ref[0])) +
                    SPL::sqr(round(ycbcr_value[1]) - round(ycbcr_ref[1])) +
                    SPL::sqr(round(ycbcr_value[2]) - round(ycbcr_ref[2]));
        err = std::sqrt(se);
      } break;
      default:
        std::cerr << "invalid error metric\n";
        std::exit(1);
    }
    sum_err += err;
  };

  // update_stat type 2 : this will not compute the error
  auto update_stat_disabled = [](const Real_pixel_value& value,
                                 const Int_pixel_value& ref_value) {};

  auto get_buf_iter = [&ref_image](int x, int y) {
    return rasterize::Stat_input_iterator<3>(ref_image, x, y);
  };

  // get_stat_iter type 1: this will compute the error and update the stats
  rasterize::Stat_collector<decltype(update_stat_enabled),
                            decltype(get_buf_iter), Real, 3>
    get_stat_iter_enable_update(update_stat_enabled, get_buf_iter);
  // get_stat_iter type 2: this will not compute the error
  rasterize::Stat_collector<decltype(update_stat_disabled),
                            decltype(get_buf_iter), Real, 3>
    get_stat_iter_disable_update(update_stat_disabled, get_buf_iter);

  // get_scan_buf_iter type 1:
  // this will return an output iterator that point to a specific image pixel.
  auto get_scan_buf_iter_enable_write = [&recon_image, &real_to_int](int x, int y) {
    return rasterize::Scan_buf_output_iterator<Real, decltype(real_to_int), 3>(
      recon_image, x, y, real_to_int);
  };
  // get_scan_buf_iter type 2: this will return an null_output_iterator.
  auto get_scan_buf_iter_disable_write = [](int x, int y) {
    return SPL::null_output_iterator();
  };

  if (write_data_to_buffer) {
    if (write_stats) {
      SPL::rasterize::Triangle_scan_line<
        int, decltype(get_func_iter_border), decltype(get_func_iter),
        decltype(get_scan_buf_iter_enable_write),
        decltype(get_stat_iter_enable_update)>
        scan_line_a(get_func_iter_border, get_func_iter,
                    get_scan_buf_iter_enable_write,
                    get_stat_iter_enable_update);

      SPL::rasterize::scan_triangle_in_iso_rectangle_domain<int>(
        a_x, a_y, b_x, b_y, c_x, c_y, xmin, ymin, xmax, ymax,
        std::ref(scan_line_a));
    } else {
      SPL::rasterize::Triangle_scan_line<
        int, decltype(get_func_iter_border), decltype(get_func_iter),
        decltype(get_scan_buf_iter_enable_write),
        decltype(get_stat_iter_disable_update)>
        scan_line_a(get_func_iter_border, get_func_iter,
                    get_scan_buf_iter_enable_write,
                    get_stat_iter_disable_update);

      SPL::rasterize::scan_triangle_in_iso_rectangle_domain<int>(
        a_x, a_y, b_x, b_y, c_x, c_y, xmin, ymin, xmax, ymax,
        std::ref(scan_line_a));
    }
  } else {
    if (write_stats) {
      SPL::rasterize::Triangle_scan_line<
        int, decltype(get_func_iter_border), decltype(get_func_iter),
        decltype(get_scan_buf_iter_disable_write),
        decltype(get_stat_iter_enable_update)>
        scan_line_a(get_func_iter_border, get_func_iter,
                    get_scan_buf_iter_disable_write,
                    get_stat_iter_enable_update);

      SPL::rasterize::scan_triangle_in_iso_rectangle_domain<int>(
        a_x, a_y, b_x, b_y, c_x, c_y, xmin, ymin, xmax, ymax,
        std::ref(scan_line_a));
    } else {
      SPL::rasterize::Triangle_scan_line<
        int, decltype(get_func_iter_border), decltype(get_func_iter),
        decltype(get_scan_buf_iter_disable_write),
        decltype(get_stat_iter_disable_update)>
        scan_line_a(get_func_iter_border, get_func_iter,
                    get_scan_buf_iter_disable_write,
                    get_stat_iter_disable_update);

      SPL::rasterize::scan_triangle_in_iso_rectangle_domain<int>(
        a_x, a_y, b_x, b_y, c_x, c_y, xmin, ymin, xmax, ymax,
        std::ref(scan_line_a));
    }
  }

  if (write_stats) {
    stats.error_ = sum_err;
  }
}

void
scan_triangle_grayscale(int a_x, int a_y, int b_x, int b_y, int c_x, int c_y,
                        int a_f, int b_f, int c_f, const Image& ori_image,
                        Error_metric error_metric, Image& recon_image,
                        Scan_stats& stats, unsigned int flags)
{
#if defined(SCAN_CONVERSION_REAL_TYPE_USE_RATIONAL)
  using Real = SPL::rasterize::Rational<long long>;
#else
  using Real = double;
#endif

  using Int_pixel_value = pixel<int, 1>;
  using Real_pixel_value = pixel<Real, 1>;

  // If write data to scan buffer (i.g. recon_image)
  bool write_data_to_buffer = flags & (Scan_triangle_options::write_scan_buf);
  // If compute the error and write statistics to stats
  bool write_stats = flags & (Scan_triangle_options::write_statistics);

  // If original image is provided, set ref_image as original image,
  // otherwise set as recon_image
  const Image& ref_image =
    (ori_image.num_components() > 0) ? ori_image : recon_image;
  assert(ref_image.num_components() == 1);

  const int xmin = 0;
  const int ymin = 0;
  const int xmax = ref_image.width() - 1;
  const int ymax = ref_image.height() - 1;
  const int prec = ref_image.precision();
  const int zmin = 0;
  const int zmax = precision_to_value(prec);

  pixel<Real, 1> a_z{ static_cast<Real>(a_f) };
  pixel<Real, 1> b_z{ static_cast<Real>(b_f) };
  pixel<Real, 1> c_z{ static_cast<Real>(c_f) };

  auto round_operator = [](const Real& x) {
    return x;
  };

  rasterize::Pixel_triangle_linear_interpolator<Real, Real,
                                                decltype(round_operator), 1>
    func(a_x, a_y, a_z, b_x, b_y, b_z, c_x, c_y, c_z, round_operator);

  auto get_func_iter_border = [&func](int x, int y, unsigned mask) {
    return func.row_begin_border(x, y, mask);
  };
  auto get_func_iter = [&func](int x, int y) {
    return func.row_begin(x, y);
  };

  auto real_to_int = [zmin, zmax](Real x) {
    return SPL::clip<int>(round(x), zmin, zmax);
  };
  double se = 0;
  double ae = 0;
  auto update_stat_enabled = [&se, &ae, &real_to_int](
    const Real_pixel_value& value, const Int_pixel_value& ref_value) {
    se += SPL::sqr(real_to_int(value[0]) - ref_value[0]);
    ae += SPL::absVal(real_to_int(value[0]) - ref_value[0]);
  };
  auto update_stat_disabled = [](const Real_pixel_value& value,
                                 const Int_pixel_value& ref_value) {};

  auto get_buf_iter = [&ref_image](int x, int y) {
    return rasterize::Stat_input_iterator<1>(ref_image, x, y);
  };
  // get_stat_iter type 1: this will compute the error and update the stats
  rasterize::Stat_collector<decltype(update_stat_enabled),
                            decltype(get_buf_iter), Real, 1>
    get_stat_iter_enable_update(update_stat_enabled, get_buf_iter);
  // get_stat_iter type 2: this will not compute the error
  rasterize::Stat_collector<decltype(update_stat_disabled),
                            decltype(get_buf_iter), Real, 1>
    get_stat_iter_disable_update(update_stat_disabled, get_buf_iter);

  // get_scan_buf_iter type 1:
  // this will return an output iterator that point to a specific image pixel.
  auto get_scan_buf_iter_enable_write = [&recon_image, &real_to_int](int x, int y) {
    return rasterize::Scan_buf_output_iterator<Real, decltype(real_to_int), 1>(
      recon_image, x, y, real_to_int);
  };
  // get_scan_buf_iter type 2: this will return an null_output_iterator.
  auto get_scan_buf_iter_disable_write = [](int x, int y) {
    return SPL::null_output_iterator();
  };

  if (write_data_to_buffer) {
    if (write_stats) {
      SPL::rasterize::Triangle_scan_line<
        int, decltype(get_func_iter_border), decltype(get_func_iter),
        decltype(get_scan_buf_iter_enable_write),
        decltype(get_stat_iter_enable_update)>
        scan_line_a(get_func_iter_border, get_func_iter,
                    get_scan_buf_iter_enable_write,
                    get_stat_iter_enable_update);

      SPL::rasterize::scan_triangle_in_iso_rectangle_domain<int>(
        a_x, a_y, b_x, b_y, c_x, c_y, xmin, ymin, xmax, ymax,
        std::ref(scan_line_a));
    } else {
      SPL::rasterize::Triangle_scan_line<
        int, decltype(get_func_iter_border), decltype(get_func_iter),
        decltype(get_scan_buf_iter_enable_write),
        decltype(get_stat_iter_disable_update)>
        scan_line_a(get_func_iter_border, get_func_iter,
                    get_scan_buf_iter_enable_write,
                    get_stat_iter_disable_update);

      SPL::rasterize::scan_triangle_in_iso_rectangle_domain<int>(
        a_x, a_y, b_x, b_y, c_x, c_y, xmin, ymin, xmax, ymax,
        std::ref(scan_line_a));
    }
  } else {
    if (write_stats) {
      SPL::rasterize::Triangle_scan_line<
        int, decltype(get_func_iter_border), decltype(get_func_iter),
        decltype(get_scan_buf_iter_disable_write),
        decltype(get_stat_iter_enable_update)>
        scan_line_a(get_func_iter_border, get_func_iter,
                    get_scan_buf_iter_disable_write,
                    get_stat_iter_enable_update);

      SPL::rasterize::scan_triangle_in_iso_rectangle_domain<int>(
        a_x, a_y, b_x, b_y, c_x, c_y, xmin, ymin, xmax, ymax,
        std::ref(scan_line_a));
    } else {
      SPL::rasterize::Triangle_scan_line<
        int, decltype(get_func_iter_border), decltype(get_func_iter),
        decltype(get_scan_buf_iter_disable_write),
        decltype(get_stat_iter_disable_update)>
        scan_line_a(get_func_iter_border, get_func_iter,
                    get_scan_buf_iter_disable_write,
                    get_stat_iter_disable_update);

      SPL::rasterize::scan_triangle_in_iso_rectangle_domain<int>(
        a_x, a_y, b_x, b_y, c_x, c_y, xmin, ymin, xmax, ymax,
        std::ref(scan_line_a));
    }
  }

  if (write_stats) {
    double err;
    if (error_metric == Error_metric::mean_comp_se ||
        error_metric == Error_metric::luma_se ||
        error_metric == Error_metric::ycbcr_vec_se ||
        error_metric == Error_metric::vec_se) {
      err = se;
    } else if (error_metric == Error_metric::mean_comp_ae ||
               error_metric == Error_metric::luma_ae ||
               error_metric == Error_metric::ycbcr_vec_ae) {
      err = ae;
    }

    stats.error_ = err;
  }
}
