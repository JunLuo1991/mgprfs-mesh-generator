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

// This program is used for comparing the result of double convolution
// versus single convolution

#include <iostream>
#include <string>
#include <SPL/Array2.hpp>
#include <SPL/Sequence2.hpp>
#include <SPL/Sequence1.hpp>

#include "util.hpp"

typedef SPL::Sequence2<double> RealSeq2;
typedef SPL::Sequence1<double> RealSeq1;

void usage()
{
   std::cout << "This program shows the result of running double convolution\n"
	     << "versus single convolution for derivative computation.\n"
	     << "By default, using binomial filter and smooth on both directions.\n"
	     << "Usage:"
	     << "-h: print the help info\n"
	     << "-x: x order\n"
	     << "-y: y order\n"
	     << "-o: binomial smooth order\n"
	     << "-m: smooth convolve mode (default: 1)\n"
	     << "  valid options: 1(zero_ext), 3(const_ext), 4(sym_ext)\n"
	     << "-M: derivative convolve mode (default: 1)\n"
	     << "  valid options: 1(zero_ext), 3(const_ext), 4(sym_ext)\n"
	     << "NOTE: for single convolution, convolve mode is decide by derivative convolve mode";
    exit(1); 
}

int main(int argc, char** argv)
{
	int smooth_conv_mode = SPL::ConvolveMode::sameDomainZeroExt;
	int deriv_conv_mode = SPL::ConvolveMode::sameDomainZeroExt;

	// default fxx
	int x_order = 2;
	int y_order = 0;
	int smooth_order = 3;

	int c;
	while ((c = getopt(argc, argv, "m:")) >= 0) {
		switch (c) {
		case 'h':
			usage();
		case 'x':
			x_order = atoi(optarg);
			break;
		case 'y':
			y_order = atoi(optarg);
			break;
		case 'o':
			smooth_order = atoi(optarg);
			break;
		case 'm':
			smooth_conv_mode = atoi(optarg);
			break;
		case 'M':
			deriv_conv_mode = atoi(optarg);
			break;
		default:
			usage();
			break;
		}
	}


	SPL::Array2<double> input_arr(6, 6, 1.0);
	Smooth_options smooth_options;
	smooth_options.bound_policy_ = smooth_conv_mode;
	smooth_options.smooth_order_ = smooth_order;
	smooth_options.smooth_operator_ = Smooth_operator::binomial;
	smooth_options.smooth_direction_ = Smooth_direction::both;
	const int derivative_bound_policy = deriv_conv_mode;

	RealSeq1 horz_smooth_filter;
	RealSeq1 vert_smooth_filter;
	RealSeq1 horz_deriv_filter;
	RealSeq1 vert_deriv_filter;
	get_smoothing_and_derivative_filters(x_order, y_order, smooth_options,
	  horz_smooth_filter, vert_smooth_filter, horz_deriv_filter, vert_deriv_filter);

	RealSeq2 input(input_arr);

	// double convolution
     	auto smoothed_seq = SPL::convolveSeparable(input, horz_smooth_filter, vert_smooth_filter,
                	                           smooth_options.bound_policy_);
     	auto double_conv_seq = SPL::convolveSeparable(smoothed_seq, horz_deriv_filter,
             	                            	      vert_deriv_filter, derivative_bound_policy);
	auto smoothed_arr = smoothed_seq.getArray();
	auto double_conv_fxx = double_conv_seq.getArray();

	// single convolution
	auto horz_filter = SPL::convolve(horz_smooth_filter, horz_deriv_filter);
    	auto vert_filter = SPL::convolve(vert_smooth_filter, vert_deriv_filter);
    	auto single_conv_seq = SPL::convolveSeparable(input, horz_filter, vert_filter,
        	                                      derivative_bound_policy);
	auto single_conv_fxx = single_conv_seq.getArray();


	std::cout << "input array: " << input_arr << "\n\n";
	std::cout << "x_order = " << x_order << ", y_order = " << y_order << "\n"
		  << "binomial smooth order = " << smooth_order << "\n"
		  << "smooth convolve mode: " << smooth_conv_mode << "\n"
		  << "deriv convolve mode: " << deriv_conv_mode << "\n\n";

	std::cout << "double convolution:\n";
	std::cout << "horz smooth filter: " << horz_smooth_filter << "\n"
		  << "vert smooth filter: " << vert_smooth_filter << "\n";
	std::cout << "smoothed array:\n" << smoothed_arr << "\n";
	std::cout << "horz deriv filter: " << horz_smooth_filter << "\n"
		  << "vert deriv filter: " << vert_smooth_filter << "\n";
	std::cout << "final derivative:\n" << double_conv_fxx << "\n\n";

	std::cout << "single convolution:\n";
	std::cout << "horz filter: " << horz_filter << "\n"
		  << "vert filter: " << vert_filter << "\n";
	std::cout << "final derivative:\n" << single_conv_fxx << "\n\n";

	return 0;

}

