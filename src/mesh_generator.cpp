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


#include "mesh_generator.hpp"
#include "pixel.hpp"
#include "util.hpp"
#include <SPL/math.hpp>
#include <boost/program_options.hpp>
#include <map>
#include <stdexcept>

Mesh_generator::Mesh_generator()
{
  clear();
}

int Mesh_generator::operator()(std::istream& in, std::ostream& out, int argc,
                               const char* const* argv)
{
  return main(in, out, argc, argv);
}

int Mesh_generator::main(std::istream& in, std::ostream& out, int argc,
                     const char* const* argv)
{
  clear();

  // process options
  if (auto status = process_options(in, out, argc, argv); status) {
    return status;
  }
  log_message(1, "process options complete.\n\n");

  log_message(1, "start making mesh...\n");
  // generate the mesh.
  if (auto status = make_mesh(); status) {
    return status;
  }
  log_message(1, "make mesh complete.\n\n");

  log_message(1, "output mesh...\n");
  // write mesh to output stream in header + OFF format.
  if (auto status = output_mesh(out); status) {
    return status;
  }
  log_message(1, "output mesh complete.\n");
  return 0;
}

void Mesh_generator::clear()
{
  if (history_stream_) {
    history_stream_.reset();
  }
}

int Mesh_generator::process_options(std::istream& in, std::ostream& out,
                                    int argc, const char* const* argv)
{
  std::string history_file_name;
  int history_level = 0;
  int debug_level = 0;
  int target_size = -1;
  double target_error = -1.0;
  double sampling_density = -1.0;
  int initial_size = -1;
  double initial_density = 100;
  double relative_initial_density = -1.0;
  std::string mesh_generation_policy;
  std::string initial_generator;
  int ed_size_tolerance;
  int smooth_order;
  std::string smooth_operator;
  std::string smooth_direction;
  std::string smooth_bound_policy;
  std::string derivative_bound_policy;
  int use_double_convolution;
  std::string ed_strategy;
  std::string ed_scan_order;
  std::string ed_leaky_mode;
  std::string ed_initial_error_mode;
  std::string error_metric;
  int enable_bad_point_removal;

  namespace po = boost::program_options;
  po::options_description desc{ "Options" };
  desc.add_options()
    ("help,h", "print help information only.\n")
    ("debug-level", po::value<int>(&debug_level)->default_value(0),
     "set debug level.\n")
    ("history-file", po::value<std::string>(&history_file_name), 
     "set history file.\n")
    ("history-level", po::value<int>(&history_level)->default_value(0),
     "set history level.\n")
    ("density", po::value<double>(&sampling_density),
     "specify target mesh sampling density (in percent).\n")
    ("size", po::value<int>(&target_size),
     "specify target mesh size (in number of samples).\n")
    ("error", po::value<double>(&target_error),
     "specify target mesh error threshold.\n")
    ("relative-initial-density", po::value<double>(&relative_initial_density),
     "set the initial density that relative to the target density. this will "
     "set the initial mesh density as relative_initial_density/ (density * "
     "100).\n")
    ("initial-density", po::value<double>(&initial_density)->default_value(100.0),
     "set initial mesh density (in percent).\n")
    ("initial-size", po::value<int>(&initial_size),
     "specify initial mesh size (in number of samples).\n")
    ("mesh-generation-policy",
     po::value<std::string>(&mesh_generation_policy)->default_value("mcmg"),
     "specify the method used for generaring the mesh. Note: This option can "
     "take influence only when the input image is a color image\n."
     "valid options include:\n"
     "1. mcmg ... multi component mesh generation. This will process multiple"
     " components together by using the color image mesh generator\n"
     "2. lsmg ... luma space mesh generation. This will convert the color image"
     " into luma space first, then call the grayscale image mesh generator to"
     " get a mesh, finally replace all sample point values with color values in"
     " the mesh.\n\n")
    ("initial-generator",
     po::value<std::string>(&initial_generator)->default_value("ed"),
     "specify the method used for generating initial mesh point set.\n"
     "valid options include:\n"
     "1. all ........ Set all sample points on the sampling grid as the"
     " initial point set.\n"
     "2. random ..... Randomly select points on the sampling grid as the"
     " initial point set.\n"
     "3. uniform .... Uniformly select points on the sampling grid as"
     " the initial point set.\n"
     "4. ed ......... Use error diffusion method to select the initial"
     " point set.\n")
    ("ed-size-tolerance", po::value<int>(&ed_size_tolerance)->default_value(10),
     "specify the size tolerance of the mesh that generated by the error"
     " diffusion method. NOTE: This option can take influence only when the"
     " initial-generator is ed.\n")
    ("smooth-order", po::value<int>(&smooth_order)->default_value(3),
     "specify the smooth order of the smoothing filter. NOTE: This option can "
     "take influence only when the initial-generator is ed. The smooth order must "
     "be a positive odd integer value. If one, no smoothing will be applied.\n")
    ("smooth-operator",
     po::value<std::string>(&smooth_operator)->default_value("binomial"),
     "specify the smoothing operator. Note: this option can take influence only "
     "when the initial-generator is ed.\n"
     "Valid options include:\n"
     "1. binomial ... use binomial filter.\n"
     "2. mean  ...... use mean filter.\n")
    ("smooth-direction",
     po::value<std::string>(&smooth_direction)->default_value("both"),
     "specify the smoothing direction. Note: this option can take influence "
     "only when the initial-generator is ed.\n"
     "Valid options include:\n"
     "1. orthogonal ... smooth on orthogonal direction.\n"
     "2. both ... smooth on both horizontal and vertical directions.\n")
    ("smooth-bound-policy",
     po::value<std::string>(&smooth_bound_policy)->default_value("zero_ext"),
     "specify the boundary handling policy used for smoothing. NOTE: This option "
     " can take influence only when the initial-generator is ed.\n"
     "Valid options include:\n"
     "1. zero_ext ...... zero extension.\n"
     "2. const_ext ..... constant extension.\n"
     "3. sym_ext ....... symmetric extension.\n")
    ("derivative-bound-policy",
     po::value<std::string>(&derivative_bound_policy)->default_value("zero_ext"),
     "specify the boundary handling policy used for partial derivative "
     "computation. NOTE: This option can take influence only when the "
     "initial-generator is ed.\n"
     "Valid options include:\n"
     "1. zero_ext ...... zero extension.\n"
     "2. const_ext ..... constant extension.\n"
     "3. sym_ext ....... symmetric extension.\n")
    ("use-double-convolution",
     po::value<int>(&use_double_convolution)->default_value(0),
     "specify if do smoothing and derivative operations together in one "
     "convolution or seperately in two convolution. valid options include 0 and 1.\n\n"
     "If 0, the program will do smoothing and derivative operations in a single "
     "convolution. The boundary policy will be based on the option "
     "derivative-bound-policy at this time.\n\n"
     "If 1, the program will compute the partial derivative by using two "
     "convolution operations. Basically it will smooth the image first "
     "using one convolution and then do another convolution with derivative filters.\n\n"
     "NOTE: This option can take influence only when the initial-generator is "
     "ed.\n")
    ("ed-strategy",
     po::value<std::string>(&ed_strategy)->default_value("max_comp_mmsodd_ed"),
     "specify the error diffusion strategy. NOTE: this option can take influence "
     "only when the initial-generator is ed. When the input image is a color "
     "image, each option corresponds to one color error diffusion strategy; When "
     "the input image is a grayscale image, if gray_comp_laplacian_ed, or "
     "mean_comp_laplacian_ed, or max_comp_laplacian_ed is selected it will "
     "take laplacian of the image as the density function, otherwise it will "
     "take mmsodd of the image as the density function.\n"
     "Valid options include:\n\n"
     "1. gray_comp_mmsodd_ed ... convert the color image to a grayscale image "
     "first, then use mmsodd of the grayscale image as the density function "
     "and run Floyd-steinberg error diffusion method.\n\n"
     "2. gray_comp_laplacian_ed ... convert the color image to a grayscale "
     "image first, then use laplacian of the grayscale image as the density "
     "function and run Floyd-steinberg error diffusion method.\n\n"
     "3. vector_space_comb .... compute the maximum directional graident of the "
     "vector space of R,G,B components as the density function and run "
     "Floyd-steinberg error diffusion method.\n\n"
     "4. max_comp_mmsodd_ed ... compute the max mmsodd of all components as the "
     "density function and run Floyd-steinberg error diffusion method.\n\n"
     "5. mean_comp_mmsodd_ed .. compute the mean mmsodd of all components as the "
     "density function and run Floyd-steinberg error diffusion method.\n\n"
     "6. max_comp_laplacian_ed ... compute the max laplacian of all components "
     "as the density function and run Floyd-steinberg error diffusion method.\n\n"
     "7. mean_comp_laplacian_ed .. compute the mean laplacian of all components "
     "as the density function and run Floyd-steinberg error diffusion method.\n\n"
     "8. three_comps_union_mmsodd .... compute mmsodd of each "
     "component separately and take as three density functions; and run "
     "Floyd-steinberg error diffusion method on each of R,G,B components "
     "and get the union set. The result mesh will have desired number of points.\n\n"
     "9. three_comps_union_laplacian .... compute laplacian "
     "of each component separately and take as three density functions; and run "
     "Floyd-steinberg error diffusion method on each of R,G,B components and get "
     "the union set. The result mesh will have desired number of points.\n\n")
    ("ed-scan-order",
     po::value<std::string>(&ed_scan_order)->default_value("serpentine"),
     "specify the error diffusion scan order. NOTE: This opiton can take "
     "influence only when the initial-generator is ed.\n"
     "Valid options include:\n"
     "1. raster ... raster scan order.\n"
     "2. serpentine ... serpentine scan order.\n")
    ("ed-leaky-mode",
     po::value<std::string>(&ed_leaky_mode)->default_value("no_leaky"),
     "specify the leaky mode for error diffusion. NOTE: This can take"
     " influence only when the initial-generator is ed.\n"
     "Valid options include:\n"
     "1. leaky ...... leaky mode\n"
     "2. no_leaky ... no_leaky mode\n")
    ("ed-initial-error-mode",
     po::value<std::string>(&ed_initial_error_mode)->default_value("mirror"),
     "specify the initial error mode for error diffusion. NOTE: this can"
     " take influence only when the initial-generator is ed.\n"
     "Valid options include: \n"
     "1. zero ..... initial errors are all zero.\n"
     "2. random ... initial errors are randomly generated.\n"
     "3. mirror ... initial errors are based on doing Floyd-steinberg"
     " error diffusion first on the mirror image of the original image."
     " Then take the diffused errors on the last row as the inital errors\n")
    ("error-metric",
     po::value<std::string>(&error_metric)->default_value("mean_comp_se"),
     "specify the error metric used to compute face error.\n"
     "valid options include:\n"
     "1. mean_comp_se ... mean squared error of all components.\n"
     "2. mean_comp_ae ... mean absolute error of all components.\n"
     "3. luma_se ........ luma space squared error.\n"
     "4. luma_ae ........ luma space absolute error.\n"
     "5. ycbcr_vec_se ... YCbCr vector space squared error.\n"
     "6. ycbcr_vec_ae ... YCbCr vector space absolute error.\n"
     "7. vec_se ... vector absolute error.\n"
     "Note: for grayscale image, there is no difference between every "
     "option of *_se, they all have the same result (i.g se). Same for *_ae.\n")
    ("bad-point-removal",
     po::value<int>(&enable_bad_point_removal)->default_value(1),
     "0/1 to specify if to enable bad point removal or not\n");
 
  po::variables_map vm;

  // Parse the mesh generation options.
  try {
    po::store(po::parse_command_line(argc, argv, desc), vm);
    po::notify(vm);
  } catch (const po::error& e) {
    // An error occurred during option parsing.
    std::cerr << e.what() << '\n';
    return -1;
  }

  // If the help option was specified, print help information and
  // return without attempting mesh generation.
  if (vm.count("help")) {
    std::cerr << desc << '\n';
    return 1;
  }

  // Parse the debug level
  debug_level_ = (vm.count("debug-level")) ? debug_level : 0;

  // Parse the history level
  history_level_ = (vm.count("history-level")) ? history_level : 0;

  if (debug_level >= 1) {
    std::cerr << "debug level: " << debug_level << '\n'
              << "history level: " << history_level << '\n'
              << "history file name: " << history_file_name << '\n'
              << "number of samples: " << target_size << '\n'
              << "sampling density: " << sampling_density << '\n'
              << "target error: " << target_error << '\n'
              << "relative initial density: " << relative_initial_density << '\n'
              << "initial density: " << initial_density << '\n'
              << "initial size: " << initial_size << '\n'
              << "mesh generation policy: " << mesh_generation_policy << '\n'
              << "initial generator: " << initial_generator << '\n'
              << "ed size tolerance: " << ed_size_tolerance << '\n'
              << "smooth operator: " << smooth_operator << '\n'
              << "smooth order: " << smooth_order << '\n'
              << "smooth direction: " << smooth_direction << '\n'
              << "smooth bound policy: " << smooth_bound_policy << '\n'
              << "derivative bound policy: " << derivative_bound_policy << '\n'
              << "use double convolution: " << use_double_convolution << '\n'
              << "ed strategy: " << ed_strategy << '\n'
              << "ed scan order: " << ed_scan_order << '\n'
              << "ed leaky mode: " << ed_leaky_mode << '\n'
              << "error metric: " << error_metric << '\n'
              << "bad point removal: " << enable_bad_point_removal << '\n';
  }

  // Parse the mesh generation policy
  if (mesh_generation_policy == "mcmg") {
    mesh_generation_policy_ = Mesh_generation_policy::multi_comp_mesh_generation;
  } else if (mesh_generation_policy == "lsmg") {
    mesh_generation_policy_ = Mesh_generation_policy::luma_space_mesh_generation;
  } else {
    std::cerr << "COMMAND ERROR: invalid argument for mesh-generation-policy."
              << mesh_generation_policy
              << ". Valid options include: mcmg, lsmg\n";
    return -1;
  }

  // Read image from input stream in PNM/PPM/PGM format.
  if (auto status = ori_input_image_.input(in); status) {
    std::cerr << "Read input image failed\n";
    return status;
  }

  // Set input_image_
  if (ori_input_image_.num_components() == 3 &&
      (mesh_generation_policy_ == Mesh_generation_policy::luma_space_mesh_generation ||
      ed_strategy == "gray_comp_mmsodd_ed" || ed_strategy == "gray_comp_laplacian_ed"))
  {
    rgb_to_gray(ori_input_image_, input_image_);
  } else {
    input_image_ = ori_input_image_;
  }

  int image_size = input_image_.width() * input_image_.height();

  // Ensure that some termination criterion for mesh generation has
  // been specified.
  if (!vm.count("density") && !vm.count("size") && !vm.count("error")) {
    std::cerr << "COMMAND ERROR: no termination criterion specified for "
                 "mesh generation\n";
    return -1;
  }

  // Ensure not both density and size are specified
  if (vm.count("density") && vm.count("size")) {
    std::cerr << "COMMAND ERROR: density and size can not be specified "
                 "at the same time\n";
    return -1;
  }

  // Parse the target size
  if (vm.count("size")) {
    if (target_size <= 0 || target_size > image_size) {
      std::cerr << "COMMAND ERROR: target size for the current mesh should be "
                   " in the range 4~" << image_size << ".\n";
      return -1;
    }
    target_size_ = target_size;
  } else if (vm.count("density")) {
    if (sampling_density <= 0 || sampling_density > 100) {
      std::cerr << "COMMAND ERROR: density should be in the range (0,100].\n";
      return -1;
    }
    target_size_ = 0.01 * sampling_density * image_size;
  } else {
    target_size_ = -1;
  }

  // Parse the target error
  if (vm.count("error")) {
    if (target_error <= 0) {
      std::cerr << "COMMAND ERROR: error should be a positive value.\n";
      return -1;
    }
    target_error_ = target_error;
  } else {
    target_error_ = -1.0;
  }

  if (target_size_ > 0 && target_error_ > 0) {
    std::cerr << "Warning: target size and target error are both specified. "
              << "Program will terminate once any of them is satisfied.\n";
  }

  if (vm.count("relative-initial-density") + vm.count("initial-size") +
      vm.count("initial-density") > 1) {
    std::cerr << "COMMAND ERROR: multiple versions of initial size is specified! "
              "Note: options relative-initial-density, initial-density and "
              "initial-size can not be specified at the same time.\n";
    return -1;
  }

  // Parse the initial size
  if (vm.count("relative-initial-density")) {
      if (relative_initial_density <= 0) {
        std::cerr << "COMMAND ERROR: relative-initial-density should be "
                     "a positive value.\n";
        return -1;
      }
      initial_size_ = 0.01 * relative_initial_density * target_size_;
      initial_size_ = SPL::clip<int>(initial_size_, 0, image_size);
  } else if (vm.count("initial-size")) {
    if (initial_size < 4 || initial_size > image_size) {
      std::cerr << "COMMAND ERROR: initial-size should be an integer that is "
                   "in the range 4~" << image_size << ".\n";
      return -1;
    }
    initial_size_ = initial_size;
  } else if (vm.count("initial-density")) {
    if ((initial_density <= 0) || (initial_density > 100)) {
      std::cerr << "COMMAND ERROR: initial-density should be in the range (0,100].\n";
      return -1;
    }
    initial_size_ = 0.01 * initial_density * image_size;
  } else {
    std::cerr << "COMMAND ERROR: no initial mesh size criteria is specified. "
                 "Note: initial mesh size can be specified by options: "
                 "--relative-initial-size, --initial-size, --initial-density\n";
    return -1;
  }

  // Ensure initial_size_ is greater than target_size_
  if (target_size_ > 0 && initial_size_ < target_size_) {
    std::cerr << "COMMAND ERROR: specified initial size is less than target size. "
              "(initial size = " << initial_size_ << ", target size = " 
              << target_size_ << ").\n";
    return -1;
  }

  // Parse the initial point seleciton policy
  if (initial_generator == "all") {
    initial_point_selection_policy_ = Initial_point_selection_policy::all;
  } else if (initial_generator == "random") {
    initial_point_selection_policy_ = Initial_point_selection_policy::random;
  } else if (initial_generator == "uniform") {
    initial_point_selection_policy_ = Initial_point_selection_policy::uniform;
  } else if (initial_generator == "ed") {
    initial_point_selection_policy_ = Initial_point_selection_policy::ed;
  } else {
    std::cerr << "COMMAND ERROR: invalid option for initial-generator. "
                 "valid options include: all, random, uniform, ed.\n";
    return -1;
  }

  // Parse the ed size tolerance
  ed_size_tolerance_ = 
    vm.count("ed-size-tolerance") ? std::abs(ed_size_tolerance) : 10;

  // Pase the smooth order
  if ((smooth_order <= 0) || (smooth_order % 2 == 0)) {
    std::cerr << "COMMAND ERROR: smooth-order should be "
                 "a positive odd number.\n";
    return -1;
  }
  smooth_options_.smooth_order_ = vm.count("smooth-order") ? smooth_order : 3;

  // Parse the smoothing operator
  if (smooth_operator == "binomial") {
    smooth_options_.smooth_operator_ = Smooth_operator::binomial;
  } else if (smooth_operator == "mean") {
    smooth_options_.smooth_operator_ = Smooth_operator::mean;
  } else {
    std::cerr << "COMAMND ERROR: invalid option for smooth-operator: "
              << smooth_operator << ". Valid options include: binomial, mean.\n";
    return -1;
  }

  // Parse the smoothing direction
  if (smooth_direction == "orthogonal") {
    smooth_options_.smooth_direction_ = Smooth_direction::orthogonal;
  } else if (smooth_direction == "both") {
    smooth_options_.smooth_direction_ = Smooth_direction::both;
  } else {
    std::cerr << "COMMAND ERROR: invalid smooth-direction argument: "
              << smooth_direction << ". valid options include: orthogonal, both.\n";
    return -1;
  }

  // Parse the smoothing boundary handling policy.
  if (smooth_bound_policy == "zero_ext") {
    smooth_options_.bound_policy_ = SPL::ConvolveMode::sameDomainZeroExt;
  } else if (smooth_bound_policy == "const_ext") {
    smooth_options_.bound_policy_ = SPL::ConvolveMode::sameDomainConstExt;
  } else if (smooth_bound_policy == "sym_ext") {
    smooth_options_.bound_policy_ = SPL::ConvolveMode::sameDomainSymExt0;
  } else {
    std::cerr << "COMAMND ERROR: invalid option for smooth-bound-policy: "
              << smooth_bound_policy
              << ". valid options include: zero_ext, const_ext, sys_ext.\n";
    return -1;
  }

  // Parse the derivative boundary handling policy.
  if (derivative_bound_policy == "zero_ext") {
    derivative_bound_policy_ = SPL::ConvolveMode::sameDomainZeroExt;
  } else if (derivative_bound_policy == "const_ext") {
    derivative_bound_policy_ = SPL::ConvolveMode::sameDomainConstExt;
  } else if (derivative_bound_policy == "sym_ext") {
    derivative_bound_policy_ = SPL::ConvolveMode::sameDomainSymExt0;
  } else {
    std::cerr << "COMAMND ERROR: invalid option for derivative-bound-policy: "
              << derivative_bound_policy
              << ". valid options include: zero_ext, const_ext, sys_ext.\n";
    return -1;
  }

  // Parse use_double_convolution
  if (vm.count("use-double-convolution")) {
    use_double_convolution_ = use_double_convolution ? true: false; 
  }

  // create an ed_strategy lookup table
  std::map<std::string, Ed_strategy> ed_strategy_table;
  {
    ed_strategy_table["gray_comp_mmsodd_ed"] = Ed_strategy::gray_comp_mmsodd_ed;
    ed_strategy_table["gray_comp_laplacian_ed"] = Ed_strategy::gray_comp_laplacian_ed;
    ed_strategy_table["vector_space_comb"] = Ed_strategy::vector_space_combination;
    ed_strategy_table["max_comp_mmsodd_ed"] = Ed_strategy::max_component_mmsodd_ed;
    ed_strategy_table["max_comp_laplacian_ed"] = Ed_strategy::max_component_laplacian_ed;
    ed_strategy_table["mean_comp_mmsodd_ed"] = Ed_strategy::mean_component_mmsodd_ed;
    ed_strategy_table["mean_comp_laplacian_ed"] = Ed_strategy::mean_component_laplacian_ed;
    ed_strategy_table["three_comps_union_mmsodd"] = Ed_strategy::three_components_union_mmsodd;
    ed_strategy_table["three_comps_union_laplacian"] = Ed_strategy::three_components_union_laplacian;
  }
  // Parse the error diffusion strategy option.
  if(ed_strategy_table.find(ed_strategy) != ed_strategy_table.end()) {
    ed_strategy_ = ed_strategy_table[ed_strategy];
  } else {
    std::cerr << "COMMAND ERROR: invalid option for ed-strategy: "
              << ed_strategy << '\n';
    return -1;
  }

  // Parse the error diffusion scan order.
  if (ed_scan_order == "raster") {
    fs_ed_options_.scan_order_ = Scan_order::raster;
  } else if (ed_scan_order == "serpentine") {
    fs_ed_options_.scan_order_ = Scan_order::serpentine;
  } else {
      std::cerr << "COMMAND ERROR: invalid option for ed-scan-order: "
                << ed_scan_order
                << ". valid options include: raster, serpentine.\n";
      return -1;
  }

  // Parse if to use leaky / no_leaky mode for error diffusion
  if (ed_leaky_mode == "leaky") {
    fs_ed_options_.is_leaky_ = true;
  } else if (ed_leaky_mode == "no_leaky") {
    fs_ed_options_.is_leaky_ = false;
  } else {
    std::cerr << "COMMAND ERROR: invalid option for ed-leaky-mode: "
              << ed_leaky_mode
              << ". valid options include: leaky, no_leaky.\n";
    return -1;
  }

  // Parse ed initial error mode for error diffusion
  if (ed_initial_error_mode == "zero") {
    fs_ed_options_.initial_error_mode_ = Initial_error_mode::zero;
  } else if (ed_initial_error_mode == "random") {
    fs_ed_options_.initial_error_mode_ = Initial_error_mode::random;
  } else if (ed_initial_error_mode == "mirror") {
    fs_ed_options_.initial_error_mode_ = Initial_error_mode::mirror;
  } else {
    std::cerr << "COMAMND ERROR: invalid argument for ed-initial-error-mode:"
              << ed_initial_error_mode
              << ". valid options include: zero, random, mirror.\n";
    return -1;
  }

  // Parse the history file name
  if (history_file_name.size()) {
    history_stream_ =
      std::make_unique<std::ofstream>(history_file_name.c_str());
    if (!*history_stream_) {
      return -1;
    }
  }

  // create the error metric lookup table
  std::map<std::string, Error_metric> error_metric_table;
  {
    error_metric_table["mean_comp_se"] = Error_metric::mean_comp_se;
    error_metric_table["mean_comp_ae"] = Error_metric::mean_comp_ae;
    error_metric_table["luma_se"] = Error_metric::luma_se;
    error_metric_table["luma_ae"] = Error_metric::luma_ae;
    error_metric_table["ycbcr_vec_se"] = Error_metric::ycbcr_vec_se;
    error_metric_table["ycbcr_vec_ae"] = Error_metric::ycbcr_vec_ae;
    error_metric_table["vec_se"] = Error_metric::vec_se;
  }
  // Parse the error metric option
  if (error_metric_table.find(error_metric) != error_metric_table.end()) {
    error_metric_ = error_metric_table[error_metric];
  } else {
    std::cerr << "COMMAND ERROR: invalid error metric option: "
              << error_metric << '\n';
    return -1;
  }

  // Parse if to enable bad point removal option
  if (vm.count("bad-point-removal")) {
    enable_bad_point_removal_ = enable_bad_point_removal ? true : false;
  }

  // set is_removable_
  {
    int width = input_image_.width();
    int height = input_image_.height();
    is_removable_ = SPL::Array2<bool>(width, height, true);
    *(is_removable_.rowBegin(0) + 0) = false;
    *(is_removable_.rowBegin(0) + width - 1) = false;
    *(is_removable_.rowBegin(height - 1) + 0) = false;
    *(is_removable_.rowBegin(height - 1) + width - 1) = false;
  }

  return 0;
}

int Mesh_generator::make_mesh()
{
  log_message(1, "start to make initial mesh...\n");
  timer_.start();
  if (int status = make_initial_mesh(); status) {
    return status;
  }
  timer_.stop();
  log_message(1, "make initial mesh complete. Time: ", timer_.get(), '\n');

  log_message(1, "start to simplify mesh...\n");
  timer_.start();
  if (int status = simplify_mesh(); status) {
    return status;
  }
  timer_.stop();
  log_message(1, "simplify mesh complete.. Time: ", timer_.get(), '\n');

  if (enable_bad_point_removal_) {
    log_message(1, "start bad point removal...\n");
    postprocess_mesh();
    log_message(1, "bad point removal complete.\n");
  }

  return 0;
}

int Mesh_generator::make_initial_mesh()
{
  std::vector<Int_point> points;

  log_message(1, "selecting initial points...\n");
  // Use ED or some appropriate scheme to select points.
  if (int status = select_initial_points(points); status) {
    return status;
  }
  log_message(1, "select initial points complete.\n");
  log_message(1, "target initial point size = ", initial_size_,
              ", actual initial point size = ", points.size(), '\n');

  // Randomly permute the points prior to insertion for efficiency.
  unsigned seed = 0;
  random_shuffle(points.begin(), points.end(), seed);

  // Insert points in mesh.
  for (auto&& p : points) {
    Triangulation::Point point(p.x, p.y);
    tri_.insert(point);
  }

  invalidate_face_errs();

  log_message(1, "scan all faces and set the face error...\n");
  // Scan all faces and sum error to overall_error_
  for (auto face = tri_.faces_begin(); face != tri_.faces_end(); ++face) {
    double err = get_face_err(&*face);
    overall_error_ = combine_error(overall_error_, err);
  }
  log_message(1, "scan faces complete.\n");

  log_message(1, "configure vertex significance and the priority queue...\n");
  // Iterate through each vertex in the mesh
  Triangulation::Vertex_iterator vertIter = tri_.vertices_begin();
  for (int remaining_count = tri_.number_of_vertices(); remaining_count;) {
    // if the vertex is not removable
    if (Int_point p(vertIter->point().x(), vertIter->point().y());
        !is_removable(p)) {
      ++vertIter;
      --remaining_count;
      continue;
    }

    // if the vertex is already on the priority queue
    if (vertIter->on_pri_que()) {
      ++vertIter;
      continue;
    }

    Triangulation::Vertex_iterator next_vert_iter = std::next(vertIter, 1);

    // update the vertex's priority
    // WARNING: any references to this vertex will be invalidated
    update_vertex_priority(&*vertIter);

    vertIter = next_vert_iter;
    --remaining_count;
  }
  log_message(1, "configure vertex significance and the priority queue complete.\n");

  return 0;
}

int Mesh_generator::simplify_mesh()
{
  // if target size is specified and initial size equals to target size, then
  // return. this is to ensure no point removal happens in the ED method 
  // (due to the influence of size tolerance)
  if (target_size_ > 0 && initial_size_ == target_size_) {
    return 0;
  }

  while (1) {
    // If target size is specified and reached, break
    if (target_size_ > 0 && tri_.number_of_vertices() <= target_size_) {
      break;
    }

    // If tareget error is specified and reached, break
    if (target_error_ > 0 && overall_error_ <= target_error_) {
      break;
    }

    // If select point succeed
    if (Triangulation::Vertex_handle v; select_point_for_removal(v)) {
      log_message(2, "remove point ", v->point(), ", sig = ", v->get_sig(), '\n');
      remove_point(v);
    } else {
      return 1;
    }
  }

  return 0;
}

void Mesh_generator::postprocess_mesh()
{
  Triangulation::Vertex_handle v;
  while (select_point_for_removal(v)) {
    // if error increase, break
    if (v->get_sig() > 0) {
      break;
    } else {
      log_message(2, "remove point ", v->point(), ", sig = ", v->get_sig(), '\n');
      remove_point(v);
    }
  }
}

int Mesh_generator::output_mesh(std::ostream& out) const
{
#if 0
	mesh format is as following:
		WIDTH HEIGHT
		NUM_COMPONENTS PRECISION
		{OFF_DATA}

#endif

  // output header
  out << ori_input_image_.width() << " " << ori_input_image_.height() << '\n'
      << ori_input_image_.num_components() << " "
      << ori_input_image_.precision() << '\n';

  // if input image is grayscale image
  if (ori_input_image_.num_components() == 1) {
    auto get_vertex_z = [& image = this->ori_input_image_](
                         Triangulation::Vertex_const_handle v) 
    {
      Triangulation::Point p = v->point();
      return image[0](p.x(), p.y());
    };
    // output OFF data
    if (tri_.output_off(out, get_vertex_z)) {
      std::cerr << "write OFF data failed.\n";
      return 1;
    }

    // if input image is color image
  } else if (ori_input_image_.num_components() == 3) {
    auto get_vertex_z = [](Triangulation::Vertex_const_handle v) { return 0; };
    auto get_vertex_color = [& image = this->ori_input_image_](
      Triangulation::Vertex_const_handle v, std::vector<double> & color,
      bool& as_int)
    {
      Triangulation::Point p = v->point();
      color = { (double)image[0](p.x(), p.y()), (double)image[1](p.x(), p.y()),
                (double)image[2](p.x(), p.y()), 0 };
      as_int = true;
      return true;
    };

    bool provide_vertex_zs = true;
    bool provide_vertex_colors = true;

    // output OFF data
    if (tri_.output_off(out, get_vertex_z, nullptr, get_vertex_color, nullptr,
                        nullptr,
                        provide_vertex_zs | (provide_vertex_colors << 2))) {
      std::cerr << "write OFF data failed.\n";
      return 1;
    }
  }
  return 0;
}

int Mesh_generator::select_initial_points(std::vector<Int_point>& points)
{

  // if select all points on the sampling grid
  switch (initial_point_selection_policy_) {
    case Initial_point_selection_policy::all:
      select_initial_points_gpr(points);
      break;
    case Initial_point_selection_policy::random:
      select_initial_points_random(points);
      break;
    case Initial_point_selection_policy::uniform:
      select_initial_points_uniform(points);
      break;
    case Initial_point_selection_policy::ed:
      if (select_initial_points_ed(points)) {
        std::cerr
          << "select initial points failed due to size tolerance not being "
             "satisfied! The current error diffusion size tolerance is "
          << ed_size_tolerance_
          << ". Please increase it to a larger value by using the command "
             "line option: --ed-size-tolerance\n";
        return -1;
      }
      break;
    default:
      std::cerr << "invalid initial strategy option\n";
      return 1;
  }
  return 0;
}

void Mesh_generator::select_initial_points_gpr(std::vector<Int_point>& points)
{
  log_message(1, "select initial points policy: all\n");
  int width = input_image_.width();
  int height = input_image_.height();
  points.reserve(width * height);

  // select all points on the sampling grid
  for (int y = 0; y < height; ++y) {
    for (int x = 0; x < width; ++x) {
      points.push_back(Int_point(x, y));
    }
  }
}

void Mesh_generator::select_initial_points_random(std::vector<Int_point>& points)
{

  log_message(1, "select initial points policy: random\n");
  int width = input_image_.width();
  int height = input_image_.height();

  std::vector<Int_point> grid;
  grid.reserve(width * height);

  // push back points on the first and last row of the image grid
  // (except extreme convex hull points) seperately for efficiency
  for (int x = 1; x < width - 1; ++x) {
    grid.push_back(Int_point(x, 0));
    grid.push_back(Int_point(x, height - 1));
  }

  // push back points on the middle rows of the image grid
  for (int y = 1; y < height - 1; ++y) {
    for (int x = 0; x < width; ++x) {
      grid.push_back(Int_point(x, y));
    }
  }

  // NOTE: using a fixed seed here can make the shuffled elements
  // in grid having the same order every time we run the program.
  // This is helpful for testing purpose.
  unsigned int seed = 0;
  // random shuffle all grid points
  random_shuffle(grid.begin(), grid.end(), seed);

  // move the front (initial_size_ - 4) elements from grid into points
  std::move(grid.begin(), grid.begin() + initial_size_ - 4,
            std::back_inserter(points));

  // insert the extreme convex hull points of the image
  points.push_back(Int_point(0, 0));
  points.push_back(Int_point(0, height - 1));
  points.push_back(Int_point(width - 1, 0));
  points.push_back(Int_point(width - 1, height - 1));
}

void Mesh_generator::select_initial_points_uniform(std::vector<Int_point>& points)
{

  log_message(1, "select initial points policy: uniform\n");
  int width = input_image_.width();
  int height = input_image_.height();

  points.reserve(initial_size_);

  // get the sampling intervals
  double w = static_cast<double>(width) / (std::sqrt(initial_size_) - 1.0);
  double h = (static_cast<double>(height) / static_cast<double>(width)) * w;

  log_message(1, "width interval = ", w, '\n');
  log_message(1, "height interval = ", h, '\n');

  // push back points on the first row of the uniform grid
  // (except extreme convex hull points)
  for (double x = w; x < width - 1; x += w) {
    int x_value = SPL::clip<int>(rint(x), 0, width - 1);
    points.push_back(Int_point(x_value, 0));
  }

  // push back points on the middle rows of the uniform grid
  for (double y = h; y < height - 1; y += h) {
    int y_value = SPL::clip<int>(rint(y), 0, height - 1);
    if (y_value == height - 1) {
      continue;
    }
    for (double x = 0; x < width; x += w) {
      int x_value = SPL::clip<int>(rint(x), 0, width - 1);
      points.push_back(Int_point(x_value, y_value));
    }
  }

  // push back points on the last row of the uniform grid
  // (except extreme convex hull points)
  for (double x = w; x < width - 1; x += w) {
    int x_value = SPL::clip<int>(rint(x), 0, width - 1);
    points.push_back(Int_point(x_value, height - 1));
  }

  // insert the extreme convex hull points of the image
  points.push_back(Int_point(0, 0));
  points.push_back(Int_point(0, height - 1));
  points.push_back(Int_point(width - 1, 0));
  points.push_back(Int_point(width - 1, height - 1));
}

int Mesh_generator::select_initial_points_ed(std::vector<Int_point>& points)
{

  log_message(1, "select initial points policy: ed\n");
  int width = input_image_.width();
  int height = input_image_.height();

  SPL::Array2<int> ed_grid(width, height, 0);
  std::vector<SPL::Array2<double>> comps;
  for (int i = 0; i < input_image_.num_components(); ++i) {
    comps.push_back(input_image_[i]);
  }

  log_message(1, "selecting ed points...\n");
  try {
    if (make_mesh_ed(initial_size_, ed_size_tolerance_, comps, ed_grid)) {
      return -1;
    }
  } catch(const std::invalid_argument& ia) {
    std::cerr << "Invalid argument: " << ia.what() << '\n';
    return -1;
  }
  log_message(1, "select ed points complete.\n");

  // add extreme convex hull points. This will increase points size.
  ed_grid(0, 0) = 1;
  ed_grid(width - 1, 0) = 1;
  ed_grid(0, height - 1) = 1;
  ed_grid(width - 1, height - 1) = 1;

  // reserve size first for efficiency when doing push_back
  int size = std::count(ed_grid.begin(), ed_grid.end(), true);
  points.reserve(size);

  // convert ed_grid into points
  for (int y = 0; y < height; ++y) {
    for (int x = 0; x < width; ++x) {
      if (ed_grid(x, y)) {
        points.push_back(Int_point(x, y));
      }
    }
  }

  int seed = 0;
  // erase some elements in points until desire size reached
  while (points.size() > initial_size_ + ed_size_tolerance_) {
    // generate a random number between 0 ~ (points.size() - 1)
    int i = generate_random_number<double>(0.0, points.size() - 1, seed);

    // if points[i] is removable, erase it from points
    if (is_removable(points[i])) {
      points.erase(points.begin() + i);
    } else {
      // change the seed to make sure the same index won't be selected again
      // this is need to aviod dead loop
      ++seed;
    }
  }

  return 0;
}

int Mesh_generator::make_mesh_ed(int target_size, int size_tolerance,
                             const std::vector<SPL::Array2<double>>& components,
                             SPL::Array2<int>& result)
{

  int status = 0;
  const int width = input_image_.width();
  const int height = input_image_.height();

  // the density function to be used for error diffusion
  SPL::Array2<double> dens_func(width, height);

  // if grayscale image
  if (components.size() == 1) {

    // if related to laplacian of the image
    if ((ed_strategy_ == Ed_strategy::max_component_laplacian_ed) ||
        (ed_strategy_ == Ed_strategy::mean_component_laplacian_ed) ||
        (ed_strategy_ == Ed_strategy::gray_comp_laplacian_ed))
    {
      get_mag_laplacian_grayscale(components[0], derivative_bound_policy_,
        smooth_options_, use_double_convolution_, dens_func);
    } else {
      get_mmsodd_grayscale(components[0], derivative_bound_policy_,
        smooth_options_, use_double_convolution_, dens_func);
    }

    // normalize the density function
    dens_func /= dens_func.max();
    status = perform_fs_error_diffusion_by_size(dens_func, target_size, 
               size_tolerance, fs_ed_options_, result);

    // if color image
  } else if (components.size() == 3) {

    // If use vector_space_combination strategy
    if (ed_strategy_ == Ed_strategy::vector_space_combination) {
      std::array<SPL::Array2<double>, 3> comps{ { components[0], 
                                                  components[1],
                                                  components[2] } };
      get_mdg_color(comps, derivative_bound_policy_, smooth_options_,
                    use_double_convolution_, dens_func);
      dens_func /= dens_func.max();

      status = perform_fs_error_diffusion_by_size(
        dens_func, target_size, size_tolerance, fs_ed_options_, result);

      // if use three components union related strategy
    } else if (ed_strategy_ == Ed_strategy::three_components_union_mmsodd ||
               ed_strategy_ == Ed_strategy::three_components_union_laplacian)
    {
      std::array<SPL::Array2<double>, 3> ds;  // density func on each component
      double d_max = 0;
      // get the density function on each component
      for (int i = 0; i < 3; ++i) {
        if (ed_strategy_ == Ed_strategy::three_components_union_laplacian) {
          get_mag_laplacian_grayscale(components[i], derivative_bound_policy_,
             smooth_options_, use_double_convolution_, ds[i]);
        } else {
          get_mmsodd_grayscale(components[i], derivative_bound_policy_,
            smooth_options_, use_double_convolution_, ds[i]);
        }
        double max = ds[i].max();
        d_max = (max > d_max) ? max : d_max;
      }
      // normalize the density function
      for (auto&& dens_f : ds) {
        dens_f /= d_max;
      }

      status = perform_three_comps_union_ed_by_size(
        ds, target_size, size_tolerance, fs_ed_options_, result);

      // if use mean_component_mmsodd_ed or mean_component_laplacian_ed strategy
    } else if (ed_strategy_ == Ed_strategy::mean_component_mmsodd_ed ||
               ed_strategy_ == Ed_strategy::mean_component_laplacian_ed) {

      dens_func.fill(0.0);
      std::array<SPL::Array2<double>, 3> ds;  // density function on each component
      for (int i = 0; i < 3; ++i) {
        if (ed_strategy_ == Ed_strategy::mean_component_mmsodd_ed) {
          get_mmsodd_grayscale(components[i], derivative_bound_policy_,
            smooth_options_, use_double_convolution_, ds[i]);
        } else {
          get_mag_laplacian_grayscale(components[i], derivative_bound_policy_,
            smooth_options_, use_double_convolution_, ds[i]);
        }
        dens_func += ds[i];
      }
      dens_func /= dens_func.max();

      status = perform_fs_error_diffusion_by_size(
        dens_func, target_size, size_tolerance, fs_ed_options_, result);

      // if use max_component_mmsodd_ed or max_component_laplacian_ed strategy
    } else if (ed_strategy_ == Ed_strategy::max_component_mmsodd_ed ||
               ed_strategy_ == Ed_strategy::max_component_laplacian_ed) {

      std::array<SPL::Array2<double>, 3> ds;
      for (int i = 0; i < 3; ++i) {
        if (ed_strategy_ == Ed_strategy::max_component_mmsodd_ed) {
          get_mmsodd_grayscale(components[i], derivative_bound_policy_,
            smooth_options_, use_double_convolution_, ds[i]);
        } else {
          get_mag_laplacian_grayscale(components[i], derivative_bound_policy_,
            smooth_options_, use_double_convolution_, ds[i]);
        }
      }

      // for each point, get its maximum value among all components
      for (int y = 0; y < height; ++y) {
        for (int x = 0; x < width; ++x) {
          double max = ds[0](x, y);
          for (int i = 1; i < 3; ++i) {
            max = (ds[i](x, y) > max ? ds[i](x, y) : max);
          }
          dens_func(x, y) = max;
        }
      }
      dens_func /= dens_func.max();

      status = perform_fs_error_diffusion_by_size(
        dens_func, target_size, size_tolerance, fs_ed_options_, result);

    } else {
      std::cerr << "Invalid ed strategy\n";
      return -1;
    }
  }

  return status;
}

int Mesh_generator::perform_three_comps_union_ed_by_size(
  const std::array<SPL::Array2<double>, 3>& density_funcs, int num_points,
  int tolerance, const Fs_error_diffusion_options& options,
  SPL::Array2<int>& result)
{
  const std::array<SPL::Array2<double>, 3>& ds = density_funcs;

  const int width = ds[0].getWidth();
  const int height = ds[0].getHeight();

  tolerance = std::abs(tolerance);
  int lower_bound = num_points - tolerance;
  int upper_bound = num_points + tolerance;
  const int n = num_points;

  // the startup error on each component
  std::array<std::vector<double>, 3> initial_errors;

  // set startup errors in each component
  switch (options.initial_error_mode_) {
    case Initial_error_mode::zero: {
      for (int i = 0; i < 3; ++i) {
        initial_errors[i] = std::vector<double>(width, 0.0);
      }
    } break;
    case Initial_error_mode::random:
      // range -2t to 1.0 - 2t
      for (int i = 0; i < 3; ++i) {
        double t = std::accumulate(ds[i].begin(), ds[i].end(), 0.0) / (2.0 * n / 3);
        initial_errors[i].resize(width);
        for (int x = 0; x < width; ++x) {
          double r = generate_random_number<double>(-2.0 * t, 1.0 - 2.0 * t, 0);
          initial_errors[i][x] = r;
        }
      }
      break;
    case Initial_error_mode::mirror: {
      for (int i = 0; i < 3; ++i) {
        SPL::Array2<int> b;
        // the errors diffused to the last row of each component
        std::array<std::vector<double>, 3> output_errors;
        SPL::Array2<double> mirror_d = ds[i];
        mirror_d.flipud();
        double t = std::accumulate(mirror_d.begin(), mirror_d.end(), 0.0) /
                   (2.0 * n / 3);
        perform_fs_error_diffusion_by_threshold(
          mirror_d, t, options, initial_errors[i], output_errors[i], b);
        initial_errors[i] = output_errors[i];
      }
    } break;
    default: {
      std::cerr << "invalid initial error mode\n";
      return -1;
    }
  }

  double t = std::accumulate(ds[0].begin(), ds[0].end(), 0.0) / (2.0 * n);

  int size = perform_three_comps_union_ed_by_threshold(ds, initial_errors, 
             options, t, result);

  double t_high = 0;
  double t_low = 0;
  // Find a t_high and t_low
  if (size > upper_bound) {
    while (size > upper_bound) {
      t_low = t;
      t *= 2;
      size = perform_three_comps_union_ed_by_threshold(ds, initial_errors,
             options, t, result);
    }
    t_high = t;
  } else if (size < lower_bound) {
    while (size < lower_bound) {
      t_high = t;
      t /= 2;
      size = perform_three_comps_union_ed_by_threshold(ds, initial_errors,
             options, t, result);
    }
    t_low = t;
  } else {
    t_low = t;
    t_high = t;
  }

  log_message(2, "t_high = ", t_high, ", t_low = ", t_low, '\n');
  int min_size_diff = std::abs(size - n);

  // perform the binary search until optimal size is achieved
  while (t_low <= t_high) {
    t = (t_low + t_high) / 2.0;
    size = perform_three_comps_union_ed_by_threshold(ds, initial_errors,
           options, t, result);

    if (size > upper_bound) {
      t_low = t;
    } else if (size < lower_bound) {
      t_high = t;
    } else {
      break;
    }

    // this is for giving a recommend size tolerance when binary search failed
    if (int diff = std::abs(size - n); diff < min_size_diff) {
      min_size_diff = diff;
    }

    if ((t_low + t_high) / 2.0 == t) {
      std::cerr << "binary search failed! recommend minimum ed size tolerance: "
                << min_size_diff << "\n";
      return -1;
    }
  }

  return 0;
}

int Mesh_generator::perform_three_comps_union_ed_by_threshold(
    const std::array<SPL::Array2<double>, 3>& density_funcs,
    const std::array<std::vector<double>, 3>& startup_errors,
    const Fs_error_diffusion_options& options, double threshold,
    SPL::Array2<int>& result_bitmap)
{

  const int width = density_funcs[0].getWidth();
  const int height = density_funcs[0].getHeight();
  std::array<SPL::Array2<int>, 3> bitmap_comps;
  std::array<std::vector<double>, 3> output_errors;
  int count = 0; // count the number of points selected

  // for each component, perform fs-ed by threshold 
  for (int i = 0; i < 3; ++i) {
    perform_fs_error_diffusion_by_threshold(density_funcs[i], threshold, 
      options, startup_errors[i], output_errors[i], bitmap_comps[i]);
  }

  // union the bitmap on all components
  for (int y = 0; y < height; ++y) {
    for (int x = 0; x < width; ++x) {
      result_bitmap(x, y) =
        bitmap_comps[0](x, y) || bitmap_comps[1](x, y) || bitmap_comps[2](x, y);

      if(result_bitmap(x, y)) {
        ++count;
      }
    }
  }

  // return the number of points selected in the bitmap
  return count;
}


bool Mesh_generator::select_point_for_removal(Triangulation::Vertex_handle& vertex)
{
  if (vertex_queue_.empty()) {
    return false;
  } else {
    vertex = vertex_queue_.top();
    return true;
  }
}

void Mesh_generator::remove_point(Triangulation::Vertex_handle vertex)
{
#if 0
 pseudocode:
	record 1-ring neighbours
	delete vertex
	for each affected vertex
		update_vertex_priority
#endif

  // record 1-ring affected vertices
  std::vector<Triangulation::Vertex_handle> affected_vertices;
  Triangulation::Halfedge_around_vertex_circulator h_circ =
    vertex->vertex_begin();
  do {
    auto v = h_circ->opposite()->vertex();
    if (v->on_pri_que()) {
      affected_vertices.push_back(v);
    }
  } while (++h_circ != vertex->vertex_begin());

  auto pri_que_pos = vertex->get_pri_que_pos();
  vertex->clear_pri_que_pos();
  // erase the vertex from the priority queue
  vertex_queue_.erase(pri_que_pos);

  // update overall error
  overall_error_ = combine_error(overall_error_, vertex->get_sig());

  // remove the vertex from the mesh
  tri_.remove(vertex->halfedge());

  // for each affected vertex, update its priority
  for (auto&& v : affected_vertices) {
    log_message(2, "update priority of affected vertex ", v->point(), '\n');
    update_vertex_priority(v);
  }
}

double Mesh_generator::compute_mesh_error() const
{
  double sum_err = 0.0;
  for (auto face = tri_.faces_begin(); face != tri_.faces_end(); ++face) {
    assert(face->is_err_valid());
    sum_err = combine_error(sum_err, face->get_err());
  }

  return sum_err;
}

double Mesh_generator::get_face_err(Triangulation::Face_handle face)
{
  double err;
  if (face->is_err_valid()) {
    err = face->get_err();
  } else {
    err = calc_face_err(face);
    face->set_err(err);
  }

  return err;
}

double Mesh_generator::calc_face_err(Triangulation::Face_handle face)
{

  auto h = face->halfedge();
  std::array<Triangulation::Point, 3> points{
    { h->vertex()->point(), h->next()->vertex()->point(),
      h->prev()->vertex()->point() }
  };

  Scan_stats stats;
  bool write_data_to_buffer = false;
  bool write_stats = true;

  // set the write_statistics flag on
  int flags = Scan_triangle_options::write_statistics;
  Image recon_image;

  // if the input is a grayscale image
  if (input_image_.num_components() == 1) {
    // function value of each vertex
    std::array<int, 3> fs{ { input_image_[0](points[0].x(), points[0].y()),
                             input_image_[0](points[1].x(), points[1].y()),
                             input_image_[0](points[2].x(), points[2].y()) } };

    scan_triangle_grayscale(points[0].x(), points[0].y(), points[1].x(),
                            points[1].y(), points[2].x(), points[2].y(), fs[0],
                            fs[1], fs[2], input_image_, error_metric_,
                            recon_image, stats, flags);

    // if the input is a rgb image
  } else {
    // function values of the first vertex
    pixel<int, 3> a_fs{ input_image_[0](points[0].x(), points[0].y()),
                        input_image_[1](points[0].x(), points[0].y()),
                        input_image_[2](points[0].x(), points[0].y()) };
    // function values of the second vertex
    pixel<int, 3> b_fs{ input_image_[0](points[1].x(), points[1].y()),
                        input_image_[1](points[1].x(), points[1].y()),
                        input_image_[2](points[1].x(), points[1].y()) };
    // function values of the third vertex
    pixel<int, 3> c_fs{ input_image_[0](points[2].x(), points[2].y()),
                        input_image_[1](points[2].x(), points[2].y()),
                        input_image_[2](points[2].x(), points[2].y()) };

    scan_triangle_color(points[0].x(), points[0].y(), points[1].x(),
                        points[1].y(), points[2].x(), points[2].y(), a_fs, b_fs,
                        c_fs, input_image_, error_metric_, recon_image, stats,
                        flags);
  }

  return stats.error_;
}

void Mesh_generator::invalidate_face_errs()
{
  for (auto face = tri_.faces_begin(); face != tri_.faces_end(); ++face) {
    face->invalidate_err();
  }
  overall_error_ = 0;
}

void Mesh_generator::update_vertex_priority(Triangulation::Vertex_handle vertex)
{
#if 0
pseudocode:
	record edges incident on hole to be formed
	save old errors + other states
	remove from priority queue
	delete vertex
	scan convert new faces
	compute new error
	add vertex back
#endif

  double err0 = 0;
  double err1 = 0;
  std::set<Triangulation::Halfedge_handle> hole_edges;
  std::vector<std::pair<Triangulation::Halfedge_handle, double>>
    cached_face_errs;

  // if the vertex is a border vertex, this is used to track the border halfedge
  Triangulation::Halfedge_handle border_edge = nullptr;
  Triangulation::Halfedge_around_vertex_circulator h_circ =
    vertex->vertex_begin();
  do {
    // if this is a border halfedge, record its prev() halfedge
    if (h_circ->is_border()) {
      border_edge = h_circ->prev();
      continue;
    }

    // record the hole halfedge
    hole_edges.insert(h_circ->prev());
    // get its incident face error
    double err = get_face_err(h_circ->face());
    err0 = combine_error(err0, err);

    // record the hole halfedge and its incident face error
    // as a pair into cached_face_errs
    cached_face_errs.push_back(std::make_pair(h_circ->prev(), err));

  } while (++h_circ != vertex->vertex_begin());

  // if the vertex is already on the priority queue, remove it from
  // the priority queue
  if (vertex->on_pri_que()) {
    Vertex_priority_queue::handle_type pri_que_pos = vertex->get_pri_que_pos();
    vertex->clear_pri_que_pos();
    vertex_queue_.erase(pri_que_pos);
  }

  Triangulation::Point p = vertex->point();
  // record a stable halfedge, used to find the vertex after point been added
  // back
  Triangulation::Halfedge_handle stable_halfedge = vertex->halfedge()->prev();

  // temporary remove the vertex from the mesh
  tri_.remove(vertex->halfedge());

  log_message(2, "point been removed from the triangulation: ", p, '\n');

  // if border edge exists, insert the new halfedge to hole_edges
  // this ensures the hole is enclosed
  if (border_edge) {
    hole_edges.insert(border_edge->next()->opposite());
  }

  auto update_face_stat = [&err1, this](Triangulation::Face_handle face) {
    double err = this->get_face_err(face);
    err1 = this->combine_error(err1, err);
  };

  // computes the sum of new face errors within the hole
  for_each_face_in_region<Triangulation, decltype(update_face_stat)>(
    tri_, hole_edges, update_face_stat);

  // add the vertex back to the mesh
  tri_.insert(p);

  // restore all the face errors from cached_face_errs
  for (int i = cached_face_errs.size() - 1; i >= 0; --i) {
    auto face = (cached_face_errs[i].first)->face();
    face->invalidate_err();
    face->set_err(cached_face_errs[i].second);
  }

  // get the vertex handle.
  // Note: this is not the original vertex, though it referes to the same point
  Triangulation::Vertex_handle new_vertex = stable_halfedge->next()->vertex();

  // set significance of the vertex
  new_vertex->set_sig(err1 - err0);

  // push the new vertex into the priority queue
  new_vertex->clear_pri_que_pos();
  new_vertex->set_pri_que_pos(vertex_queue_.push(new_vertex));
}

bool Mesh_generator::is_removable(const Int_point& p) const
{
  return *(is_removable_.rowBegin(p.y) + p.x);
}

double Mesh_generator::combine_error(double a, double b) const
{
  return a + b;
}
