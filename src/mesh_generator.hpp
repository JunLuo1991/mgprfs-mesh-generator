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


#ifndef mesh_generator_hpp
#define mesh_generator_hpp

#include "error_diffusion.hpp"
#include "image.hpp"
#include "scan_triangle.hpp"
#include "util.hpp"
#include <CGAL/Cartesian.h>
#include <CGAL/Filtered_kernel.h>
#include <SPL/Array2.hpp>
#include <SPL/Timer.hpp>
#include <SPL/cgal/Delaunay_triangulation_2.hpp>
#include <SPL/cgal/Triangulation_hierarchy_2.hpp>
#include <algorithm>
#include <boost/heap/fibonacci_heap.hpp>
#include <iostream>
#include <random>

/*
@brief A structure representating points in 2D space.
*/
struct Int_point
{
  Int_point(int x_value, int y_value)
    : x(x_value)
    , y(y_value)
  {
  }
  int x;
  int y;
};

class Mesh_vertex;
class Mesh_face;

struct Mesh_traits
{
  using Geom_kernel = CGAL::Filtered_kernel<CGAL::Cartesian<double>>;
  using Vertex = Mesh_vertex;
  using Edge = SPL::Triangulation_edge_base_2<Mesh_traits>;
  using Face = Mesh_face;
  using Halfedge = SPL::Triangulation_halfedge_base_2<Mesh_traits>;
};

/*!
 * @brief A functor to compare the priority of the priority queue.
 */
struct Vertex_priority_queue_compare
{
  bool operator()(const Mesh_vertex*, const Mesh_vertex*) const;
};

/*
 @brief A class representing the mesh vertex
*/
class Mesh_vertex
  : public SPL::Triangulation_hierarchy_vertex_base_2<Mesh_traits>
{
public:

  // The vertext priority queue
  using Vertex_priority_queue = boost::heap::fibonacci_heap<
    Mesh_vertex*, boost::heap::compare<Vertex_priority_queue_compare>,
    boost::heap::stable<true>>;

  /*
   @brief The default constructor
   */
  Mesh_vertex();

  /*
   @brief Get the significance of the vertex.
   @detail
   Get the significance of the vertex with respect to deletion.
   i.e., the increase in triangulation cost if this vertex is removed.
   */
  double get_sig() const;

  /*
   @brief Set the significance of the vertex.
   */
  void set_sig(double sig);

  /*
   @brief Check if the vertex is on the vertex priority queue or not.
   */
  bool on_pri_que() const;

  /*
   @brief Get the priority queue position associated with the vertex.
   @pre The vertex must be on the queue.
   */
  Vertex_priority_queue::handle_type get_pri_que_pos() const;

  /*
   @brief Set the priority queue position.
   @pre The vertex must not be on the queue.
  */
  void set_pri_que_pos(Vertex_priority_queue::handle_type pos);

  /*
   @brief Clear the priority queue position.
   @post The vertex is marked as no longer being on the queue.
   */
  void clear_pri_que_pos();

private:
  //! Flag indicating if this vertex is on the priority queue.
  bool on_pri_que_;

  //! The entry on the priority queue associated with this vertex.
  Vertex_priority_queue::handle_type pri_que_pos_;

  //! The significance (with respect to deletion) of the vertex.
  // (i.e., the amount by which the triangulation cost increases
  // when the vertex is deleted).
  // This value is used to determine the priority of the vertex
  // on the vertex priority queue.
  // Changing this value must always be accompanied by a change to
  // the priority queue.
  double sig_;
};

/*
 @brief A class representing the mesh face
*/
class Mesh_face : public SPL::Triangulation_face_base_2<Mesh_traits>
{
public:
  Mesh_face();

  /*
   @brief Get the error for this face.
   @pre The error must be valid.
   */
  double get_err() const;

  /*
   @brief Set the error for this face.
   @pre The error must be invalid?
   */
  void set_err(double err);

  /* 
   @brief Is the error valid or not.
   */
  bool is_err_valid() const;

  /*
   @brief Mark the current error as invalid.
   */
  void invalidate_err();

private:
  //! Is the cached error valid?
  bool err_valid_;

  //! The error currently associated with this face.
  double err_;
};


/*
 @brief A functor representating the mesh generator
 */
class Mesh_generator
{
public:
  /*
   @brief Default constructor
   */
  Mesh_generator();

  // Not copyable.
  Mesh_generator(const Mesh_generator&) = delete;
  Mesh_generator& operator=(const Mesh_generator&) = delete;

  // Not movable.
  Mesh_generator(Mesh_generator&&) = delete;
  Mesh_generator& operator=(Mesh_generator&&) = delete;

  /*
   @brief Generate a tirangle mesh based on specified image and options.
   @param in  The stream which specifies the image data to read from.
   @param out  The stream which specifies the mesh model to write to.
   @param argc  Number of command-line arguments
   @param argv  Array of character pointers listing all the arguments.
   @return
   upon success, return zero; otherwise, a non-zero value is returned.
   */
  int operator()(std::istream& in, std::ostream& out, int argc,
                 const char* const* argv);

private:
  /*!
  @brief Enumerate the method to generate a color image mesh model
  @detail
  This lists all avaliable methods for color image mesh genration.

  multi_comp_mesh_generation: process multiple components together by
  using the color image mesh generator.
  luma_space_mesh_generator: convert the color image into luma space first,
  then call the grayscale image mesh generator to get a mesh, finally 
  replace all sample point values with color values in the mesh.
  */
  enum class Mesh_generation_policy
  {
    multi_comp_mesh_generation,
    luma_space_mesh_generation
  };

  /*!
  @brief The initial point selection policy.
  */
  enum class Initial_point_selection_policy
  {
    all,
    random,
    uniform,
    ed
  };

  /*!
  @brief The error diffusion strategy for selecting initial sample points.
  @details
  1. gray_comp_mmsodd_ed: convert the image to a grayscale image, use mmsodd
  of the converted image to perform Floyd-steinberg error diffusion, 
  eventually replace all gray value to rgb value for each mesh point.
  2. gray_comp_laplacian_ed: convert the image to a grayscale image, use
  laplacian of the converted image to perform Floyd-steinberg error diffusion, 
  eventually replace all gray value to rgb value for each mesh point.
  3. vector_space_combination: compute the vector gradient based on Di Zenzo's
  method as the density function, and perform Floyd-steinberg error diffusion.
  4. max_components_mmsodd_ed: compute mmsodd on each channel, then pixel-wisely
  take the maximum value and compose a density function, then perform 
  Floyd-steinberg error diffusion
  5. max_components_laplacian_ed: compute laplacian on each channel, then
  pixel-wisely take the maximum value and compose a density function, then 
  perform Floyd-steinberg error diffusion.
  6. mean_components_mmsodd_ed: compute mmsodd on each channel, then 
  pixel-wisely take the mean value and compose a density function, then
  perform Floyd-steinberg error diffusion.
  7. mean_components_laplacian_ed: compute laplacian on each channel, then 
  pixel-wisely take the mean value and compose a density function, then
  perform Floyd-steinberg error diffusion.
  8. three_components_union_mmsodd: on each component take mmsodd as density
  function, perform Floyd-steinberg error diffusion, and get the union set.
  By applying binary search method, the mesh will have desired (or closed
  to desired) number of points.
  9. three_components_union_laplacian: on each component take laplacian as
  density function, perform Floyd-steinberg error diffusion, and get the 
  union set. By applying binary search method, the mesh will have desired 
  (or closed to desired) number of points.
  */
  enum class Ed_strategy
  {
    gray_comp_mmsodd_ed,
    gray_comp_laplacian_ed,
    vector_space_combination,
    max_component_mmsodd_ed,
    max_component_laplacian_ed,
    mean_component_mmsodd_ed,
    mean_component_laplacian_ed,
    three_components_union_mmsodd,
    three_components_union_laplacian,
  };

#if defined(USE_TRIANGULATION_HIERARCHY)
  using Triangulation =
    SPL::Triangulation_hierarchy_2<SPL::Delaunay_triangulation_2<Mesh_traits>>;
#else
  using Triangulation = SPL::Delaunay_triangulation_2<Mesh_traits>;
#endif

  using Vertex_priority_queue = Mesh_vertex::Vertex_priority_queue;

  /*!
  @brief Clear the state of the mesh generator.
  @detail
  This function is needed in case the user of the class
  tries to generate multiple meshes with the same mesh-generator
  object.
  */
  void clear();

  /*!
  @brief Generate a mesh.
  @param in  The stream which specifies the image data to read from.
  @param out  The stream which specifies the mesh model to write to.
  @param argc  Number of command-line arguments
  @param argv  Array of character pointers listing all the arguments.
  @return
  upon success, return zero; otherwise, a non-zero value is returned.
  @detail
  This function is the main entry point for mesh generation.
  */
  int main(std::istream& in, std::ostream& out, int argc,
           const char* const* argv);

  /*!
  @brief Process the mesh-generation options and initialize the
  state of the mesh generator based on these options.
  @param in  The stream which specifies the image data to read from.
  @param out  The stream which specifies the mesh model to write to.
  @param argc  Number of command-line arguments.
  @param argv  Array of character pointers listing all the arguments.
  @return
  Upon success, return zero; otherwise, a non-zero value is returned.
  */
  int process_options(std::istream& in, std::ostream& out, int argc,
                      const char* const* argv);

  /*!
  @brief Generate mesh.
  @detail
  This function actually performs mesh generation.
  @pre
  The mesh generator state is assumed to already have been initialized.
   */
  int make_mesh();

  /*!
  @brief Generate the initial mesh.
  @detail
  This function generates the initial mesh used for mesh generation.
  Several strategies can be used to select the initial mesh.
  For GPR, the entire sampling grid is used.
  For GPRFS-ED, error diffusion is used.
   */
  int make_initial_mesh();

  /*!
  @brief Simplify the mesh (by removing vertices).
  @detail
  This function implements the mesh simplification process.
  Vertices are removed from the mesh until the target mesh size
  is acheived or an error threshold would be violated.
   */
  int simplify_mesh();

  /*!
  @brief Postprocess the mesh (with bad-point removal).
  @detail
  This function applies bad-point removal (BPR).
   */
  void postprocess_mesh();

  /*!
  @brief Output current mesh to specified stream.
  @detail
  This function outputs the mesh model to the specified stream.
  The format is header plus OFF data.
   */
  int output_mesh(std::ostream& out) const;

  /*!
  @brief Select the sampling points for the initial mesh.
  @param points  The resulting selected points.
  @detail
  This function selects the sampling points for the initial mesh.
  This can be done in one of several ways.
  If the initial mesh density is 100%, all of the points on the
  sampling grid should be chosen.
  If the initial mesh density is less than 100%, the points should
  be chosen using error diffusion or some other method
  (such as random or uniform).
  */
  int select_initial_points(std::vector<Int_point>& points);

  /*!
  @brief
  Use the greedy point removal method for selecting the initial mesh points.
  @param points  The resulting selected points.
  @detail
  This funciton will select all sampling points on the sampling grid as
  the initial mesh sampling points.
  */
  void select_initial_points_gpr(std::vector<Int_point>& points);

  /*!
  @brief
  Randomly select points on the sampling grid as the intial mesh points.
  @param points  The resulting selected points.
  */
  void select_initial_points_random(std::vector<Int_point>& points);

  /*!
  @brief
  Uniformly select points on the sampling grid as the intial mesh points.
  @param points  The resulting selected points.
  */
  void select_initial_points_uniform(std::vector<Int_point>& points);

  /*!
  @brief
  Use error diffusion method to select the intial mesh points.
  @param points  The resulting selected points.
  @details
  This will call make_mesh_ed to get the error diffusion bitmap
  result. It will then add extreme convex hull points into the bitmap,
  and finally convert the bitmap into points.
  */
  int select_initial_points_ed(std::vector<Int_point>& points);

  /*!
  @brief 
  Perform the error diffusion method based on number of components
  and generate a bitmap result.
  @param target_size  The target mesh size.
  @param size_tolerance  The tolerance of the result size.
  @param components  The image components.
  @param result  A 2d array to store the bitmap result.
  @return
  Upon success, return zero; otherwise, a non-zero value is returned.

  @details
  This function implements the user defined error diffusion method.
  It will smooth the components, get the density function and 
  Floyd-Steinberg error diffusion based on different cases.

  In the bitmap result, 1 means the point is selected, 0 means not
  being selected.
   */
  int make_mesh_ed(int target_size, int size_tolerance,
                   const std::vector<SPL::Array2<double>>& components,
                   SPL::Array2<int>& result);

  /*!
  @brief
  Perform the three_components_union error diffusion method by a given size.

  @param density_funcs  The sample point density function on each component.
  @param num_points  The desired number of sampling points in the result bitmap.
  @param tolerance  The tolerance in terms of the num_points.
  @param options
  The options used for the Floyg-steinberg error diffusion method.
  @param result  A 2d array to store the bitmap result.
  @return 
  Upon success, return zero; otherwise, a non-zero value is returned.

  @details
  This function will run the Floy-steinberg error diffusion method on each of
  the image components and set the union as the result bitmap.

  For a point (x, y), if any of the result bitmap in R, G, B
  components is one, result(x,y) will be set as one, otherwise set as zero.

  The result bitmap will have a desired (or close to desired) number of points.
  Internally it will do binary search to find an optimal threshold first and
  perform error diffusion. The result size will be in the range:
  (num_points - tolerance, num_points + tolerance).

  The function can fail in the case that the desired size is not achievable
  while doing binary search, a non-zero integer will be returned in this case.
  */
  int perform_three_comps_union_ed_by_size(
    const std::array<SPL::Array2<double>, 3>& density_funcs, int num_points,
    int tolerance, const Fs_error_diffusion_options& options,
     SPL::Array2<int>& result_bitmap);

  /*!
  @brief
  Perform the three_components_union error diffusion method by a threshold.

  @param density_funcs  The sample point density function on each component.
  @param startup_errors  The startup errors on each component
  @param options
  The options used for the Floyg-steinberg error diffusion method.
  @param threshold  The threshold of error diffusion.
  @param result_bitmap  A 2d array to store the bitmap result.
  @return  The number of points selected in result_bitmap.

  @details
  This function will run the Floy-steinberg error diffusion method on each of
  the image components to get three sub bitmaps, and then set the union as the 
  result bitmap.

  For a point (x, y), if any of the sub bitmap in R, G, B components is one,
  result_bitmap(x,y) will be set as one, otherwise set as zero.
  */
  int perform_three_comps_union_ed_by_threshold(
    const std::array<SPL::Array2<double>, 3>& density_funcs,
    const std::array<std::vector<double>, 3>& startup_errors,
    const Fs_error_diffusion_options& options, double threshold,
     SPL::Array2<int>& result_bitmap);

  /*!
  @brief Select the next point for removal from the mesh during
  mesh simplification.
  @param vertex  The result selected vertex.
  @detail
  This function selects the next point for removal during mesh
  simplification (but does not remove it).
  If a point is successfully selected, @c vertex is set to the
  selected point and true is returned.
  If point selection fails (i.e., due to no removable points being
  left), false is returned.
  */
  bool select_point_for_removal(Triangulation::Vertex_handle& vertex);

  /*!
  @brief Remove the specified vertex from the mesh.
  @detail
  This function removes the specified point from the mesh.
  The significance values of the affected vertices are updated.
  WARNING: This will as a side effect invalidate numerous vertex handles.
  */
  void remove_point(Triangulation::Vertex_handle vertex);

  /*!
  @brief Compute the error associated with the mesh model.
  @detail
  This function is extremely slow and it only intended for testing
  and debugging purposes.
  */
  double compute_mesh_error() const;

  /*!
  @brief Rasterize a face and save various resulting statistics in the
  face itself.
  @detail
  If the face has a cache error, it will return the cache error directly.
  Otherwise this function will call calc_face_err to rasterize and
  calculate the face error, and save the error value in the face itself.
  */
  double get_face_err(Triangulation::Face_handle face);

  /*!
  @brief Rasterize a face and return its face error.
  @detail
  This function uses scan_triangle to do most of its work.
  */
  double calc_face_err(Triangulation::Face_handle face);

  /*!
  @brief Invalidate the entire face-error cache.
  */
  void invalidate_face_errs();

  /*!
  @brief Update the priority of a particular vertex.
  @detail
  This function updates the priority of a vertex whose priority
  value may be out of date.
  WARNING!  WARNING!  WARNING!
  This function deletes the specified vertex (although the point for
  the vertex is reinserted afterwards).
  Any references to this vertex will be invalidated.
  Be extremely careful how this function is used!
  */
  void update_vertex_priority(Triangulation::Vertex_handle vertex);

  /*!
  @brief Test is a point is allowed to be removed from the mesh.
  @param p  The input point
  @detail
  If the point is allowed to be removed from the mesh, true is returned;
  otherwise, false is returned.
  The extreme convex hull points should never be removed, since this
  will cause many unnecessary complications to the mesh-generation
  method.
  */
  bool is_removable(const Int_point& p) const;

  /*!
  @brief Combine two error values.
  @detail
  This function would be useful for supporting non-additive errors
  such as an infinity norm.
  For an additive error, returns a+b.
  For infinity norm, returns max(a, b).
  */
  double combine_error(double a, double b) const;

  /*!
  @brief Generate a debugging message of the specified severity.
  
  @param severity  The severity of the message.
  @param args  The data to be output.

  @details
  If the severity @c severity is less than or equal to the
  debug level setting for the mesh generator, the specified
  message is output to the standard error stream.  Otherwise,
  this function does nothing.
  (Essentially, this function discards any error message
  deemed not to be sufficiently important to print.)
  */
  template <class... T>
  void log_message(int severity, const T&... args) const
  {
    if (severity <= debug_level_) {
      // fold expression (c++17)
      ((std::cerr << args), ...);
    }
  }

  /*!
  @brief Generate a history message of the specified level.

  @param level  The level of the message.
  @param args  The data to be output.

  @details
  If the level @c level is less than or equal to the history level
  setting of the mesh generator, the specified message is output
  to the mesh-generator's history stream.  Otherwise, this function
  does nothing.
  */
  template <class... T>
  void history_message(int level, const T&... args) const
  {
    if (level <= history_level_) {
      // fold expression (c++17)
      ((*history_stream_ << args), ...);
    }
  }

  /*!
  @brief
  Is the mesh generator data structure appears to be sane?
  @detail
  This function is intended for debugging.
  */
  bool is_sane() const;

  //! Mesh generation policy
  Mesh_generation_policy mesh_generation_policy_;

  //! The original input image.
  // NOTE: This is used as a backup of the input image.
  // This is added because for color image case for some Ed_strategy, we need 
  // to convert the input image into grayscale first.
  Image ori_input_image_;

  //! The image for which a mesh model is to be generated.
  // NOTE: The mesh generator will internally use this image object to generate
  // the mesh model, not ori_input_image_. But the finaly mesh points function
  // values are outputed based on ori_input_image_.
  // When the input image is a color image and the the Ed_strategy is
  // gray_comp_laplacian_ed or gray_comp_mmsodd_ed, this object will be set
  // as the converted grayscale image that converted from ori_input_image_. 
  // In any other cases, it will be a copy of ori_input_image_;
  Image input_image_;

  //! The target size (in vertices) for the mesh to be generated.
  int target_size_;

  //! The target error criterion.
  // If negative, no error criterion.
  double target_error_;

  //! The size of the initial mesh (in samples).
  int initial_size_;

  //! The level of debugging messages.
  // If zero, no debugging message are generated.
  int debug_level_;

  //! The level of history message generation.
  // If zero, no history messages are generated.
  int history_level_;

  //! The stream to which history messages should be written.
  std::unique_ptr<std::ostream> history_stream_;

  //! The triangulation used to represent the mesh.
  Triangulation tri_;

  //! The vertex priority queue.
  Vertex_priority_queue vertex_queue_;

  // Array indicating if every grid point is removable from the mesh or not.
  // ture means yes, false means no.
  SPL::Array2<bool> is_removable_;

  std::vector<double> per_comp_error_;
  double overall_error_;

  // The error metric
  Error_metric error_metric_;

  // Enable bad point removal or not
  bool enable_bad_point_removal_;

  // Initial point selection strategy
  Initial_point_selection_policy initial_point_selection_policy_;

  //! The size tolerance for the error diffusion method.
  int ed_size_tolerance_;
  //! The smooth options
  Smooth_options smooth_options_;
  //! The derivative boundary handling policy.
  int derivative_bound_policy_;
  //! A bool value to indicate if do smooth and derivative together
  bool use_double_convolution_;
  //! The error diffusion strategy.
  Ed_strategy ed_strategy_;
  //! The Floyd-steinberg error diffusion options.
  Fs_error_diffusion_options fs_ed_options_;

  //! A timer used to track execution time
  SPL::Timer timer_;
};

#include "mesh_generator_impl.hpp"

#endif
