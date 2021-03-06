
# The ``make_model`` Program

## Synopsis

``make_model [options] <input_file> output_file``

## Description

This program reads an image in PNM/PGM/PPM format from
standard input, generates a mesh model in MODEL (Header + OFF) format,
and writes the model to standard output. The ``make_model`` program
consists of three main functional blocks: the initial mesh generation
block, the mesh simplification block, and a postprocessing block.
The initial mesh-generation block will select the initial mesh points
based on the specified initial generator (i.e., ``all``, ``random``,
``ed``). The mesh simplification block will remove sample points from
the mesh until the desired mesh size is achieve. The postprocessing
block provides the ability to keep removing points if this will reduce
the overall approximation error. In order to run the ``make_model``
program, an image in PNM/PPM/PGM format must be specified from standard
input, and the desired sampling density (or size, or error) must be
specified by the option ``--density`` (or ``--size `` or ``--error``).

Note: by default, the initial mesh density is 100%. User can specify the
initial density (or size) through the option:
``--initial-density`` (or ``--initial-size``, or
 ``--relative-initial-density``).



## Details of Each Command Line Option

**NOTE**: if an option has a default value, the default value is
shown in the round bracket after the option.
e.g., ``--debug-level $arg (=0)`` indicates ``$arg`` has a default value ``0``.

Options:

  ``-h [ --help ] ``<br />
  Print a help message listing all free parameters and valid arguments for
  each parameter.

  ``--debug-level $arg (=0)``<br />
  Set the debug level of the program to ``$arg``. Valid values for ``$arg`` are:

  * ``0``: no debugging message are generated.
  * ``1``: the basic mesh generation process log will be generated.
  Information regarding the completion of each major step will be displayed.
  * ``2``: more detailed information of the mesh generation will be
  displayed. e.g., each point that is removed from the mesh while
  performing mesh simplification.


  ``--history-file $arg ``<br />
  Set the history file to ``$arg``. This parameter is optional.
  If specified, the history message will be output to ``$arg``.

  ``--history-level $arg (=0) ``<br />
  Set the history level to ``$arg``. If ``$arg=0``, no history message will
  be generated; if ``$arg`` is greater than zero, with the increase of
  ``$arg``, more detailed information will be generated.

  ``--density $arg``<br />
  Set the desired mesh sampling density as a percentage of the total number
  of samples in the original image to ``$arg``. Valid values for ``$arg``
  are as follows: ``0 < $arg <= 100``.

  ``--size $arg``<br />
  Set the desired output mesh size (in number of samples) to ``$arg``.
  If specified, the value is required to be greater than 4 (at least 4
  extreme convex hull points).

  ``--error $arg``  <br />
  Set the desired mesh error threshold to ``$arg``.

  ``--relative-initial-density $arg``<br />
  Set the initial density ``d`` relative to the target density ``D``
  to ``d = $arg / (D * 100)``.

  ``--initial-density $arg (=100)``<br />
  Set the initial mesh density (in percent) to ``$arg``.

  ``--initial-size $arg``<br />
  Set the initial mesh size (in number of samples) to ``$arg``.

  ``--mesh-generation-policy $arg (=mcmg)``<br />
  Set the method used for generating the mesh to ``$arg``.
  Note: This option can take influence only when the input image is a
  color image.
  Valid values for ``$arg`` are:
  * ``mcmg`` : multi component mesh generation. This will process multiple
  components together by using the color image mesh generator.
  * ``lsmg`` : luma space mesh generation. This will convert the color
  image into luma space first, then call the grayscale image mesh generator
  to get a mesh, finally replace all sample point values with color values
  in the mesh.


  ``--initial-generator $arg (=ed) ``<br />
  Set the method used for seleting the initial mesh points to ``$arg``.
  Valid values for ``$arg`` are:
  * ``all``: set all sample points on the sampling grid as the initial
  mesh points.
  * ``random``: randomly select points on the sampling grid as the initial
  mesh points.
  * ``uniform``: uniformly select points on the sampling grid as the
  initial mesh points.
  * ``ed ``: use the error diffusion method to select the initial mesh points.


  ``--ed-size-tolerance $arg (=10) ``<br />
  Set the size tolerance of the mesh that generated by the error
  diffusion method to ``$arg``.

  ``--smooth-order $arg (=3)  ``<br />
  Set the smooth order of the smoothing filter fo ``$arg``.
  The value of ``$arg`` must be an odd positive integer.
  If one, no smoothing will be applied.
  NOTE: This option can take influence only when ``initial-generator``
  is ``ed``. <br>


  ``--smooth-operator $arg (=binomial) ``<br />
  Set the smoothing operator to ``$arg``. Note: this option can take influence
  only when ``initial-generator`` is ``ed``.
  Valid values for ``$arg`` are:
  * ``binomial``: use a binomial filter as the smoothing filter.
  * ``mean``: use a mean filter as the smoothing filter.


  ``--smooth-direction $arg (=both) ``<br />
  Set the smoothing direction to ``$arg``. Note: this option can take
  influence only when ``initial-generator`` is ``ed``.
  Valid values for ``$arg`` are:
  * ``orthogonal``: smooth on the orthogonal direction.
  * ``both``: smooth on both horizontal and vertical directions.


  ``--smooth-bound-policy $arg (=zero_ext)``<br />
   Set the boundary handling policy used for smoothing to ``$arg``.
   NOTE: This option can take influence only when ``initial-generator``
   is ``ed``.
   Valid values for ``$arg`` are:
  * ``zero_ext``: zero extension.
  * ``const_ext``: constant extension.
  * ``sym_ext``: symmetric extension.


  ``--derivative-bound-policy $arg (=zero_ext)``<br />
  Set the boundary handling policy used for partial derivative
  computation to ``$arg``. NOTE: This option can take influence only when
  ``initial-generator`` is ``ed``.
  Valid values for ``$args`` are:
  * ``zero_ext``: zero extension.
  * ``const_ext``: constant extension.
  * ``sym_ext``: symmetric extension.


  ``--use-double-convolution $arg (=0)``<br />
  Set if to do smoothing and derivative operations together in
  one convolution or saperately in two convolutions.
  NOTE: This option can take influence only when ``initial-generator``
  is ``ed``.
  Valid values for ``$args`` are:

  * ``0``: the program will do smoothing and derivative operations in
   a single convolution. The boundary policy will correspond to
  ``derivative-bound-policy``.
  * ``1``: the program will compute the partial derivative by using two
  convolution operations. Basically, it will smooth the image first using one
  convolution and then do another convolution with derivative filters.



  ``--ed-strategy $arg (=max_comp_mmsodd_ed)``<br />
  Set the error diffusion strategy to ``$arg``.
  If the input image is a color image, each of the valid values
  corresponds to a color image error diffusion strategy shown below;
  If the input image is a grayscale image, if ``$arg`` is set as one of
  ``gray_comp_laplacian_ed``, ``mean_comp_laplacian_ed``, or
  ``max_comp_laplacian_ed``, the program will use laplacian of the image
  as the density function for the error diffusion algorithm; otherwise,
  the program will use MMSODD of the image as the density function.
  Valid values for ``$arg`` are:

  * `` gray_comp_mmsodd_ed ``: convert the color image to a grayscale
   image first, then use mmsodd of the grayscale image as the density
  function and run Floyd-steinberg error diffusion method.
  * ``gray_comp_laplacian_ed``: convert the color image to a grayscale
   image first, then use laplacian of the grayscale image as the density
  function and run Floyd-steinberg error diffusion method.
  * ``vector_space_comb``: compute the maximum directional gradient of the
  vector space of R,G,B components as the density function and run
  Floyd-steinberg error diffusion method.
  * ``max_comp_mmsodd_ed``: compute the max mmsodd of all components as
   the density function and run Floyd-steinberg error diffusion method.
  * ``mean_comp_mmsodd_ed``: compute the mean mmsodd of all components
   as the density function and run Floyd-steinberg error diffusion method.
  * `` max_comp_laplacian_ed ``: compute the max laplacian of all
  components as the density function and run Floyd-steinberg error
  diffusion method.
  * ``mean_comp_laplacian_ed ``: compute the mean laplacian of all
  components as the density function and run Floyd-steinberg error
  diffusion method.
  * ``three_comps_union_mmsodd``: compute mmsodd of each component
  separately and take as three density functions; and run Floyd-steinberg
  error diffusion method on each of R,G,B components and get the union
  set. The result mesh will have desired number of points.
  * ``three_comps_union_laplacian``: compute laplacian of each component
  separately and take as three density functions; and run Floyd-steinberg
  error diffusion method on each of R,G,B components and get the union set.
  The result mesh will have desired number of points.


  ``--ed-scan-order $arg (=serpentine) ``<br />
  Set the error diffusion scan order to ``$arg``.
  Valid values for ``$arg`` are:
  * ``raster``: raster scan order.
  * ``serpentine``: serpentine scan order.


  ``--ed-leaky-mode $arg (=no_leaky) ``<br />
  Set the boundary leaky/no_leaky mode for error diffusion to ``$arg``.
  Valid values for ``$arg`` are:
  * ``leaky``: leaky mode.
  * ``no_leaky``: no_leaky mode.


  ``--ed-initial-error-mode $arg (=mirror)``<br />
  Set the error diffusion startup error condition to ``$arg``.
  Valid values for ``$arg`` are:
  * ``zero``: initial errors are all zero.
  * ``random ``: initial errors are randomly generated.
  * ``mirror``: initial errors are generated by performing Floyd-Steinberg error
  diffusion on a vertical flipped image of the original image first, then take
  the diffused errors to the last row as the startup errors.


  ``--error-metric $arg (=mean_comp_se) ``<br />
  Set the error metric for computing face error to ``$arg``.
  Valid values for ``$arg`` are:
  * ``mean_comp_se``: mean squared error of all components.
  * `` mean_comp_ae ``: mean absolute error of all components.
  * ``luma_se ``: luma space squared error.
  * ``luma_ae``: luma space absolute error.
  * ``ycbcr_vec_se``: YCbCr vector space squared error.
  * ``ycbcr_vec_ae ``: YCbCr vector space absolute error.
  * ``vec_se``: vector absolute error.


  ``--bad-point-removal $arg (=1)  ``<br />
  Sset if bad point removal is enabled.
  Valid values for ``$arg`` are:
  * ``0``: not enable bad point removal.
  * ``1``: enable bad point removel.


  ## Examples of Program Usage:

  1) Given a rgb or grayscale image ``input_image.pnm``, one can generate
  a mesh model with a sampling density of 2% using the MED method, by using
  command:

```
      make_model <input_image.pnm> output_model.model \
      --density 2 --initial-density 2 --initial-generator ed \
      --bad-point-removal 0 --ed-strategy-error-mode mirror \
      --ed-strategy max_comp_mmsodd_ed
```

  2) Given a rgb or grayscale image ``input_image.pnm``, one can generate
  a mesh model with a sampling density of 2% using the MGPRFS method by
  using command:

  ```
     make_model <input_image.pnm> output_model.model \
      --density 2 --relative-initial-density 400 --initial-generator ed \
      --bad-point-removal 0 --ed-initial-error-mode mirror \
      --ed-strategy max_comp_mmsodd_ed
  ```

  3) Given a grayscale image ``input_image.pnm``, one can generate
  a mesh model with a sampling density of 2% using the ED method by
  using command:

  ```
      make_model <input_image.pnm> output_model.model \
      --density 2 --initial-density 2 --initial-generator ed \
      --bad-point-removal 0 --ed-initial-error-mode zero
  ```

  4) Given a grayscale image ``input_image.pnm``, one can generate
  a mesh model with a sampling density of 2% using the GPR method by
  using command:

  ```
     make_model <input_image.pnm> output_model.model \
      --density 2 --initial-density 100 --initial-generator all \
      --bad-point-removal 0
  ```

  5) Given a grayscale image ``input_image.pnm``, one can generate
  a mesh model with a sampling density of 2% using the GPRFS-ED method by
  using command:

  ```
      make_model <input_image.pnm> output_model.model \
      --density 2 --relative-initial-density 400 --initial-generator ed \
      --bad-point-removal 0 --ed-initial-error-mode zero
  ```

# The ``rasterize_model`` Program
## Synopsis
``rasterize_model <input_file> output_file``

## Description
The ``rasterize_model`` program is used to reconstruct an image from
a mesh model. The program reads a mesh model from standard input,
reconstructs an image from the mesh model, and writes the image to standard
output. The reconstructed image is written in PGM/PNM format for grayscale
images, and in PPM/PNM format for color images.

## Examples of Program Usage:
Given a mesh model ``input.model``, one can generate the reconstructed image
by using command:
```
      rasterize_model <input.model> reconstructed_image.pnm
```
