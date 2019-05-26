# An Effective Triangle Mesh Generation Software For Grayscale or RGB Color Images

This software was developed by Jun Luo as part of his Master of Applied Science
thesis project. 

This software can be used to effectively generate a Delaunay triangle mesh model
of a grayscale or RGB color image. The source code can be found in the directory src.
The detailed documentation of the software can be found in the directory doc.

A thesis wrote based on this software can be found in
[Jun Luo MASc Thesis](doc/Jun_Luo_MASc_Thesis.pdf).
For people who have no background related to this area of research, reviewing the 
thesis can be helpful. Some of the example results produced by this software can be
viewed in the thesis.

Jun Luo can be reached at the following email:
* junluo152@gmail.com

## Installation
### Prerequisite Libraries.

To compile, build, and run this program, the compiler should support C++17.
We recommend the use of GCC with a version 7.2.0 or higher as the C++ compiler.
In addition, the software implementation makes heavy use of some C++ libraries.
Ensure that all of the following libraries are installed before building the
program:
```
    1. Computational Geometry Algorithm Library (CGAL)
    http://www.cgal.org

    2. Boost Library
    https://www.boost.org

    3. OpenGL Utility Toolkit (GLUT)
    http://www.opengl.org/resources/libraries/glut/
    http://freeglut.sourceforge.net

    4. Signal Processing Library (SPL)
    http://www.ece.uvic.ca/~mdadams/SPL
   
    5. Signal Processing Library Extension Library (SPLEL)
    NOTE: THE SPLEL LIBRARY IS CURRENTLY A PRIVATE LIBRARY OWNED
    BY Michael D.Adams. THIS LIBRARY IS NOT INTENDED TO BE RELEASED
    TO PUBLIC AT THIS TIME.
```

### Build the Program
 To build the program, run following commands in order:
```
    1. mkdir -p $BINARY_DIR
    2. cmake -H. -B$BINARY_DIR
    3. cmake --build $BINARY_DIR
```


## Using the Software
### The `make_model` program.

In order to run the `make_model` program, the following options must be provided:
   * --density  (or --size  or --error): specify the target density
   * the input image file in PNM format (read from standard input)

All other options have default values and can be found by the *-h* option.
The result mesh model is output to standard output stream.

Examples of using the software:

1. use the MED method to generate the mesh for a grayscale or rgb image:
```
   $BINARY_DIR/src/make_model <input_image.pnm> output_model.model \
    --density 2 --initial-density 2 --initial-generator ed \
    --bad-point-removal 1 --ed-strategy-error-mode mirror \
    --ed-strategy max_comp_mmsodd_ed
```

2. use the MGPRFS method to generate the mesh for a grayscale or rgb image:
```
   $BINARY_DIR/src/make_model <input_image.pnm> output_model.model \
    --density 2 --relative-initial-density 400 --initial-generator ed \
    --bad-point-removal 0 --ed-initial-error-mode mirror \
    --ed-strategy max_comp_mmsodd_ed
```

More detailed information can be seen in 
[Software Command Line Interface](doc/command_line_interface.txt).

### The `rasterize_model` program

In order to run the `rasterize_model` program, you can run the following command:
```
   $BINARY_DIR/src/rasterize_model <input_model.model> recon_image.pnm
```