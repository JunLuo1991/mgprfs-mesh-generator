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
http://www.opengl.org/resources/libraries/glut
http://freeglut.sourceforge.net

4. Signal Processing Library (SPL)
http://www.ece.uvic.ca/~mdadams/SPL

5. Signal Processing Library Extension Library (SPLEL)
NOTE: THE SPLEL LIBRARY IS CURRENTLY A PRIVATE LIBRARY
OWNED BY Michael D.Adams. THIS LIBRARY IS NOT INTENDED
TO BE RELEASED TO PUBLIC AT THIS TIME.
```

### Install the Software
In what follows, let $TOP_DIR denote the top-level directory of the mesh generator software distribution (i.e., the directory containing the
README.md file); let $BUILD_DIR denote a new directory to be created for
building the software; and let $INSTALL_DIR denote the directory in which
to install the software.

To install the mesh generator software with CMake, run following
  commands in order:

1. ``cmake -H$TOP_DIR -B$BINARY_DIR -DCMAKE_INSTALL_PREFIX=$INSTALL_DIR``<br \>
2. `` cmake --build $BINARY_DIR --target install``


## Using the Software
Detailed information on the use of the software can be found in [Software Command Line Interface](doc/command_line_interface.md).
