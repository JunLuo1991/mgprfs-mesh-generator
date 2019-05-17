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

#ifndef debug_hpp
#define debug_hpp

////////////////////////////////////////////////////////////////////////////////
// Debugging Support
////////////////////////////////////////////////////////////////////////////////

#include <cassert>

/*!
DEBUG_LEVEL determines the debug level. It should be an integer that is
greater than or equal to zero.
If zero, no debugging message are generated.
If one, the basic mesh generation process log will be generated.
  Information regarding each running step will be displayed.
  (e.g. process options, make initial mesh, initial points selection, etc.)
If two, some more detailed info of the mesh generation will be displayed.
  (e.g. each point that is been removed from the mesh while doing the mesh
   simplification process, each affected vertex that influenced by removing
   a vertex from the mesh, etc.)

*/

#if !defined(DEBUG_LEVEL)
#define DEBUG_LEVEL 0     /* not debugging */
//#define DEBUG_LEVEL 1   /* debugging, show basic mesh generation process log*/
//#define DEBUG_LEVEL 2   /* debugging, show more detaied info of mesh generation*/
#endif


/*!
expensive_assert is used to perform assertion checks that are
computationally expensive
*/
#if (DEBUG_LEVEL > 0)
// Enable expensive assertion checks (assuming that assertions are enabled).
#define expensive_assert(x) assert(x)
#else
// Disable expensive assertion checks.
#define expensive_assert(x)
#endif

#endif
