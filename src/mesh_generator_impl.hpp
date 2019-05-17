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

inline bool
Vertex_priority_queue_compare::operator()(const Mesh_vertex* v1,
                                          const Mesh_vertex* v2) const
{
  return (v1->get_sig() != v2->get_sig()) ? 
    (v1->get_sig() > v2->get_sig()) : (v1->point() > v2->point());
}

inline Mesh_vertex::Mesh_vertex()
  : on_pri_que_(false)
  , sig_(0)
  , pri_que_pos_(Vertex_priority_queue::handle_type())
{
}

inline double
Mesh_vertex::get_sig() const
{
  return sig_;
}

inline void
Mesh_vertex::set_sig(double sig)
{
  sig_ = sig;
}

inline bool
Mesh_vertex::on_pri_que() const
{
  return on_pri_que_;
}

inline void
Mesh_vertex::set_pri_que_pos(Vertex_priority_queue::handle_type pos)
{
  assert(on_pri_que_ == false);
  assert(pos != Vertex_priority_queue::handle_type());
  pri_que_pos_ = pos;
  on_pri_que_ = true;
}

inline Mesh_vertex::Vertex_priority_queue::handle_type
Mesh_vertex::get_pri_que_pos() const
{
  assert(on_pri_que_ == true);
  return pri_que_pos_;
}

inline void
Mesh_vertex::clear_pri_que_pos()
{
  pri_que_pos_ = Vertex_priority_queue::handle_type();
  on_pri_que_ = false;
}

inline Mesh_face::Mesh_face()
  : err_valid_(false)
  , err_(0)
{
}

inline double
Mesh_face::get_err() const
{
  assert(err_valid_ == true);
  return err_;
}

inline void
Mesh_face::set_err(double err)
{
  err_ = err;
  err_valid_ = true;
}

inline bool
Mesh_face::is_err_valid() const
{
  return err_valid_;
}

inline void
Mesh_face::invalidate_err()
{
  err_ = 0;
  err_valid_ = false;
}
