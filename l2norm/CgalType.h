// Copyright (c) 2012, 2013
// Ravi Peters -- r.y.peters@tudelft.nl
// All rights reserved
// 
// This file is part of Surfonoi.
// 
// Surfonoi is free software: you can redistribute it and/or modify
// it under the terms of the GNU General Public License as published by
// the Free Software Foundation, either version 3 of the License, or
// (at your option) any later version.
// 
// Surfonoi is distributed in the hope that it will be useful,
// but WITHOUT ANY WARRANTY; without even the implied warranty of
// MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
// GNU General Public License for more details.
// 
// You should have received a copy of the GNU General Public License
// along with Surfonoi.  If not, see <http://www.gnu.org/licenses/>.

#ifndef PhDVis_CgalType_h
#define PhDVis_CgalType_h

/* stupid OS X macro names */
#ifdef check
#undef check
#endif
#ifdef nil
#undef nil
#endif
#ifdef Nil
#undef Nil
#endif
#ifdef handle
#undef handle
#endif
#ifdef Handle
#undef Handle
#endif

#include <CGAL/Plane_3.h>
#include <CGAL/intersections.h>

#include <CGAL/Exact_predicates_inexact_constructions_kernel.h>
#include <CGAL/Triangulation_vertex_base_with_info_2.h>
#include <CGAL/Triangulation_face_base_with_info_2.h>
#include <CGAL/Projection_traits_xy_3.h>
#include <CGAL/Delaunay_triangulation_2.h>

#include <CGAL/Triangulation_hierarchy_2.h>

//#include <CGAL/Voronoi_diagram_2.h>
//#include <CGAL/Delaunay_triangulation_adaptation_traits_2.h>
//#include <CGAL/Delaunay_triangulation_adaptation_policies_2.h>

#include <CGAL/natural_neighbor_coordinates_2.h>

//#include <CGAL/Simple_cartesian.h>
#include <CGAL/AABB_tree.h>
#include <CGAL/AABB_traits.h>
#include <CGAL/AABB_triangle_primitive.h>
#include <list>

struct vertexInfo {
    double metricL2;
    bool metricSafety;
    //bool tooBig;
    bool regionSmooth;
    vertexInfo():metricL2(0),metricSafety(0),regionSmooth(0){}
};

struct faceInfo {
    bool tooBig;
    faceInfo():tooBig(0){}
};

//typedef CGAL::Simple_cartesian<double> K;
typedef CGAL::Exact_predicates_inexact_constructions_kernel K;

typedef K::Triangle_3 Triangle;

typedef std::list<Triangle>::iterator AABB_Iterator;
typedef CGAL::AABB_triangle_primitive<K,AABB_Iterator> AABB_Primitive;
typedef CGAL::AABB_traits<K, AABB_Primitive> AABB_triangle_traits;
typedef CGAL::AABB_tree<AABB_triangle_traits> AABB_Tree;

typedef CGAL::Projection_traits_xy_3<K>                     Gt;

typedef CGAL::Triangulation_vertex_base_with_info_2<vertexInfo,Gt>       Vbb;
typedef CGAL::Triangulation_hierarchy_vertex_base_2<Vbb>    Vb;
typedef CGAL::Triangulation_face_base_with_info_2<faceInfo,Gt>                 Fb;
typedef CGAL::Triangulation_data_structure_2<Vb,Fb>         Tds;

typedef CGAL::Delaunay_triangulation_2<Gt,Tds>              Dtt;
typedef CGAL::Triangulation_hierarchy_2<Dtt>                Dt;

typedef Dt::Vertex_handle                                 Vertex_handle;
typedef Dt::Face_handle                                   Face_handle;
typedef Dt::Vertex_iterator                               Vertex_iterator;
typedef Dt::Face_iterator                                 Face_iterator;

typedef Dt::Point                                           PointDt;
typedef K::Point_3                                          Point3D;
typedef K::Point_2                                          Point2D;
typedef CGAL::Segment_3<K>                                  Segment3D;

#endif
