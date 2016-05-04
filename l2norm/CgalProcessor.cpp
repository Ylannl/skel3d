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

#include "CgalProcessor.h"

CgalProcessor::CgalProcessor(const char *inputFile){
    if (dt.number_of_vertices() != 0)
        clear();
    
    std::string inputFile_bounds = inputFile;
    inputFile_bounds.append(".bounds");
    
    std::ifstream in_bounds(inputFile_bounds.c_str());
    in_bounds >> minx >> maxx >> miny >> maxy >> minz >> maxz;
    in_bounds.close();

    std::ifstream in(inputFile);
    std::string inputFileStr = inputFile;
    double x,y,z;
    bool e;

    if(inputFileStr.find("xyze")!=std::string::npos){
        std::cout << "Detected xyze!" << std::endl;
        while (in >> x >> y >> z >> e)
        {
            Vertex_handle v = dt.insert(PointDt(x,y,z));
            v->info().regionSmooth = e;
        }
    } else {
        std::istream_iterator<Point3D> begin(in), end;
        dt.insert(begin, end);
    }
    in.close();
    
    std::cout << "Succesfully opened " << inputFile << std::endl;
    std::cout << "  X range: " << minx << " to " << maxx << " ("<<maxx-minx<<" m)" << std::endl;
    std::cout << "  Y range: " << miny << " to " << maxy << " ("<<maxy-miny<<" m)" << std::endl;
    std::cout << "  Z range: " << minz << " to " << maxz << " ("<<maxz-minz<<" m)" << std::endl;
    std::cout << "  # Points: " << dt.number_of_vertices() << std::endl;
    std::cout << "  # Triangles: " << dt.number_of_faces() << std::endl;
}

CgalProcessor::CgalProcessor(float *pts, int N){
    if (dt.number_of_vertices() != 0)
        clear();
    
    double x,y,z;
    minx = maxx = pts[0];
    miny = maxy = pts[1];
    minz = maxz = pts[2];

    for(int i=0; i<N; i++){
        x = pts[i*3+0];
        y = pts[i*3+1];
        z = pts[i*3+2];
        dt.insert( Point3D(x,y,z) );
        if(x < minx) minx=x;
        if(x > maxx) maxx=x;
        if(y < miny) miny=y;
        if(y > maxy) maxy=y;
        if(z < minz) minz=z;
        if(z > maxz) maxz=z;
    }
    
    std::cout << "  X range: " << minx << " to " << maxx << " ("<<maxx-minx<<" m)" << std::endl;
    std::cout << "  Y range: " << miny << " to " << maxy << " ("<<maxy-miny<<" m)" << std::endl;
    std::cout << "  Z range: " << minz << " to " << maxz << " ("<<maxz-minz<<" m)" << std::endl;
    std::cout << "  # Points: " << dt.number_of_vertices() << std::endl;
    std::cout << "  # Triangles: " << dt.number_of_faces() << std::endl;
}

CgalProcessor::CgalProcessor(CgalProcessor& otherCp){
    minx = otherCp.minx;
    maxx = otherCp.maxx;
    miny = otherCp.miny;
    maxy = otherCp.maxy;
    minz = otherCp.minz;
    maxz = otherCp.maxz;
    
    dt = otherCp.dt;
    
}

void CgalProcessor::printTags(){
    for( Face_iterator ib = dt.faces_begin();
        ib != dt.faces_end(); ++ib) {
        std::cout << "this one is: " << ib->info().tooBig << std::endl;
    }
    std::cout << "end" << std::endl;
}

void CgalProcessor::markBigTriangles_mel(double maxEdgeLength){

    for (Dt::All_edges_iterator ei = dt.all_edges_begin(); ei != dt.all_edges_end(); ++ei) {
        // check if this face hasnt been marked as too big already
        if (ei->first->info().tooBig == false) {
            // compute edge length
            Point3D p1 = ei->first->vertex(ei->first->cw(ei->second))->point();
            Point3D p2 = ei->first->vertex(ei->first->ccw(ei->second))->point();
            double edgeLength = CGAL::squared_distance(p1, p2);
    //        std::cout << edgeLength << std::endl;
            // compare to userdefined maximum and if needed mark this face
            if (edgeLength > maxEdgeLength){
                ei->first->info().tooBig = true;
                ei->first->neighbor(ei->second)->info().tooBig = true;
            }
        }
    }
    
}

void CgalProcessor::markBigTriangles_mta(double maxDoubleArea){
    // iterate over all triangle faces
    for( Face_iterator ib = dt.finite_faces_begin();
        ib != dt.finite_faces_end(); ++ib) {

        // shorthand notations for the 3 triangle vertices and their position w.r.t. the contouring depth
        PointDt p1 = ib->vertex(0)->point();
        PointDt p2 = ib->vertex(1)->point();
        PointDt p3 = ib->vertex(2)->point();
        
        double doubleArea = p1.x()*p2.y() + p1.y()*p3.x() + p2.x()*p3.y() - p1.x()*p3.y() - p1.y()*p2.x() - p2.y()*p3.x();
        
        if (doubleArea > maxDoubleArea){
            ib->info().tooBig = true;
        }
    }
}

// Evaluates whether a vertex is position under, above or on (within range +/- e) the (contour) depth.
inline int CgalProcessor::cntrEvalVertex(Vertex_handle v, double depth) {
    double e = 1e-7;
    double z = v->point().z();

    if (z < depth - e)
        return -1;
    else if (depth + e < z)
        return 1;
    else
        return 0;
}

// Calculates and returns the point of intersection between the (contour) depth and the edge formed by two vertices v1 and v2.
inline PointDt CgalProcessor::cntrIntersectEdge(Vertex_handle v1, Vertex_handle v2, double depth) {
    PointDt p1 = v1->point();
    PointDt p2 = v2->point();
    
    double lambda = (depth - p1.z()) / (p2.z() - p1.z());
    double x = (1-lambda) * p1.x() + lambda * p2.x();
    double y = (1-lambda) * p1.y() + lambda * p2.y();

    return PointDt(x,y,depth);
}

// Extracts the complete set of contour lines from the this triangulation for a given isoDepth.
void CgalProcessor::extractContour(contourSegmentVec& segmentVec, double isoDepth) {
    
    // faceCache is used to ensure line segments are outputted only once. It will contain faces that have an edge exactly on the contouring depth.
    std::set<Face_handle> faceCache;
    
    // iterate over all triangle faces
    for( Face_iterator ib = dt.finite_faces_begin();
        ib != dt.finite_faces_end(); ++ib) {
        
        // shorthand notations for the 3 triangle vertices and their position w.r.t. the contouring depth
        Vertex_handle v0 = ib->vertex(0);
        Vertex_handle v1 = ib->vertex(1);
        Vertex_handle v2 = ib->vertex(2);
        int v0_ = cntrEvalVertex(v0, isoDepth);
        int v1_ = cntrEvalVertex(v1, isoDepth);
        int v2_ = cntrEvalVertex(v2, isoDepth);
        
        // following is a big if-else-if statement to identify the basic triangle configuration (wrt the contouring depth)
        
        //its on a horizontal plane: skip it
        if (v0_ == v1_ && v1_ == v2_)
            continue;
        
        //one edge is equal to isodepth: extract that edge. Use faceCache to check if this segment hasn't been extracted earlier.
        else if (v0_ == 0 && v1_ == 0) {
            faceCache.insert(ib);
            if( faceCache.find(ib->neighbor(2)) == faceCache.end() )
                segmentVec[isoDepth].push_back(CGAL::Segment_3<K>(v0->point(), v1->point()));
        } else if (v1_ == 0 && v2_ == 0) {
            faceCache.insert(ib);
            if( faceCache.find(ib->neighbor(0)) == faceCache.end() )
                segmentVec[isoDepth].push_back(CGAL::Segment_3<K>(v1->point(), v2->point()));            
        } else if (v2_ == 0 && v0_ == 0) {
            faceCache.insert(ib);
            if( faceCache.find(ib->neighbor(1)) == faceCache.end() )
                segmentVec[isoDepth].push_back(CGAL::Segment_3<K>(v2->point(), v0->point()));
        
        //there is an intersecting line segment in between the interiors of 2 edges: calculate intersection points and extract that edge
        } else if ( (v0_ == -1 && v1_ == 1 && v2_ == 1) or (v0_ == 1 && v1_ == -1 && v2_ == -1) ){
            PointDt p1 = cntrIntersectEdge(v0, v1, isoDepth);
            PointDt p2 = cntrIntersectEdge(v0, v2, isoDepth);
            segmentVec[isoDepth].push_back(CGAL::Segment_3<K>(p1, p2));
        } else if ( (v0_ == 1 && v1_ == -1 && v2_ == 1) or (v0_ == -1 && v1_ == 1 && v2_ == -1) ) {
            PointDt p1 = cntrIntersectEdge(v1, v0, isoDepth);
            PointDt p2 = cntrIntersectEdge(v1, v2, isoDepth);
            segmentVec[isoDepth].push_back(CGAL::Segment_3<K>(p1, p2));
        } else if ( (v0_ == 1 && v1_ == 1 && v2_ == -1) or (v0_ == -1 && v1_ == -1 && v2_ == 1) ) {
            PointDt p1 = cntrIntersectEdge(v2, v0, isoDepth);
            PointDt p2 = cntrIntersectEdge(v2, v1, isoDepth);
            segmentVec[isoDepth].push_back(CGAL::Segment_3<K>(p1, p2));
        
        // one vertex is on the isodepth the others are above and below: return segment, consisting out of the vertex on the isodepth and the intersection on the opposing edge
        } else if ( v0_ == 0 && v1_ != v2_ ) {
            PointDt p = cntrIntersectEdge(v1, v2, isoDepth);
            segmentVec[isoDepth].push_back(CGAL::Segment_3<K>(v0->point(), p));
        } else if ( v1_ == 0 && v0_ != v2_ ) {
            PointDt p = cntrIntersectEdge(v0, v2, isoDepth);
            segmentVec[isoDepth].push_back(CGAL::Segment_3<K>(v1->point(), p));
        } else if ( v2_ == 0 && v0_ != v1_ ) {
            PointDt p = cntrIntersectEdge(v0, v1, isoDepth);
            segmentVec[isoDepth].push_back(CGAL::Segment_3<K>(v2->point(), p));
        }
    }
}

// extract contour sets for a range of iso-depths
contourSegmentVec CgalProcessor::extractContours(std::vector<double> isoDepths) {
    contourSegmentVec segmentVec;
    for(std::vector<double>::iterator iD = isoDepths.begin(); iD != isoDepths.end(); ++iD) {
        extractContour(segmentVec, *iD);
    }
    return segmentVec;
}

contourSegmentVec CgalProcessor::extractContoursCgalInt(std::vector<double> isoDepths) {
    // this method might extract some line segments twice: if we are contouring exactly along a triangle edge. Possibly this gives problems later with the geos linemerger, resulting in broken contours. Especially problematic for gridded input (e.g. de Waal dataset)
    contourSegmentVec segmentVec;
    
    for(std::vector<double>::iterator isoDepth = isoDepths.begin(); isoDepth != isoDepths.end(); ++isoDepth) {
        
        CGAL::Plane_3<K> plane(0,0,-1,*isoDepth);
        
        for( Face_iterator ib = dt.faces_begin();
            ib != dt.faces_end(); ++ib) {
            if(!ib->info().tooBig){
                CGAL::Object result = CGAL::intersection(plane, dt.triangle(ib));
                //std::cout << CGAL::do_intersect(plane, dt.triangle(ib));
                if (const CGAL::Segment_3<K> *isegment = CGAL::object_cast<CGAL::Segment_3<K> >(&result)) {
                    segmentVec[*isoDepth].push_back(*isegment);
                }
            }
        }
    }
    return segmentVec;
}

//segmentVec CgalProcessor::extractVoronoi() {
//    segmentVec sVec;
//    
//    for( Face_iterator ib = dt.finite_faces_begin();
//        ib != dt.finite_faces_end(); ++ib) {
//    
//        
//    }
//    
//    return sVec;
//}


void CgalProcessor::clear(){
    dt.clear();
}

// dump all triangle vertices' coordinates to plain ascii
void CgalProcessor::dumpXYZ(const char * outFile){
    
    std::string outFile_bounds = outFile;
    outFile_bounds.append(".bounds");
    
    std::ofstream ofs_bounds(outFile_bounds.c_str());
    ofs_bounds << minx << std::endl << maxx << 
     std::endl << miny << std::endl << maxy << 
     std::endl << minz << std::endl << maxz << std::endl;
    
    ofs_bounds.close();
    
    std::ofstream ofs(outFile);
    
    
    for( Dt::Finite_vertices_iterator vit= dt.finite_vertices_begin(); vit != dt.finite_vertices_end() ; ++vit) {
        ofs <<std::setprecision(2)<<std::fixed;
        
        ofs << " " << (vit->point().cartesian(0));
        ofs << " " << (vit->point().cartesian(1));
        ofs << " " << (vit->point().cartesian(2));
        ofs << std::endl;
    }
    
    ofs.close();
}

//dump triangulation to simple wavefront .obj format
void CgalProcessor::dumpOBJ(const char * outFile, double zExageration){
   
    std::ofstream ofs(outFile);
    
    std::map<Vertex_handle,int> V;

    int jnum = 0;
    double x_offset = (maxx-minx)/2 + minx;
    double y_offset = (maxy-miny)/2 + miny;
    double z_offset = (maxz-minz)/2 + minz;
    for( Dt::Finite_vertices_iterator vit= dt.finite_vertices_begin(); vit != dt.finite_vertices_end() ; ++vit) {
        ofs << "v"<<std::setprecision(6)<<std::showpoint;

        ofs << " " << (vit->point().cartesian(0) - x_offset) ;
        ofs << " " << (vit->point().cartesian(1) - y_offset) ;
        ofs << " " << zExageration*(vit->point().cartesian(2) - z_offset);

        ofs << std::endl;
        V[vit] = ++jnum;
    }
    
    // vertices of the faces
    
    for( Face_iterator ib = dt.faces_begin();
        ib != dt.faces_end(); ++ib) {
        ofs << "f";
        for(int i = 0; i < 3 ; i++) {
            ofs << " " << V[ib->vertex(i)];
        }
        ofs << std::endl;
    }
    
    ofs.close();
}

double CgalProcessor::MaxX(){
    return maxx;
}
double CgalProcessor::MinX(){
    return minx;
}
double CgalProcessor::MaxY(){
    return maxy;
}
double CgalProcessor::MinY(){
    return miny;
}
double CgalProcessor::MaxZ(){
    return maxz;
}
double CgalProcessor::MinZ(){
    return minz;
}
Dt& CgalProcessor::t(){
    Dt &dtr = dt;
    return dtr;
}

// calculates the weight for Laplace Interpolant
inline double CgalProcessor::LaplaceWeight(double tsx, double tsy, double ttx, double tty, double vsx, double vsy, double vtx, double vty){
    return sqrt( pow( ( vtx - vsx ),2) +
                 pow( ( vty - vsy ),2)    ) /
            sqrt( pow( ( ttx - tsx ),2) +
                 pow( ( tty - tsy ),2)    );
}

// estimate the depth on (x,y) position of vertex v after removal of that vertex, using Linear interpolation.
double CgalProcessor::estimateZ_LIN(Vertex_handle v) throw(OutsideConvexHullException)
{
    PointDt p1,p2,p3,q = v->point();
    
    Dt t2;
    Dt::Vertex_circulator vc = dt.incident_vertices(v), done(vc);
    do{
        if(!dt.is_infinite(vc))
            t2.insert(vc->point());
    } while(++vc!=done);
    
    Dt::Face_handle face = t2.locate(q);
    if (face==NULL) throw OutsideConvexHullException();
        
    p1 = face->vertex(0)->point();
    p2 = face->vertex(1)->point();
    p3 = face->vertex(2)->point();
    CGAL::Plane_3<K> plane(p1,p2,p3);
    
    return - plane.a()/plane.c() * q.x() - plane.b()/plane.c()*q.y() - plane.d()/plane.c();
}

// estimate the depth on (x,y) position of vertex v after removal of that vertex, using Natural Neighbour interpolation.
double CgalProcessor::estimateZ_NN(Vertex_handle v) throw(OutsideConvexHullException)
{
    typedef std::vector< std::pair< Dt::Vertex_handle, Dt::Geom_traits::FT > > Point_coordinate_vector;
    
    Dt t2;
    Dt::Vertex_circulator vc = dt.incident_vertices(v),
    done(vc);
    do{
        if(!dt.is_infinite(vc))
            t2.insert(vc->point());
        else throw OutsideConvexHullException(); // meaning: this vertex is on the convexhull
            
    } while(++vc!=done);
    
    Point_coordinate_vector coords;
    
    CGAL::Triple<
    std::back_insert_iterator<Point_coordinate_vector>,
    K::FT, bool> result =
        CGAL::natural_neighbor_coordinates_vertex_2(t2, v->point(), std::back_inserter(coords));
    
    double newZ = 0;

    if(result.third)
        for(Point_coordinate_vector::iterator it = coords.begin(); it != coords.end(); ++it)
            newZ += ( it->first->point().z() * (it->second/result.second) );

    return newZ;
}


// estimate the depth on (x,y) position of vertex v, using Laplace Interpolant.
double CgalProcessor::estimateZ(smoothAlg alg, PointDt p) throw(OutsideConvexHullException)
{
    try {
        //temporarily insert vertex at p:
        Vertex_handle v = dt.insert(p);
        
        double newZ;
        
        // estimate value
        newZ = estimateZ(alg, v);
        
        //and remove the vertex again:
        dt.remove(v);
        
        return newZ;
        
    } catch (OutsideConvexHullException& e) {
        throw e;
    }
}

// estimate the depth on (x,y) position of vertex v after removal of that vertex, using Laplace Interpolant.
double CgalProcessor::estimateZ_LP(Vertex_handle v) throw(OutsideConvexHullException)
{
    std::vector< double > rawWeights, zDepth;
    
    Dt::Edge_circulator iEdges = dt.incident_edges(v), done(iEdges);
    do {
        
        if( (!dt.is_infinite(iEdges->first)) &&  (!dt.is_infinite(iEdges->first->neighbor(iEdges->second))) )
        {
            PointDt vs = dt.circumcenter(iEdges->first);
            PointDt vt = dt.circumcenter(iEdges->first->neighbor(iEdges->second));
            PointDt ts = iEdges->first->vertex(Dt::ccw(iEdges->second))->point();
            PointDt tt = iEdges->first->vertex(Dt::cw(iEdges->second))->point();
            
            rawWeights.push_back(LaplaceWeight(ts.x(), ts.y(), tt.x(), tt.y(), vs.x(), vs.y(), vt.x(), vt.y()));
            zDepth.push_back(ts.z());
        } else throw OutsideConvexHullException(); // meaning: this vertex is on the convexhull
        
    } while(++iEdges != done);
    
    double sum = std::accumulate(rawWeights.begin(), rawWeights.end(), 0.0);
    double newZ=0;

    for(size_t i=0; i<rawWeights.size(); ++i)
        newZ += (rawWeights[i]/sum)*zDepth[i];
    
    return newZ;
}

// estimate the depth on (x,y) position of vertex v after removal of that vertex, using specified interpolation method (smoothAlg).
double CgalProcessor::estimateZ(smoothAlg algorithm, Vertex_handle v) throw(OutsideConvexHullException)
{
    try {
        if(algorithm == NN)
            return estimateZ_NN(v);
        else if(algorithm == LP)
            return estimateZ_LP(v);
        else if(algorithm == LIN)
            return estimateZ_LIN(v);
        else {
            std::cerr << "Undefined interpolation algorithm";
            exit(1);
        }
    }
    catch (OutsideConvexHullException& e) {
        throw e;
    }
}

void CgalProcessor::smoothRegion(smoothAlg algorithm, bool upOnly)
{
    
    std::vector< Vertex_handle > tmpVertexVec;
    typedef std::pair< PointDt, double > pdPair;
    std::vector< pdPair > regionPoints;
    
    for( Dt::Finite_vertices_iterator vit=dt.finite_vertices_begin() ; vit != dt.finite_vertices_end(); ++vit ) {
        
        // temporarily remove all marked vertices and push them in a vector
        if (vit->info().regionSmooth) {
            regionPoints.push_back(std::make_pair(vit->point(), vit->point().z()));
            tmpVertexVec.push_back(vit);
        }
        
    }
    
    for( std::vector<Vertex_handle>::iterator it = tmpVertexVec.begin(); it != tmpVertexVec.end(); ++it ) {
        dt.remove(*it);
    }
    
    
    // estimate new depths
    for( std::vector<pdPair>::iterator it = regionPoints.begin(); it != regionPoints.end(); ++it ) {
        //dt.remove(it->first);
        Vertex_handle v = dt.insert(it->first);
        
        try {
            double newZ = estimateZ(algorithm, v);
            
            if (it->second < newZ)
                it->second = newZ;
            
            //            std::cerr << oldZ << " | " << newZ << std::endl;
        } catch (OutsideConvexHullException& e) {
            //            std::cerr << "ohoh - outside convex hull";
        }
        
        dt.remove(v);
    }
    
    // re-insert regionpoints to triangulation
    for( std::vector<pdPair>::iterator it = regionPoints.begin(); it != regionPoints.end(); ++it ) {
        dt.insert(PointDt(it->first.x(), it->first.y(), it->second));
    }
    
    std::cout << "Smoothing success" << std::endl;
}

void CgalProcessor::smooth(smoothAlg algorithm, bool upOnly)
{
    typedef std::pair< Vertex_handle, double > vdPair;
    std::vector< vdPair > newDepths;
    
    for( Dt::Finite_vertices_iterator vit=dt.finite_vertices_begin() ; vit != dt.finite_vertices_end(); ++vit ) {
        
        double oldZ = vit->point().z();
        
        try {
            double newZ = estimateZ(algorithm, vit);
            
            if (!upOnly) {
                newDepths.push_back(std::make_pair(vit, newZ));
            } else if (oldZ < newZ) {
                newDepths.push_back(std::make_pair(vit, newZ));
            }
            
//            std::cerr << oldZ << " | " << newZ << std::endl;
        } catch (OutsideConvexHullException& e) {
//            std::cerr << "ohoh - outside convex hull";
        }
    
    }
    
    for( std::vector<vdPair>::iterator it = newDepths.begin(); it != newDepths.end(); ++it ) {
        //dt.remove(it->first);
        it->first->set_point(PointDt(it->first->point().x(), it->first->point().y(),it->second ));
    }
    std::cout << "Smoothing success" << std::endl;
}

//void CgalProcessor::smoothDrop(smoothAlg algorithm, double treshold)
//{
//    typedef std::pair< Vertex_handle, double > vdPair;
//    std::vector< vdPair > newDepths;
//    
//    for( Dt::Finite_vertices_iterator vit=dt.finite_vertices_begin() ; vit != dt.finite_vertices_end(); ++vit ) {
//        
//        double oldZ = vit->point().z();
//        
//        double newZ = estimateZ(algorithm, vit);
//        
//        //if (oldZ < newZ){
//        if (newZ > oldZ)
//            if ((newZ-oldZ) < treshold)
//                newDepths.push_back(std::make_pair(vit, newZ));
//    }
//    
//    for( std::vector<vdPair>::iterator it = newDepths.begin(); it != newDepths.end(); ++it ) {
//        if(dropPoints)
//            dt.remove(it->first);
//        else
//            it->first->set_point(PointDt(it->first->point().x(), it->first->point().y(),it->second ));
//    }
//}

void CgalProcessor::densify(smoothAlg algorithm, bool tooBigOnly)
{
    typedef std::vector < std::pair< PointDt, Dt::Face_handle > > pfPair;
    pfPair newPoints;
    for( Dt::Finite_faces_iterator ib = dt.finite_faces_begin();
        ib != dt.finite_faces_end(); ++ib) {
        if ((ib->info().tooBig && tooBigOnly) || (!tooBigOnly) ){
            
            // Wikipedia says this is how to calculate the coordinates of the incenter
            // http://en.wikipedia.org/wiki/Incenter#Coordinates_of_the_incenter
//            PointDt pa = ib->vertex(0)->point();
//            PointDt pb = ib->vertex(1)->point();
//            PointDt pc = ib->vertex(2)->point();
//            
//            double a = CGAL::squared_distance(pb, pc);
//            double b = CGAL::squared_distance(pa, pc);
//            double c = CGAL::squared_distance(pa, pb);
//            
//            double P = a+b+c;
//            
//            double ic_x = (a*pa.x() + b*pb.x() + c*pc.x())/P;
//            double ic_y = (a*pa.y() + b*pb.y() + c*pc.y())/P;
//            
//            newPoints.push_back(std::make_pair(PointDt(ic_x,ic_y,0), ib) );

            // but circumcenter should be better than incircle
            newPoints.push_back(std::make_pair(dt.circumcenter(ib), ib) );
            
        }
    }
                            
    for (pfPair::iterator p = newPoints.begin(); p != newPoints.end(); ++p){        
        try{
            Vertex_handle v = dt.insert(p->first, p->second);
            
            double newZ = estimateZ(algorithm, v);
            
            v->set_point(PointDt(v->point().x(), v->point().y(), newZ ));
        } catch (OutsideConvexHullException& e) {
            // do something?
        }
    }
    std::cout << "Densification success" << std::endl;
}

double CgalProcessor::vertError(smoothAlg algorithm, Vertex_handle v) throw(OutsideConvexHullException)
{
    try {
        return v->point().z() - estimateZ(algorithm, v);
    } catch (OutsideConvexHullException& e) {
        throw e;
    }
    
}
    
double CgalProcessor::verticalError(PointDt q)
{
    PointDt p1,p2,p3;
    
    // what if no triangle "contains" this point?
    Dt::Face_handle face = dt.locate(q);
    
    p1 = face->vertex(0)->point();
    p2 = face->vertex(1)->point();
    p3 = face->vertex(2)->point();
    CGAL::Plane_3<K> plane(p1,p2,p3);
    
    double interpolatedZ = - plane.a()/plane.c() * q.x() - plane.b()/plane.c()*q.y() - plane.d()/plane.c();
    
    return q.z()-interpolatedZ; //std::fabs()
}

void CgalProcessor::updateCache(SimplificationCache &cache, smoothAlg algorithm, Vertex_handle v, bool upOnly)
{
    try {
        // following function call might throw OutsideConvexHullException
        double e = vertError(algorithm, v);
        
        // Are we performing 'safe' (only going up in this particular iteration) or regular simplification?
        if (upOnly){
            // negative e means it is 'safe' to drop this points. Attempt to account for floating point error
            if (e < -1E-5){
                cache.verts[v] = std::fabs(e);
                
            // remove point from simplification cache if its about to be moved downwards
            } else if (cache.verts.find(v) == cache.verts.end()) {
                cache.verts.erase(v);
            }
        } else {
            cache.verts[v] = std::fabs(e);
        }
        
        // keep record of the minimal error
        if (cache.minError.second > e){
            cache.minError = std::make_pair(v,e);
        }
    // skip nasty points that might cause code to crash
    } catch (OutsideConvexHullException& e) {
        std::cerr << "Simpflification: Skipping point on boundary";
    }
    
}

void CgalProcessor::simplify(smoothAlg algorithm, double treshold, bool upOnly, bool dumpIntermediateStepsToOBJ)
{
    SimplificationCache cache;
    
    Vertex_handle vh = dt.finite_vertices_begin();
    cache.minError = std::make_pair(vh,vertError(algorithm, vh));
    
    for( Dt::Finite_vertices_iterator vit=dt.finite_vertices_begin() ; vit != dt.finite_vertices_end(); ++vit ) {
        updateCache(cache, algorithm, vit, upOnly);
    }
    
    int numberOfDroppedPoints = 0;
    while(cache.minError.second < treshold) {
        vdMap::iterator currentPosIt = cache.verts.begin();
        cache.minError = *currentPosIt;
        
        for (vdMap::iterator it = cache.verts.begin(); it != cache.verts.end(); ++it) {
            if (cache.minError.second > it->second) {
                cache.minError = *it;
                currentPosIt = it;
            }
        }
        std::cout << cache.minError.first->point() << " | " << cache.minError.second << "\n";        

        Dt::Vertex_circulator iVertex = dt.incident_vertices(cache.minError.first), done(iVertex);
                
        std::vector<Vertex_handle> neighbours;
        do{
            if(!dt.is_infinite(iVertex)) 
                neighbours.push_back(iVertex);
        } while(++iVertex != done);
        
        dt.remove(cache.minError.first);
        cache.verts.erase(currentPosIt);
        numberOfDroppedPoints++;
        
        for(std::vector<Vertex_handle>::iterator iVertex = neighbours.begin(); iVertex != neighbours.end(); ++iVertex){
            updateCache(cache, algorithm, *iVertex, upOnly);
        }
        if(dumpIntermediateStepsToOBJ) {
            std::stringstream outPath;
            outPath << "/Volumes/macData/mScpp/testdata/intermediate_output/dump_" << numberOfDroppedPoints << ".obj";
            dumpOBJ(outPath.str().c_str());
        }
    }
    std::cout << "# points dropped for simplify: " << numberOfDroppedPoints << std::endl;
    
}

void CgalProcessor::metricL2(CgalProcessor &otherSurface)
{
    Dt & otherDt = otherSurface.t();
    AABB_Tree TriTree;

    std::list<Triangle> triangleList;
    for( Dt::Finite_faces_iterator ib = otherDt.finite_faces_begin();
        ib != otherDt.finite_faces_end(); ++ib) {
        triangleList.push_back(otherDt.triangle(ib));
    }
    
    TriTree.insert(triangleList.begin(), triangleList.end());
    if (TriTree.accelerate_distance_queries()) std::cout << "build AABB_tree succes\n"; else std::cout << "build AABB_tree fail\n";
    
    double min=1, max=0;
    for( Dt::Finite_vertices_iterator vit=dt.finite_vertices_begin() ; vit != dt.finite_vertices_end(); ++vit ) {
        double d = TriTree.squared_distance(vit->point());
        vit->info().metricL2 = d;
        if (min>d) min = d;
        if (max<d) max = d;
        std::cout << "point-distance: " << d  << std::endl;
    }
    maxL2error = max;
    minL2error = min;
    std::cout << "metricL2: max="<<max<<", min="<<min<<std::endl;
}

void CgalProcessor::metricL2potri(float* pts, int N, float* out)
{
    // build index for this triangulation
    AABB_Tree TriTree;

    std::list<Triangle> triangleList;
    for( Dt::Finite_faces_iterator ib = dt.finite_faces_begin();
        ib != dt.finite_faces_end(); ++ib) {
        triangleList.push_back(dt.triangle(ib));
    }
    
    TriTree.insert(triangleList.begin(), triangleList.end());
    if (TriTree.accelerate_distance_queries()) std::cout << "build AABB_tree succes\n"; else std::cout << "build AABB_tree fail\n";
    
    // double min=1, max=0;
    // for( Dt::Finite_vertices_iterator vit=dt.finite_vertices_begin() ; vit != dt.finite_vertices_end(); ++vit ) {
    for( int i=0; i<N; i++ ) {
        Point3D p = Point3D(pts[3*i+0], pts[3*i+1], pts[3*i+2]);
        Point3D q = TriTree.closest_point(p);
        // double d = sqrt(TriTree.squared_distance(p));
        out[3*i+0] = q.x();
        out[3*i+1] = q.y();
        out[3*i+2] = q.z();
        // if (min>d) min = d;
        // if (max<d) max = d;
    }
    // maxL2error = max;
    // minL2error = min;
    // std::cout << "metricL2: max="<<max<<", min="<<min<<std::endl;
}

void CgalProcessor::metricSafety(CgalProcessor &otherSurface)
{
    double NofUnsafePoints=0;
    for( Dt::Finite_vertices_iterator vit=dt.finite_vertices_begin() ; vit != dt.finite_vertices_end(); ++vit ) {
        double d = otherSurface.verticalError(vit->point());
        
        bool notSafe = d>1E-5;
        vit->info().metricSafety = notSafe;
        if (notSafe) {
            std::cout << "Vertical distance: " << d << " : is marked unsafe!" << std::endl;

            NofUnsafePoints++;
        } 

    }
    std::cout << "# unsafe points: " << NofUnsafePoints << std::endl;
}

void CgalProcessor::dumpDiffs(CgalProcessor &otherSurface, const char * outFile)
{
    
    std::ofstream ofs(outFile);
    ofs <<std::setprecision(2)<<std::fixed;
    
    for( Dt::Finite_vertices_iterator vit=dt.finite_vertices_begin() ; vit != dt.finite_vertices_end(); ++vit ) {
        double d = otherSurface.verticalError(vit->point());
        
        if (d<1E-5 and d>-1E-5)
            d=0;

        ofs << " " << (vit->point().cartesian(0));
        ofs << " " << (vit->point().cartesian(1));
        ofs << " " << (vit->point().cartesian(2));
        ofs << " " << d;
        ofs << std::endl;
        
    }
    
    ofs.close();
}