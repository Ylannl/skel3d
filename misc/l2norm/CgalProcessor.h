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

#ifndef PhDVis_CgalProcessor_h
#define PhDVis_CgalProcessor_h

//#define REAL double
//#define ANSI_DECLARATORS
//#define VOID void

#include "CgalType.h"


#include <map>
#include <vector>
#include <set>
#include <numeric>
#include <fstream>
#include <sstream>
#include <string>

enum smoothAlg {LP, NN, LIN};
typedef std::vector< CGAL::Segment_3<K> > segmentVec;
typedef std::map<double, segmentVec > contourSegmentVec;

// Simplfication
typedef std::map<Vertex_handle, double> vdMap; // a vertex and an error
typedef std::pair<Vertex_handle, double> vdPair;

struct SimplificationCache {
    vdMap verts;
    vdPair minError;
};

class OutsideConvexHullException: public std::exception
{
public:
    virtual const char* what() const throw()
    {
        return "Vertex lays outside triangulation convex hull";
    }
};

class CgalProcessor {
private:
    Dt dt;
    double minx, maxx, miny, maxy, minz, maxz;
    double minL2error, maxL2error;
//    const char* inputFile;
    
    inline double LaplaceWeight(double tsx, double tsy, double ttx, double tty, double vsx, double vsy, double vtx, double vty);

    //contouring
    inline int cntrEvalVertex(Vertex_handle v, double depth);
    inline PointDt cntrIntersectEdge(Vertex_handle v1, Vertex_handle v2, double depth);
    void extractContour(contourSegmentVec& segmentVec, double isoDepth);
    
    double estimateZ(smoothAlg, PointDt) throw(OutsideConvexHullException);
    double estimateZ_LP(Vertex_handle) throw(OutsideConvexHullException);
    double estimateZ_NN(Vertex_handle) throw(OutsideConvexHullException);
    double estimateZ_LIN(Vertex_handle) throw(OutsideConvexHullException);
    double estimateZ(smoothAlg, Vertex_handle) throw(OutsideConvexHullException);
    double vertError(smoothAlg, Vertex_handle) throw(OutsideConvexHullException);
    double verticalError(PointDt q);
    double verticalError(Vertex_handle);
    void updateCache(SimplificationCache &cache, smoothAlg algorithm, Vertex_handle v, bool upOnly);
    
public:
    CgalProcessor(const char *inputFile);
    CgalProcessor(CgalProcessor &);
    CgalProcessor(float *pts, int N);
    // Load data
//    void load(const char* inputFile);
    void clear();

    // export data
    void dumpOBJ(const char * outfile, double zExageration=1);
    void dumpXYZ(const char * outFile);
    void dumpDiffs(CgalProcessor &otherSurface, const char * outFile);
    contourSegmentVec extractContours(std::vector<double> isoDepths);
    contourSegmentVec extractContoursCgalInt(std::vector<double> isoDepths);

    double MaxX();
    double MinX();
    double MaxY();
    double MinY();
    double MaxZ();
    double MinZ();
    Dt& t();
    
    void markBigTriangles_mel(double maxEdgeLength);
    void markBigTriangles_mta(double maxDoubleArea);
    void printTags();

    void smooth(smoothAlg, bool upOnly=true);
    void smoothRegion(smoothAlg, bool upOnly=true);
    void densify(smoothAlg, bool tooBigOnly=false);
    void simplify(smoothAlg, double treshold, bool upOnly=true, bool dumpIntermediateStepsToOBJ=false);
        
    void metricL2(CgalProcessor &otherSurface);
    void metricL2potri(float* pts, int N, float* out);
    void metricSafety(CgalProcessor &processedSurface);
    
    friend class visVoronoi;
    friend class visContour;
    friend class visTriangulation;
    
//    const Dt& getDt();
//    Vd getVd();
};

#endif
