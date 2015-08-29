#ifndef dgp_common_h
#define dgp_common_h

// global header file, 所有加入到这个common project的header file都需要在此register.
//
//DGPCommon.h -- 
//			   -| OpenMeshAll.h
//
//
//

// 目标是只跟OpenMesh // mesh, math, system util.
//           OpenGL/Glut/glew相关。

// 跟openmesh相关的.
#include "OpenMeshAll.h"

// 跟opengl相关的. //这跟gl,glut相关的,要是跟其它库想冲突的话,可能不能包含在这里.
#include "gl.hh"
#include "gl.h"
#include "TrackballManipulator.h"
//#include "GlutWindow.h"  // renc removed file at 20150829.

#include "geom3d.h"

// 跟Mesh相关的
//#include "GeometryObject.h" // renc, removed file at 20150829.
#include "OpenMeshMeshType.h"

#include "barycentric.h"
#include "LaplacianT.h"
#include "GaussianCurvatureQuadT.h"

// 跟simplification相关的.
#include "QuadricT.h"
#include "VHierarchy.h"
#include "QEMDecimationT.h"

// 跟Subdivision相关的.
#include "MidpointLinearSubdT.h"
#include "CLoop.h"
#include "VertexFeatureType.h"
#include "Hoppe94LoopSubT.h"
#include "CatmullClarkSubdT.h"
#include "CompositeSqrt2SubdT.h"
#include "UnifiedSubd24dbsT.h"

// 跟Subdivision evaluation相关的.
#include "Stam98ClassicLoopEval.h"
#include "Hoppe94LoopEval.h"
#include "LoopEvaluatorT.h"


#endif // dgp_common_h

