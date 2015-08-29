#ifndef dgp_common_h
#define dgp_common_h

// global header file, ���м��뵽���common project��header file����Ҫ�ڴ�register.
//
//DGPCommon.h -- 
//			   -| OpenMeshAll.h
//
//
//

// Ŀ����ֻ��OpenMesh // mesh, math, system util.
//           OpenGL/Glut/glew��ء�

// ��openmesh��ص�.
#include "OpenMeshAll.h"

// ��opengl��ص�. //���gl,glut��ص�,Ҫ�Ǹ����������ͻ�Ļ�,���ܲ��ܰ���������.
#include "gl.hh"
#include "gl.h"
#include "TrackballManipulator.h"
//#include "GlutWindow.h"  // renc removed file at 20150829.

#include "geom3d.h"

// ��Mesh��ص�
//#include "GeometryObject.h" // renc, removed file at 20150829.
#include "OpenMeshMeshType.h"

#include "barycentric.h"
#include "LaplacianT.h"
#include "GaussianCurvatureQuadT.h"

// ��simplification��ص�.
#include "QuadricT.h"
#include "VHierarchy.h"
#include "QEMDecimationT.h"

// ��Subdivision��ص�.
#include "MidpointLinearSubdT.h"
#include "CLoop.h"
#include "VertexFeatureType.h"
#include "Hoppe94LoopSubT.h"
#include "CatmullClarkSubdT.h"
#include "CompositeSqrt2SubdT.h"
#include "UnifiedSubd24dbsT.h"

// ��Subdivision evaluation��ص�.
#include "Stam98ClassicLoopEval.h"
#include "Hoppe94LoopEval.h"
#include "LoopEvaluatorT.h"


#endif // dgp_common_h

