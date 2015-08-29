#ifndef subdi_classic_loop_h
#define subdi_classic_loop_h

// 这是对OpenMesh中(classic, no sharp features)Loop细分模板的使用，以及扩充。
#include "OpenMeshMeshType.h"

namespace DGP {

//void cloop_subdi_1(TriMesh &_m);// 对参数_m做一次classic loop细分
void cloop_subdi(TriMesh &_m, int _n);//模型_m做_n次classic loop细分.

// 求参数_m的极限坐标
void cloop_limit_surf(TriMesh &_m); 

// 求mesh中的某一个顶点的极限坐标
OpenMesh::Vec3d cloop_limit_pos(const TriMesh &_m, const TriMesh::VertexHandle &_vh);

// in mesh _m, calculate the tangent vector along edge _hh.
OpenMesh::Vec3d cloop_tangent_vector_along_edge(const TriMesh &_m, const TriMesh::HalfedgeHandle & _hh);

// 求mesh中的某一个顶点的极限法向
OpenMesh::Vec3d cloop_limit_nor(const TriMesh &_m, const TriMesh::VertexHandle &_vh);
} // end of namespace DGP
#endif