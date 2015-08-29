#ifndef subdi_classic_loop_h
#define subdi_classic_loop_h

// ���Ƕ�OpenMesh��(classic, no sharp features)Loopϸ��ģ���ʹ�ã��Լ����䡣
#include "OpenMeshMeshType.h"

namespace DGP {

//void cloop_subdi_1(TriMesh &_m);// �Բ���_m��һ��classic loopϸ��
void cloop_subdi(TriMesh &_m, int _n);//ģ��_m��_n��classic loopϸ��.

// �����_m�ļ�������
void cloop_limit_surf(TriMesh &_m); 

// ��mesh�е�ĳһ������ļ�������
OpenMesh::Vec3d cloop_limit_pos(const TriMesh &_m, const TriMesh::VertexHandle &_vh);

// in mesh _m, calculate the tangent vector along edge _hh.
OpenMesh::Vec3d cloop_tangent_vector_along_edge(const TriMesh &_m, const TriMesh::HalfedgeHandle & _hh);

// ��mesh�е�ĳһ������ļ��޷���
OpenMesh::Vec3d cloop_limit_nor(const TriMesh &_m, const TriMesh::VertexHandle &_vh);
} // end of namespace DGP
#endif