#ifndef dgp_geom3d_h
#define dgp_geom3d_h

/// Geometry Math
#include "OpenMeshAll.h" 

namespace DGP {
 
// dependent 3 vectors(not points) whether on the same plane or not
template<typename Scalar>
inline bool is_coplanar(const OpenMesh::VectorT<Scalar,3> &_v0, const OpenMesh::VectorT<Scalar,3> &_v1, const OpenMesh::VectorT<Scalar,3> &_v2) {
	Scalar project = dot(cross(_v0, _v1), _v2); //_v2��_v0��_v1�淨���ϵ�ͶӰ.
	//std::cout << project << ", ";
	if (fabs(project) < 1.0e-6) return true;
	else return false;
}
// ��������������Ҽ����, �ж�ʸ���Ƿ�Ϊ��
template<typename Scalar>
inline bool is_zero(const OpenMesh::VectorT<Scalar,3> &_v, double _epsilon = 1.0e-6) {
	if (fabs(_v[0]) < _epsilon && fabs(_v[1]) < _epsilon && fabs(_v[2]) < _epsilon) return true;
	else return false;
}

template<typename Scalar>
inline OpenMesh::VectorT<Scalar, 3>
normal_of_3points(const OpenMesh::VectorT<Scalar,3> &_p0, const OpenMesh::VectorT<Scalar,3> &_p1, 
				  const OpenMesh::VectorT<Scalar,3> &_p2) {
	OpenMesh::VectorT<Scalar,3> n0 = cross(_p2-_p1, _p0-_p1).normalize();
	OpenMesh::VectorT<Scalar,3> n1 = cross(_p0-_p2, _p1-_p2).normalize();
	OpenMesh::VectorT<Scalar,3> n2 = cross(_p1-_p0, _p2-_p0).normalize();
	//std::cout << n2 << ".\t";////�������һ���־�Ȼ��-1.#####<< n0 << "." << n1 << "." 
	//return ((n0 + n1 + n2) / (Scalar)3.0).normalize();
	return n2;//
}

// ������ʸ���ļн�, ��������߶��Ѿ�normalized��.
template<typename Scalar> // Make sure the _v0 and _v1 are normalized.
inline Scalar angle_degree_between_two_vec(const OpenMesh::VectorT<Scalar, 3> &_v0, const OpenMesh::VectorT<Scalar, 3> &_v1) {
	// Method 1.ȱ���ǵ����������н�Ϊ0ʱ���õ�-1.#IDF������ֵ
	double const gfxRADTODEG = 57.29577951308232087680;//�Ƕȵ�λ��ת�� 180/pi
	return acos(dot(_v0, _v1)) * gfxRADTODEG;//[0, 180]
	// Method 2.atan2�õ�����
	//return atan2(cross(_v0.normalize(), _v1.normalize()).norm(), dot(_v0.normalize(), _v1.normalize())) * gfxRADTODEG;
}
} // namespace 

#endif //dgp_geom3d_h