#ifndef hoppe94loopsubvft_h
#define hoppe94loopsubvft_h

namespace DGP {

	// vertex feature types ����ļ�����������, �Ǹ���hoppe94[Piecewise Smooth Surface Reconstruction]�������.
	// ���ڶ�����һ��ͷ�ļ���Ϊ�˱�ĳ���ͳһʹ��.
	// Usage: DGP::SMOOTH_VFT ������Ҫ DGP::v_feature_type::SMOOTH_VFT, DGP::v_feature_type::SMOOTH_VFT.
	enum v_feature_type { SMOOTH_VFT =0, DART_VFT = 1, CREASE_VFT = 2, CORNER_VFT = 3 }; 


	template <typename MeshT> 
	inline unsigned int identify_all_sharp_features(MeshT &_m, OpenMesh::VPropHandleT<v_feature_type> &_vp_type, 
		bool _debug = false, double _angle_tresh = 44.0) {
			// �ҳ�����������ߺ�������
			//         �ߵ�����            ��������
			// edge -> crease edge		-> dart				��ʵҲ������OpenMesh�е�feature״̬,��cornerһ��,���Ǽ�ʱ���ܶ���.
			//							-> crease,			
			//							-> corner vertex.	��ʵҲ������OpenMesh�е�feature״̬,��dartһ��,���Ǽ�ʱ���ܶ���.
			//		-> boundary edge	-> boundary vertex.
			// OpenMesh��bool is_boundary(VertexHandle/EdgeHandle/HalfedgeHandle/FaceHandle )���Բ��Ա߽�����,����������Բ��ü����

			// һ��property������valid�Ļ��ͱ�������û�б�add.
			if (_vp_type.is_valid() == false) {
				std::cout << "Error: DGP::identify_all_sharp_features: \n"
					<< "\t�жϵ�����ǰ_vp_type is not valid, make sure mesh.add_property(_vp_type) before.\n";
			}

			typename MeshT::VIter v_it, v_end(_m.vertices_end());
			double const gfxDEGTORAD =  0.01745329251994329577;//�Ƕȵ�λ��ת�� pi/180 degree to radian
			unsigned int n_feature_edges = _m.find_feature_edges(_angle_tresh * gfxDEGTORAD); 
			if (_debug) std::cout << "Num of feature edges(not include boundaries) " << n_feature_edges << std::endl; //"\t";

			int s = 0, d = 0, cr = 0, co = 0, b = 0; // for test 
			for (v_it = _m.vertices_begin(); v_it != v_end; ++v_it) {
				if (_m.status(v_it).deleted() == false) { // ���������, ��ʱ���ǵ��������ԭ��.

					if (_m.is_boundary(v_it)) {
						b++;
					} 
					unsigned int n_fe = 0; // ����������ܵ�feature edges�ĸ���
					for (typename MeshT::VEIter ve_it = _m.ve_iter(v_it.handle()); ve_it; ++ve_it) {
						if (_m.status(ve_it).feature()) n_fe++;//
						//if (_m.status(ve_it).feature() || _m.is_boundary(ve_it)) n_fe++;//�߽��Ҳ����������.
						//Ϊ�˳���ļ���, boundary edges and vertices ��������.
					}
					if (n_fe == 0) {
						_m.property(_vp_type, v_it) = DGP::SMOOTH_VFT;	s++;	
					} else if (n_fe == 1) {
						_m.property(_vp_type, v_it) = DGP::DART_VFT;	d++;	
					} else if (n_fe == 2) {
						_m.property(_vp_type, v_it) = DGP::CREASE_VFT;	cr++;
					} else if (n_fe >= 3) {
						_m.property(_vp_type, v_it) = DGP::CORNER_VFT;	co++;	
					}
				}
			}////
			for (typename MeshT::EdgeIter e_it(_m.edges_begin()), e_end(_m.edges_end()); e_it != e_end; ++e_it) {
				//simplified_mesh_.status(e_it).set_feature(_m.status(e_it).feature());//copy the _m's feature properties
			}
			if (_debug) std::cout << "Num of feature vertex: SMOOTH_VFT " << s << ", DART_VFT " << d 
				<< ", CREASE_VFT " << cr << ", CORNER_VFT " << co << "; ����BOUNDARY " << b << ". \n"; 
			//std::cout << "--\n";// ���߽��Ҳ��������crease edge������֮��, s+d+cr+co == n_vertices();
			return n_feature_edges;
	}

} // end of namespace DGP
#endif // hoppe94loopsubvft_h