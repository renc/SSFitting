#ifndef hoppe94loopsubvft_h
#define hoppe94loopsubvft_h

namespace DGP {

	// vertex feature types 顶点的尖锐特征类型, 是根据hoppe94[Piecewise Smooth Surface Reconstruction]来定义的.
	// 放在独立的一个头文件是为了别的程序统一使用.
	// Usage: DGP::SMOOTH_VFT 而不需要 DGP::v_feature_type::SMOOTH_VFT, DGP::v_feature_type::SMOOTH_VFT.
	enum v_feature_type { SMOOTH_VFT =0, DART_VFT = 1, CREASE_VFT = 2, CORNER_VFT = 3 }; 


	template <typename MeshT> 
	inline unsigned int identify_all_sharp_features(MeshT &_m, OpenMesh::VPropHandleT<v_feature_type> &_vp_type, 
		bool _debug = false, double _angle_tresh = 44.0) {
			// 找出并标记特征边和特征点
			//         边的类型            顶点类型
			// edge -> crease edge		-> dart				其实也可以用OpenMesh中的feature状态,和corner一样,都是简化时候不能动的.
			//							-> crease,			
			//							-> corner vertex.	其实也可以用OpenMesh中的feature状态,和dart一样,都是简化时候不能动的.
			//		-> boundary edge	-> boundary vertex.
			// OpenMesh的bool is_boundary(VertexHandle/EdgeHandle/HalfedgeHandle/FaceHandle )可以测试边界性了,所以这里可以不用检测了

			// 一个property还不是valid的话就表明它还没有被add.
			if (_vp_type.is_valid() == false) {
				std::cout << "Error: DGP::identify_all_sharp_features: \n"
					<< "\t判断点类型前_vp_type is not valid, make sure mesh.add_property(_vp_type) before.\n";
			}

			typename MeshT::VIter v_it, v_end(_m.vertices_end());
			double const gfxDEGTORAD =  0.01745329251994329577;//角度单位的转换 pi/180 degree to radian
			unsigned int n_feature_edges = _m.find_feature_edges(_angle_tresh * gfxDEGTORAD); 
			if (_debug) std::cout << "Num of feature edges(not include boundaries) " << n_feature_edges << std::endl; //"\t";

			int s = 0, d = 0, cr = 0, co = 0, b = 0; // for test 
			for (v_it = _m.vertices_begin(); v_it != v_end; ++v_it) {
				if (_m.status(v_it).deleted() == false) { // 后来加入的, 当时考虑到简化网格的原因.

					if (_m.is_boundary(v_it)) {
						b++;
					} 
					unsigned int n_fe = 0; // 这个顶点四周的feature edges的个数
					for (typename MeshT::VEIter ve_it = _m.ve_iter(v_it.handle()); ve_it; ++ve_it) {
						if (_m.status(ve_it).feature()) n_fe++;//
						//if (_m.status(ve_it).feature() || _m.is_boundary(ve_it)) n_fe++;//边界边也当成特征边.
						//为了程序的简明, boundary edges and vertices 单独处理.
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
				<< ", CREASE_VFT " << cr << ", CORNER_VFT " << co << "; 其中BOUNDARY " << b << ". \n"; 
			//std::cout << "--\n";// 当边界点也被当作是crease edge来计算之后, s+d+cr+co == n_vertices();
			return n_feature_edges;
	}

} // end of namespace DGP
#endif // hoppe94loopsubvft_h