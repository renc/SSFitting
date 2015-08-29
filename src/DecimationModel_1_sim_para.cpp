#include "DecimationModel.h" 
//#include "../../../src/CourseExamples/04-Fairing/TaucsSolver.hh"

#include <algorithm> // for max_element algorithm

#include <assert.h>
#include <fstream>
int fandisk_arr[] = {
	9, 269, 265,
	261, 257, 254,
	250, 246, 241,
	237, 233, 229,
	225, 221, 217,
	214, 210, 206,
	//202, 995, 9762

	9922, 9925, 9928,
	9931, 9934, 9939,
	9943, 9947, 9951,
	9954, 9958, 9962,
	9966, 9970, 9974,
	9978, 9982, 9986,
	9990, 9994, 9998,
	10002, 10006, 10011,
	10014, 784, 788

}; 
std::vector<int> fandisk(fandisk_arr, fandisk_arr+18+27);
void DecimationModel::add_required_properties() {
	// ����ֻ���mesh_.add_property().
	// used by the whole 
	mesh_.add_property(vp_type_);//���������

	// only used by the simplification
	mesh_.add_property(vp_quadric_);//�����Quadric
	mesh_.add_property(vp_cost_);//�����ڶ����е����ȼ�
	mesh_.add_property(vp_target_);

	// only convenience the vertex hierarchy
	mesh_.add_property(heovh_);
	mesh_.add_property(vp_node_handle_);
	mesh_.add_property(vp_leaf_node_handle_);

	// used by parameterization
	mesh_.add_property(fvset_);
	mesh_.add_property(vf_);
	mesh_.add_property(vf0_);mesh_.add_property(vf1_);mesh_.add_property(vf2_);
	mesh_.add_property(vbc_);
	mesh_.add_property(vuv_);

	mesh_.add_property(hep_heset_);
	mesh_.add_property(vp_feature_edge_);	
}
void DecimationModel::remove_required_properties() {
	mesh_.remove_property(vp_quadric_);//�����Quadric
	mesh_.remove_property(vp_cost_);//�����ڶ����е����ȼ�
	mesh_.remove_property(vp_target_);

	mesh_.remove_property(heovh_);
}

//-----------------------------------------------------------------------------
//��enqueue_vertex��������
// calculate and return priority: the smaller the better. ������������.
double DecimationModel::priority_qem(TriMesh::HalfedgeHandle _heh)
{	// use quadrics to estimate approximation error, approximation cost.
	// (from, to) -> to
	TriMesh::VertexHandle from(mesh_.from_vertex_handle(_heh));
	TriMesh::VertexHandle to(mesh_.to_vertex_handle(_heh));

	DGP::Quadricd q = mesh_.property(vp_quadric_, from);  
	q += mesh_.property(vp_quadric_, to);//q�����������߶�Ӧ�ĵ�Ե�Quadric֮��

	if (vertex_optimal_placement == true) { //��������һ��Ŀ��������һ�������������.
		double newvp0 = 0.0, newvp1 = 0.0, newvp2 = 0.0;
		if (q.optimal_placement(newvp0, newvp1, newvp2)) {//
			return q(TriMesh::Point(newvp0, newvp1, newvp2));
		}
	}

	// ���Կ�������֮��������ε�quality.
	// collect faces ��ߵĶ�Ӧ����
	std::vector<TriMesh::HalfedgeHandle> loop_h;
	TriMesh::HHandle hh = _heh;
	do {
		loop_h.push_back(hh);
		hh = mesh_.opposite_halfedge_handle(mesh_.prev_halfedge_handle(hh));
	} while(hh != _heh);
	unsigned int n = loop_h.size();
	TriMesh::FaceHandle fl = mesh_.face_handle(_heh);
	TriMesh::FaceHandle fr = mesh_.face_handle(mesh_.opposite_halfedge_handle(_heh));
	// backup point positions ������������ֵ�����޸�	
	TriMesh::Point p0 = mesh_.point(from);
	TriMesh::Point p1 = mesh_.point(to);
	mesh_.set_point(from, p1);

	double aspect_ratio(0); int n_count(0);// (1) ��������״�����̱ߵı���
	const double av_ang = 60*M_PI/180;
	double angle_deviation = 60*M_PI/180;//(2)ƽ���Ƕ���60���Ĳ��
	double tri_area = 0.0; // (3)trianle area
	unsigned int fvset_size = 1;//(4)��������������ϵĶ�������,ע�ⲻҪ��ʼ��Ϊ0,����������ܳ���0, ����������񲢲���ʵ��.
	double area_ratio = 1;//(5) �������������
	double valence_deviation = 0.0;//(6)���Ϊ6�Ĳ��
	for (unsigned int i = 1; i <= n-2; ++i) { // ����fl��fr�����������ε�.
		TriMesh::FHandle fh = mesh_.face_handle(loop_h[i]);
		TriMesh::FaceHalfedgeIter fh_it(mesh_, fh);
		double len0 = mesh_.calc_edge_length(fh_it); ++fh_it;
		double len1 = mesh_.calc_edge_length(fh_it); ++fh_it;
		double len2 = mesh_.calc_edge_length(fh_it);
		double maxlen = len0, minlen = len0;
		if (len1 > maxlen) maxlen = len1; if (len2 > maxlen) maxlen = len2;  
		if (len1 < minlen) minlen = len1; if (len2 < minlen) minlen = len2;  
		//aspect_ratio += maxlen / minlen; 
		n_count++;
		if (maxlen / minlen > aspect_ratio) aspect_ratio = maxlen / minlen; //get the max aspect ratio

		TriMesh::VHandle v0 = from, v1 = mesh_.to_vertex_handle(loop_h[i]), v2 = mesh_.to_vertex_handle(loop_h[i+1]);
		OpenMesh::Vec3d v1v0 = mesh_.point(v1) - mesh_.point(v0);
		OpenMesh::Vec3d v2v0 = mesh_.point(v2) - mesh_.point(v0);
		double ang = acos(dot(v1v0, v2v0) / (v1v0.norm() * v2v0.norm()));
		if (fabs(ang - av_ang) > fabs(angle_deviation - av_ang)) angle_deviation = ang;
		tri_area = cross(v1v0, v2v0).norm() * 0.5;

		fvset_size += mesh_.property(fvset_, fh).size();		
	}
	if (mesh_.face_handle(loop_h[0]).is_valid()) fvset_size += mesh_.property(fvset_, mesh_.face_handle(loop_h[0])).size();
	if (mesh_.face_handle(loop_h[n-1]).is_valid()) fvset_size += mesh_.property(fvset_, mesh_.face_handle(loop_h[n-1])).size();

	//aspect_ratio /= n_count; 
	//angle_deviation /= n_count;
	if (angle_deviation > av_ang) angle_deviation /= av_ang;
	else angle_deviation = av_ang / angle_deviation;

	valence_deviation = (mesh_.valence(from) + mesh_.valence(to));
	if (mesh_.is_boundary(_heh) || mesh_.is_boundary(mesh_.opposite_halfedge_handle(_heh))) {
		valence_deviation -= 3;
	} else valence_deviation -= 4;
	//if (valence_deviation > 6) { if (valence_deviation > 8) { valence_deviation *= 2; } valence_deviation /= 6; } 
	//else valence_deviation = 6 / valence_deviation;//�����̫���valence��һ���Ŵ���۵�����
	if (valence_deviation != 6) valence_deviation *= 100;
	//std::cout << aspect_ratio << ", " << q(mesh_.point(to)) << ", ";
	// undo simulation changes
	mesh_.set_point(from, p0);

	double weight = 1//aspect_ratio;//Ч������

		// fvset_size * (angle_deviation + aspect_ratio);
		////* valence_deviation 
		//* area_ratio//
		//* aspect_ratio 
		//* angle_deviation //һֱ����ȥ��.
		;
	//fvset_size * angle_deviation * aspect_ratio;
	//(angle_deviation + aspect_ratio);//fandisk.off��Щ�������
	//angle_deviation;//fandisk.offЧ������

	return weight  * q(mesh_.point(to));
	//if (fabs(q(mesh_.point(to))) < 1.e-6) { return aspect_ratio; }
	//else return aspect_ratio * q(mesh_.point(to)); //fandisk.off��Ч��Ҳ����.
	//return 0.5*aspect_ratio + 0.5*q(mesh_.point(to));	// +����Ч������.
	//return q(mesh_.point(to));// evaluate quadric Q at vector v: v*Q*v, ���error cost
	// ���������Ҫ�ǵ���������ɵ��µ�����Ϊ�������۵�v��,����from��to, ��q(NewPoint(from and to)), ��, ������������.
	// ��������ʹ�õ���halfedge collapse, ��QEM��ԭʼʹ�õ�vertex pair collpase�е�����.
}

//-----------------------------------------------------------------------------
//��decimate��������
//����������(��ĳһ�������)����ʱ��Ĵ���, ����Ϊ��queque�е���������,
//�����������С�Ĵ��۵��ǰ�ߴ������Ա�����.
void DecimationModel::enqueue_vertex(TriMesh::VertexHandle _vh)
{
	float	 min_prio(FLT_MAX);//floatֵԽ��Ҳ���Ǵ���Խ��, ���ȼ�ԽС
	TriMesh::HalfedgeHandle  min_hh;//�����������ص���С���۵İ��, 

	// find best out-going halfedge, _vh����Ϊ���ܵĳ���ߵ�����
	for (TriMesh::VOHIter vh_it(mesh_, _vh); vh_it; ++vh_it) // �����������ܵĳ����
	{
		if (DGP::QEMDecimationT<TriMesh>::is_collapse_legal(mesh_, vh_it, vp_type_))
		{	//��������, Ҳ�����������ߵ���������(һ���)��������
			float prio = priority_qem(vh_it);//��һ�����(��������ɵİ�߱�ʾ)����ʱ��Ĵ���
			if (prio != -1.0 && prio < min_prio)
			{
				min_prio = prio;
				min_hh   = vh_it.handle();
			}
		}
	} //��������������ص����е�����ĸ���ԵĴ�����С(���ȼ����)��Ϊmin_prio, �������ð�߱�ʾ
	//��ʵ�����ԵĴ���,Ҳ���� float priority(HalfedgeHandle )������ֻ���������ߵ�to vertex(target vertex)�Ĵ���

	// update queue
	if (mesh_.property(vp_cost_, _vh) != -1.0) 
	{
		queue.erase(_vh);
		mesh_.property(vp_cost_, _vh) = -1.0;
	}

	if (min_hh.is_valid()) 
	{
		mesh_.property(vp_cost_, _vh) = min_prio;//�������������queue�������.
		mesh_.property(vp_target_, _vh)   = min_hh;
		queue.insert(_vh);//�����������_vh, ��Ҫʹ�õ�mesh_.property(vp_cost_, _vh)���������ڶ����е�λ��
		//�����Ļ��Ӷ�����������С���۵ĵ�_vh,��mesh_.property(vp_target_, _vh)�л�����Ӧ��min_hh
	}
}

void DecimationModel::initial_parameterization(TriMesh::HalfedgeHandle _hh) {
	// (from, to)->to
	TriMesh::VertexHandle from = mesh_.from_vertex_handle(_hh);
	TriMesh::VertexHandle to = mesh_.to_vertex_handle(_hh);
	TriMesh::HalfedgeHandle o = mesh_.opposite_halfedge_handle(_hh);
	//std::cout << "Para (" << from << ", " << to <<")\n";

	// ������������������, ��һ�ְ��_hh����԰�߶����Ǳ߽��,Ҳ����˵from������interior vertex,
	// �ڶ��������_hh�Ǳ߽��, from�����Ǳ߽綥��, �����������_hh�Ķ԰���Ǳ߽��, ��ʱfrom����Ҳ�Ǳ߽��.
	if (!mesh_.is_boundary(_hh) && !mesh_.is_boundary(o)) { //�Ǳ߽����, from��interior vertex
		// ��һ��case
		//std::cout << "c1\n"; // for test 
		// (from, to)->to, ��from����һ�����ϵĵ��������һ����λԲ��, ��Mean Value Mapping. 
		// ��from���㿴����vi����, ��from�����one-ring������vj����

		std::vector<TriMesh::VertexHandle> loop; //����ѭ��ؼ�¼vi����(����Ҳ����from����)��one-ring��������Ķ���
		std::vector<TriMesh::HalfedgeHandle> loop_h;//���ڼ�¼vi����(����Ҳ����from����)��one-ring��������ĳ����
		std::vector<TriMesh::FHandle> loop_f;
		TriMesh::HalfedgeHandle h = _hh;
		do {	// ���ﲻʹ��VertexOHalfedgeIter��VertexVertexIter��ԭ������Ҫ��֤loop[j]��˳��.loop[0]=to
			loop_h.push_back(h);
			loop_f.push_back(mesh_.face_handle(h));
			loop.push_back(mesh_.to_vertex_handle(h));//std::cout << mesh_.to_vertex_handle(h) << "\t";	 ////	for test	
			h = mesh_.opposite_halfedge_handle(mesh_.prev_halfedge_handle(h));
		} while(h != _hh);
		if (loop.size() != loop_h.size()) { std::cerr << "Error: collect one-ring vertices and halfedges.\n"; };

		unsigned int n = loop_h.size();//�ж��ٸ��߽綥��
		unsigned int i = 0;
		double length = 0.0, rou_angle = 0.0;//
		std::vector<double> vec_angle;
		for (i=0 ; i<n; ++i) {//����ܵĳ���(�ܳ�), �Լ��Ƕ�֮��. 
			length += (OpenMesh::vector_cast<OpenMesh::Vec3d>(mesh_.point(loop[i])) - OpenMesh::vector_cast<OpenMesh::Vec3d>(mesh_.point(loop[(i+1)%n]))).norm();

			OpenMesh::Vec3d d1 = (OpenMesh::vector_cast<OpenMesh::Vec3d>(mesh_.point(loop[i])) - OpenMesh::vector_cast<OpenMesh::Vec3d>(mesh_.point(from))).normalize();
			OpenMesh::Vec3d d2 = (OpenMesh::vector_cast<OpenMesh::Vec3d>(mesh_.point(loop[(i+1)%n])) - OpenMesh::vector_cast<OpenMesh::Vec3d>(mesh_.point(from))).normalize();
			double angle = atan2(cross(d1, d2).norm(), dot(d1, d2)); //acos(dot(d1, d2));
			rou_angle += angle;
			vec_angle.push_back(angle);
		}  //std::cout << "angle " << rou_angle << ".\n";//std::cout << "len: " << length << std::endl; //for test	
		if (vec_angle.size() != n) { std::cout << "Error: ������һ��.\n"; } // 1000

		bool is_crease_edge = false; //���ж�Ҫ������_hh�ǲ���crease edge, Ҳ����from�����ǲ���crease vertex, 
		if (mesh_.property(vp_type_, from) == DGP::CREASE_VFT) { // if (mesh_.status(mesh_.edge_handle(_hh)).feature() == true)��һ����.
			is_crease_edge = true;
		}
		int feature_vertex_index = -1;// ��is_crease_edge == true;���loop[]��������һ����������.
		if (is_crease_edge == true) {//�ҳ�feature_vertex_index���ڶ���.
			for (i = 2; i < n-1; ++i) //��һ��crease edgeֻ������loop_h[2]��loop_h[n-2]�ⷶΧ֮��
				if (mesh_.status(mesh_.edge_handle(loop_h[i])).feature() == true) 
					feature_vertex_index = i;
			if (feature_vertex_index == -1) std::cerr << "Error: initial_parameterization(): feature vertex index.\n";
		} // ֻҪis_crease_edge == true, ��ô
		double bcv, bt;//��is_crease_edge == trueʱ���������, �������������С����, ʵ�����Ǳ���.

		//---------------------------------------------------------------------------------------------------
		/*// ��һ ��λԲ
		double l = 0.0, angle = 0.0;
		// to decide the coordiante of vertex from.
		if (is_crease_edge == false) // using the mean value coordiante
		{
		//std::vector<OpenMesh::Vec2d > h_vj(loop_h.size());//vj����Ĳ���������
		std::vector<double> wij(loop_h.size());//���Ƕ���vi��vj��ϵ��
		double sum_wij = 0.0;
		OpenMesh::Vec2d sum_wh(0, 0, 0);
		for (i = 0; i < n; ++i) { //����vi��one-ring�ϵ���n���ڽӵ�
		// fix the boundary/one-ring vertices, 
		angle = l / length * (2.0*M_PI); //std::cout << "ang: " << angle << std::endl;
		mesh_.property(vuv_, loop[i]) = OpenMesh::Vec2d(0.5*cos(angle)+0.5, 0.5*sin(angle)+0.5);
		//�����λԲ����(0.5, 0.5)ΪԲ��, �뾶��0.5. Բ��Ҫ�Ƕ���(0, 0)�Ļ�, ���籾����(0, 0.5)�������������Ϊ(e, 0.5),e��һ����С��С��ֵ.
		l += (mesh_.point(loop[i]) - mesh_.point(loop[(i+1)%n])).norm();

		// ���vi��vj�����ϵ�ϵ��, ע������ʹ�õ���mean value coordiante, ������harmonic map�е�cotangent weight
		TriMesh::VertexHandle v0 = from; // = mesh_.from_vertex_handle(loop_h[i]);
		TriMesh::VertexHandle v1 = loop[i];// =mesh_.to_vertex_handle(loop_h[i]);
		TriMesh::VertexHandle v2 = loop[(i+1)%n];//= mesh_.to_vertex_handle(mesh_.next_halfedge_handle(loop_h[i]));//������������ȷ��ǰ���ǷǱ߽����
		TriMesh::VertexHandle v3 = loop[(i+n-1)%n];//mesh_.to_vertex_handle(mesh_.next_halfedge_handle(mesh_.opposite_halfedge_handle(loop_h[i])));
		OpenMesh::Vec3d v1v0 = mesh_.point(v1) - mesh_.point(v0);	// ֮ǰ������ֵ�-1
		OpenMesh::Vec3d v2v0 = mesh_.point(v2) - mesh_.point(v0);
		OpenMesh::Vec3d v3v0 = mesh_.point(v3) - mesh_.point(v0);
		double v1v0v2_angle = acos(dot(v1v0, v2v0) / (v1v0.norm() * v2v0.norm()));//ʸ��v1v0��v2v0�ļн�.
		double v3v0v1_angle = acos(dot(v3v0, v1v0) / (v3v0.norm() * v1v0.norm()));//ʸ��v3v0��v1v0�ļн�.
		wij[i] = (tan(v1v0v2_angle / 2) + tan(v3v0v1_angle / 2)) / v1v0.norm();
		if (wij[i] < 0) { std::cout << "Error: wij is negative.\n"; return; }

		sum_wij += wij[i];
		sum_wh += wij[i] * mesh_.property(vuv_, loop[i]);
		}
		mesh_.property(vuv_, from) = sum_wh * (1.0/sum_wij);//�����vi�����ڵ�λԲ�ϵĲ���������.
		// ���������from��������һ�����ڵĲ�����uv����.
		} else {
		for (i = 0; i < n; ++i) { //����vi��one-ring�ϵ���n���ڽӵ�
		// fix the boundary/one-ring vertices, 
		angle = l / length * (2.0*M_PI); //std::cout << "ang: " << angle << std::endl;
		mesh_.property(vuv_, loop[i]) = OpenMesh::Vec2d(0.5*cos(angle)+0.5, 0.5*sin(angle)+0.5);
		l += (mesh_.point(loop[i]) - mesh_.point(loop[(i+1)%n])).norm();
		} 
		double len_to_from = (double)(mesh_.point(to) - mesh_.point(from)).norm();
		double len_creasevertex_from = (double)(mesh_.point(loop[feature_vertex_index]) - mesh_.point(from)).norm();
		bcv = len_to_from / (len_to_from + len_creasevertex_from);//barycentric coordinate������ڶ���crease vertex�ķ���.
		bt = 1- bcv;
		mesh_.property(vuv_, from) = (mesh_.property(vuv_, to) - mesh_.property(vuv_, loop[feature_vertex_index])) * bt + mesh_.property(vuv_, loop[feature_vertex_index]);
		}*/
		// ���� ѹƽ
		double angle_scale_ratio = 2 * M_PI / rou_angle; //���ű���
		double temp_sum_angle = 0.0, len = 0;
		for (int i = 0; i < n; ++i) { //һ�����Ĳ���������.
			temp_sum_angle += (vec_angle[i] * angle_scale_ratio);
			len = (OpenMesh::vector_cast<OpenMesh::Vec3d>(mesh_.point(loop[(i+1)%n])) - OpenMesh::vector_cast<OpenMesh::Vec3d>(mesh_.point(from))).norm(); 
			len *= angle_scale_ratio;
			mesh_.property(vuv_, loop[(i+1)%n]) = OpenMesh::Vec2d(len*cos(temp_sum_angle), len*sin(temp_sum_angle));
		} 		
		if (feature_vertex_index != -1) {
			double len_to_from = (double)(mesh_.point(to) - mesh_.point(from)).norm();
			double len_creasevertex_from = (double)(mesh_.point(loop[feature_vertex_index]) - mesh_.point(from)).norm();
			bcv = len_to_from / (len_to_from + len_creasevertex_from);//barycentric coordinate������ڶ���crease vertex�ķ���.
			bt = 1- bcv; //from��������Ĳ�������ֱ�Ӹ���, �����������ĸ���ͺ���ͳһ����.
			mesh_.property(vuv_, from) = mesh_.property(vuv_, to)
				+ (mesh_.property(vuv_, loop[feature_vertex_index]) - mesh_.property(vuv_, to)) * bcv;
		} else { // not crease 
			mesh_.property(vuv_, from) = OpenMesh::Vec2d(0, 0); //���һ��������� 
		} 

		//---------------------------------------------------------------------------------------------------		

		//�ռ����е�free vertices, ��Щfree vertices��������������.
		//�ֱ���from����, �Լ�from����һ�����n����������������֮ǰ�Ĳ����������а����Ķ���.
		// ����(1)����from�����Ƿ�crease vertex�����뵽free_vertics���������������䵽��һ������,
		// ����(2)from���㲻��crease vertexʱ��ż��뵽free_vertics����, ��Ϊ�Ѿ�����ֱ�Ӹ��������������Ϣ.
		// ����ʹ�õ�������(2).
		std::vector<TriMesh::VertexHandle> free_vertices;		
		if (is_crease_edge == false) {
			free_vertices.push_back(from);
		}////
		for (i = 0; i < n; ++i) //ѭ��n����
		{	
			TriMesh::FaceHandle fh = loop_f[i];//= mesh_.face_handle(loop_h[i]);
			TriMesh::VertexHandle v0 = from;//=mesh_.from_vertex_handle(loop_h[i]);//������3������, �������from����
			TriMesh::VertexHandle v1 = loop[i];//=mesh_.to_vertex_handle(loop_h[i]);
			TriMesh::VertexHandle v2 = loop[(i+1)%n];//=mesh_.to_vertex_handle(mesh_.next_halfedge_handle(loop_h[i]));

			//ȡ�����������������֮ǰ�Ĳ��������̵��в�������������ϵĶ���, �����������ѹƽ��one-ring�еĲ�����uv����
			//��Щ��������һ������crease edge����.
			std::vector<TriMesh::VertexHandle> fvset = mesh_.property(fvset_, fh);
			for (std::vector<TriMesh::VertexHandle>::iterator it(fvset.begin()), it_end(fvset.end()); it != it_end; ++it) {
				//��ʼ��ʱ��ÿһ�����ϻ�û�в�������������Ķ���, ��ʱfvsetΪ��.
				if (mesh_.status(*it).deleted() == false) { std::cerr << "Error: \n"; return; } // ��Щ��ǰ�Ѿ����������ĵ�Ӧ�ö���deleted vertices.
				free_vertices.push_back(*it);				

				// ͬʱ�����Щfree vertices�����һ��������Ĳ�������.
				TriMesh::VertexHandle vf0 = mesh_.property(vf0_, *it); //���*it��������ǰ�Ĳ��������������ڵ������������,
				TriMesh::VertexHandle vf1 = mesh_.property(vf1_, *it); //���������������ɵ���, �պ�Ӧ�þ�������ѭ���е��Ǹ���. 
				TriMesh::VertexHandle vf2 = mesh_.property(vf2_, *it); //����������֮����һ������from����.
				if (vf0 != v0 && vf0 != v1 && vf0 != v2) { std::cerr << "Error: not match vf0.\n"; return; }// 
				if (vf1 != v0 && vf1 != v1 && vf1 != v2) { std::cerr << "Error: not match vf1.\n"; return; }// 
				if (vf2 != v0 && vf2 != v1 && vf2 != v2) { std::cerr << "Error: not match vf2.\n"; return; }//
				mesh_.property(vuv_, *it) = mesh_.property(vuv_, vf0) * (mesh_.property(vbc_, *it))[0]
				+ mesh_.property(vuv_, vf1) * (mesh_.property(vbc_, *it))[1]
				+ mesh_.property(vuv_, vf2) * (mesh_.property(vbc_, *it))[2];
			}

			(mesh_.property(fvset_, fh)).clear();//���������Ķ��������, ���������¿�����Щ�����������������.
		}

		// ��initial parameterization�м�����һ��ƽ��,�����и�������ƽ���������Щ���㱻����from�����һ������֮����.
		for (std::vector<TriMesh::VertexHandle>::iterator freev_it = free_vertices.begin(), freev_end = free_vertices.end(); freev_it != freev_end; ++freev_it) 
		{ //freev_it(from����),������ƽ����. 
			if (*freev_it == from) continue;
			if (mesh2_.is_boundary(*freev_it)) continue;//����Ǳ߽��ϵĶ���(��������ʱ��Ҳ���������߽��ϵ�)�Ͳ�Ҫ����smooth parameterize.
			if (mesh_.property(vp_type_, *freev_it) == DGP::CREASE_VFT)  continue;					 

			bool one_ring_vertex_is_valid = true; //�ȼ��趥��*freev_it��һ���򶥵㶼���������Բ��, �������Ч.
			double sum_wij = 0.0; OpenMesh::Vec2d sum_wh(0, 0, 0);
			for (TriMesh::VOHIter voh_it(mesh2_, *freev_it); voh_it; ++voh_it) 
			{ 
				TriMesh::VHandle vj = mesh2_.to_vertex_handle(voh_it.handle());
				std::vector<TriMesh::FHandle>::iterator result = find(loop_f.begin(), loop_f.end(), mesh_.property(vf_, vj)); 
				if (result != loop_f.end()) { //�������vj�Ƿ����������Բ��?
					sum_wij += mesh2_.property(hmvc_, voh_it.handle());
					sum_wh += mesh2_.property(hmvc_, voh_it.handle()) * mesh_.property(vuv_, vj);
				} else { one_ring_vertex_is_valid = false; }
			} // end of for
			if (one_ring_vertex_is_valid == true) {
				mesh_.property(vuv_, *freev_it) = sum_wh * (1.0 / sum_wij);
			}
		}/**//**/

		// �����Ƕ�λ��Щfree vertices�����ڰ��_hh�������γɵ���һ��������֮��, ���_hh����֮��һ��������n��Ϊn-2��.
		// ���е�free vertices��ʼ��Ϊ��û�ж�λ��������.
		OpenMesh::VPropHandleT<bool> is_face_located;
		mesh_.add_property(is_face_located);
		for (std::vector<TriMesh::VHandle>::size_type freev_i = 0, freev_size = free_vertices.size(); freev_i < freev_size; ++ freev_i) {
			mesh_.property(is_face_located, free_vertices[freev_i]) = false;
		}

		for (i = 1 ;i <= n-2; ++i) { // ѭ��n-2����, ��n-2�����Ǽ���(from, to)->to֮���γɵ�
			// �ж���������Ҫ���¶�λ��free vertices�Ƿ������triangle(loop[0], loop[i], loop[i+1])֮��
			for (std::vector<TriMesh::VertexHandle>::iterator freev_it = free_vertices.begin(), freev_end = free_vertices.end(); freev_it != freev_end; ++freev_it) {
				if (mesh_.property(is_face_located, *freev_it) == false) {

					OpenMesh::Vec2d pp = mesh_.property(vuv_, *freev_it);//���free vertex�Ĳ�������
					OpenMesh::Vec3d bc = DGP::calc_barycentric_coordinates(mesh_.property(vuv_, loop[0]), mesh_.property(vuv_, loop[i]), mesh_.property(vuv_, loop[i+1]), pp);
					//std::cout << bc << std::endl;//for test

					//(0.5, -1.11022e-016, 0.5)//�����������������������,���е�������С��С��-1.11022e-016��5.55112e-017��Ӧ����Ϊ0��,
					//(0.5, 0.5, 5.55112e-017) //�����Ļ���һ������������ʵ�ǺϺ���һ���������ڲ�(�߽���)���ڲ����Ҫ���.
					if (fabs(bc[0]) < 1.0e-10) bc[0] = 0;//ԭ����1.0e-6,���ڸ�Ϊ1.0e-5, ������bc(-1.37927e-006 0.19194 0.808062)
					if (fabs(bc[1]) < 1.0e-10) bc[1] = 0;//�ͱ���Ϊ����Ч��. ����ʵ�ϵ�һ���Ϊ0, �������Ϊ1.
					if (fabs(bc[2]) < 1.0e-10) bc[2] = 0;//����(0.303709 -0.000569765 0.696861)

					//if (!(bc[0]<0) && !(bc[1]<0) && !(bc[2]<0)) { //Ҫ�����bc��������������>=0
					if ((bc[0] >= 0) && (bc[1] >= 0) && (bc[2] >= 0)) { //Ҫ�����bc��������������>=0

						mesh_.property(is_face_located, *freev_it) = true;

						// For interior vertex, its' barycentric coordinates is bc according to the triangle(0, i, i+1).
						//std::cout << *freev_it << ": " << loop[0] << "," << loop[i] << "," << loop[i+1] << ": " << bc << "\n"; // for test

						TriMesh::FaceHandle fh = mesh_.face_handle(mesh_.find_halfedge(loop[i], loop[i+1]));
						if (fh.is_valid() == false) std::cerr << "Error: face invalid." << std::endl;//Only for assert

						// �������о���ÿһ�β�������Ӧ�û�õĻ��Ǹ��ĵ���Ϣ.
						mesh_.property(vf_,  *freev_it) = fh;
						mesh_.property(vf0_, *freev_it) = loop[0]; mesh_.property(vf1_, *freev_it) = loop[i]; mesh_.property(vf2_, *freev_it) = loop[i+1];
						mesh_.property(vbc_, *freev_it) = bc;

						(mesh_.property(fvset_, fh)).push_back(*freev_it); 
					}
				} // ���free vertex�����������
			} // һ����free vertex �����������			
		} // ���е�free vertices ������n-2������
		for (std::vector<TriMesh::VHandle>::size_type freev_i = 0, freev_size = free_vertices.size(); freev_i < freev_size; ++ freev_i) { 
			//���������жϵ�����,���ƿ���ȥ����.
			if (mesh_.property(is_face_located, free_vertices[freev_i]) == false) { //��ʾ���free vertex�������n-2������û���ҵ�ͶӰ��������һ����, ���Ǵ����.
				std::cerr << "Error: free vertex reloacted 1. " << free_vertices[freev_i] << "; " << from << ", " << to << "\n";// 
				return;
			}  
		}
		mesh_.remove_property(is_face_located);

		if (is_crease_edge == true) {
			// ����, crease vertex���͵�from����û�м��뵽free_vertcies������, from�����������Ϣ��û������.
			mesh_.property(vf_,  from) = loop_f[feature_vertex_index - 1]; // = mesh_.face_handle(loop_h[feature_vertex_index - 1]);������������Ĭ��ѡ����Ǹ�
			mesh_.property(vf0_, from) = loop[0]; mesh_.property(vf1_, from) = loop[feature_vertex_index - 1]; mesh_.property(vf2_, from) = loop[feature_vertex_index];
			mesh_.property(vbc_, from) = OpenMesh::Vec3d(bt, 0, bcv);

			(mesh_.property(fvset_, loop_f[feature_vertex_index-1])).push_back(from);

			mesh_.property(vp_feature_edge_, from) = mesh_.edge_handle(loop_h[feature_vertex_index]);//���from�㱻���������Ǹ���������.

			// ������߼�¼����ʲô���˳�����, �����for���Ӧ����copy�㷨��������
			// ����
			for (std::vector<TriMesh::HalfedgeHandle>::iterator it(mesh_.property(hep_heset_, _hh).begin()), end(mesh_.property(hep_heset_, _hh).end()); it != end; ++it)
			{ //���������İ��_hh�ϵ�֮ǰ���������������ߵļ�¼���Ƶ��µ���������.
				mesh_.property(hep_heset_, mesh_.opposite_halfedge_handle(loop_h[feature_vertex_index])).push_back(*it);

				if (mesh2_.to_vertex_handle(*it) != to) { //������������_hh���ϵĵ�, �����µ��µı���
					mesh_.property(vp_feature_edge_, mesh2_.to_vertex_handle(*it)) = mesh_.edge_handle(loop_h[feature_vertex_index]);
				}
			}
			// ����. ������������_hh�ķ�����ϵĵ�, ����û�и��µ��µı�������Ϊ������ʵ�Ѿ��ﵽЧ����, ����һ����.
			std::vector<TriMesh::HalfedgeHandle> tmp(mesh_.property(hep_heset_, mesh_.opposite_halfedge_handle(_hh)));
			for (std::vector<TriMesh::HalfedgeHandle>::iterator it(mesh_.property(hep_heset_, loop_h[feature_vertex_index]).begin()), end(mesh_.property(hep_heset_, loop_h[feature_vertex_index]).end()); it != end; ++it)
			{
				tmp.push_back(*it);
			}
			mesh_.property(hep_heset_, loop_h[feature_vertex_index]).clear();
			for (std::vector<TriMesh::HalfedgeHandle>::iterator it(tmp.begin()), end(tmp.end()); it != end; ++it)
			{
				mesh_.property(hep_heset_, loop_h[feature_vertex_index]).push_back(*it);
			}			
		} // end of if

	} else if (mesh_.is_boundary(_hh) && !mesh_.is_boundary(o)) {
		// �ڶ���case
		//std::cout << "c2\n"; // for test // means in case 2
		std::vector<TriMesh::VertexHandle> loop; //�ռ�from�����one-ring����vj// ��ʼʱ��vj = to, loop[0] = to 
		std::vector<TriMesh::HalfedgeHandle> loop_h;
		TriMesh::HalfedgeHandle h = _hh;
		TriMesh::VertexHandle vj = mesh_.to_vertex_handle(h);
		do {
			loop_h.push_back(h);
			loop.push_back(vj); // std::cout << "vj: " << vj << "\n";// for test
			h = mesh_.next_halfedge_handle(mesh_.opposite_halfedge_handle(h));
			vj = mesh_.to_vertex_handle(h);
		} while(vj != to);
		assert(loop.size() == loop_h.size());

		// map boundary loop to unit circle 2D domain
		unsigned int i, n = loop.size();//�߽�����n������
		TriMesh::Scalar  angle, l, length;
		double len_to_from = (double)(mesh_.point(to) - mesh_.point(from)).norm();
		double len_vn_1_from = (double)(mesh_.point(loop[n-1]) - mesh_.point(from)).norm();

		for (i = 0, length = 0.0; i < n-1; ++i) //����ܵĳ���
			length += (mesh_.point(loop[i]) - mesh_.point(loop[(i+1)%n])).norm();
		length += (len_to_from + len_vn_1_from);

		for (i=0, l=0.0; i<n; ++i) { //fix the boundary vertices		
			angle = l/length * (2.0*M_PI);
			mesh_.property(vuv_, loop[i]) = OpenMesh::Vec2d(0.5*cos(angle)+0.5, 0.5*sin(angle)+0.5);//�����λԲ����(0.5, 0.5)ΪԲ��, �뾶��0.5.
			l += (mesh_.point(loop[i]) - mesh_.point(loop[(i+1)%n])).norm();
		} // ����, from�����һ���򶥵��Ѿ����������.
		for (i = 0; i < n; ++i) {
			//std::cout << mesh_.property(vuv_, loop[i]) << "\n";
		}

		// ����from����, ���� 
		double bn_1 = len_to_from / (len_to_from + len_vn_1_from);
		double bt = 1 - bn_1; // = len_vn_1_from / (len_to_from + len_vn_1_from);
		OpenMesh::Vec3d bc_from(bt, 0, bn_1);
		mesh_.property(vuv_, from) = (mesh_.property(vuv_, loop[0]) - mesh_.property(vuv_, loop[n-1])) * bt + mesh_.property(vuv_, loop[n-1]);
		//std::cout << from << ": " << to<< "," << loop[n-2] << "," << loop[n-1] << ": " << bc_from << "\n"; // for test

		// �ռ�from�����One-ring��������ϵĶ���, ��Щ��������֮ǰ�Ĳ����������в��������䵽��Ӧ���������ϵ�.
		std::vector<TriMesh::VertexHandle> free_vertices;
		OpenMesh::VPropHandleT<bool> is_face_located;
		mesh_.add_property(is_face_located);
		for (i = 1; i <= n-1; ++i) { //ѭ��from�����һ������, ��n-1��.
			TriMesh::FaceHandle fh = mesh_.face_handle(loop_h[i]);
			TriMesh::VertexHandle v0 = from;//=mesh_.from_vertex_handle(loop_h[i]);//������3������, ��һ������from����
			TriMesh::VertexHandle v1 = loop[i-1];////����������from�����һ���򶥵�,
			TriMesh::VertexHandle v2 = loop[i];//mesh_.to_vertex_handle(loop_h[i]);//����������Ĳ��������������.

			std::vector<TriMesh::VertexHandle> fvset = mesh_.property(fvset_, fh);
			for (std::vector<TriMesh::VertexHandle>::iterator it(fvset.begin()), it_end(fvset.end()); it != it_end; ++it) {
				//��ʼ��ʱ��ÿһ�����ϻ�û�в�������������Ķ���, ��ʱfvsetΪ��.
				if (mesh_.status(*it).deleted() == false) { std::cerr << "Error: \n"; return; } // ��Щ��ǰ�Ѿ����������ĵ�Ӧ�ö���deleted vertices.
				free_vertices.push_back(*it);

				TriMesh::VertexHandle vf0 = mesh_.property(vf0_, *it); //���*it��������ǰ�Ĳ��������������ڵ������������,
				TriMesh::VertexHandle vf1 = mesh_.property(vf1_, *it); //���������������ɵ���, �պ�Ӧ�þ�������ѭ���е��Ǹ���. 
				TriMesh::VertexHandle vf2 = mesh_.property(vf2_, *it);
				if (vf0 != v0 && vf0 != v1 && vf0 != v2) { std::cerr << "Error: parameterziation vf0.\n"; return; }
				if (vf1 != v0 && vf1 != v1 && vf1 != v2) { std::cerr << "Error: parameterziation vf1.\n"; return; }
				if (vf2 != v0 && vf2 != v1 && vf2 != v2) { std::cerr << "Error: parameterziation vf2.\n"; return; }
				mesh_.property(vuv_, *it) = mesh_.property(vuv_, vf0) * (mesh_.property(vbc_, *it))[0]
				+ mesh_.property(vuv_, vf1) * (mesh_.property(vbc_, *it))[1]
				+ mesh_.property(vuv_, vf2) * (mesh_.property(vbc_, *it))[2];//���������ɶ�����µĲ���������.
			}
			(mesh_.property(fvset_, fh)).clear();//���������Ķ��������, ���������¿�����Щ�����������������.
		}
		for (unsigned int freev_i = 0, freev_size = free_vertices.size(); freev_i < freev_size; ++ freev_i) {
			mesh_.property(is_face_located, free_vertices[freev_i]) = false;//���е�free vertices��ʼΪ��û�ж�λ������.
		}

		for (unsigned int j = 1; j <= n-2; ++j) { //ѭ��n-2����, ��Щ������(from, to)->to֮��ȥ��fl��fr���γɵ�.
			TriMesh::VertexHandle v0 = loop[0];//=to
			TriMesh::VertexHandle vj = loop[j];
			TriMesh::VertexHandle vj1 = loop[j+1];//�����������γ�һ��������, �����ж�free_vertices��������Ķ����Ƿ���������.
			TriMesh::FaceHandle fh = mesh_.face_handle(loop_h[j+1]);

			for (std::vector<TriMesh::VertexHandle>::iterator freev_it = free_vertices.begin(), freev_end = free_vertices.end(); freev_it != freev_end; ++freev_it) {
				if (mesh_.property(is_face_located, *freev_it) ==false) {

					OpenMesh::Vec2d pp = mesh_.property(vuv_, *freev_it);//���free vertex�Ĳ�������
					OpenMesh::Vec3d bc = DGP::calc_barycentric_coordinates(mesh_.property(vuv_,loop[0]), mesh_.property(vuv_, loop[j]), mesh_.property(vuv_, loop[j+1]), pp);
					//std::cout << bc << std::endl;//for test
					//(0.5, -1.11022e-016, 0.5)//�����������������������,���е�������С��С��-1.11022e-016��5.55112e-017��Ӧ����Ϊ0��,
					//(0.5, 0.5, 5.55112e-017) //�����Ļ���һ������������ʵ�ǺϺ���һ���������ڲ�(�߽���)���ڲ����Ҫ���.
					if (fabs(bc[0]) < 1.0e-10) bc[0] = 0;
					if (fabs(bc[1]) < 1.0e-10) bc[1] = 0;
					if (fabs(bc[2]) < 1.0e-10) bc[2] = 0;
					if (!(bc[0]<0) && !(bc[1]<0) && !(bc[2]<0)) { //Ҫ�����bc��������������>=0

						// For interior vertex, its' barycentric coordinates is bc according to the triangle(v0, vj, vj1).
						//std::cout << *freev_it << ": " << loop[0] << "," << loop[j] << "," << loop[j+1] << ": " << bc << "\n"; // for test
						// �������о���ÿһ�β�������Ӧ�û�õĻ��Ǹ��ĵ���Ϣ.
						mesh_.property(vf_, *freev_it) = fh;
						mesh_.property(vf0_, *freev_it) = loop[0]; mesh_.property(vf1_, *freev_it) = loop[j]; mesh_.property(vf2_, *freev_it) = loop[j+1];
						mesh_.property(vbc_, *freev_it) = bc;

						(mesh_.property(fvset_, fh)).push_back(*freev_it);

						mesh_.property(is_face_located, *freev_it) = true;
					}
				} // ���free vertex�����������
			} // һ����free vertex �����������			
		}
		mesh_.property(vbc_, from) = bc_from;
		mesh_.property(vf0_, from) = to; mesh_.property(vf1_, from) = loop[n-2]; mesh_.property(vf2_, from) = loop[n-1];
		mesh_.property(vf_, from) = mesh_.face_handle(loop_h[n-1]);
		(mesh_.property(fvset_, mesh_.face_handle(loop_h[n-1]))).push_back(from);

		for (unsigned int freev_i = 0, freev_size = free_vertices.size(); freev_i < freev_size; ++ freev_i) { //���������жϵ�����,���ƿ���ȥ����.
			if (mesh_.property(is_face_located, free_vertices[freev_i]) == false) { //��ʾ���free vertex�������n-2������
				std::cerr << "Error: parameterization case 2.\n";// û���ҵ�ͶӰ��������һ����, ���Ǵ����.
				return;
			}  
		}
		mesh_.remove_property(is_face_located);

	} else if (!mesh_.is_boundary(_hh) && mesh_.is_boundary(o)) {
		// ������case
		//std::cout << "c3\n"; // for test 
		std::vector<TriMesh::VertexHandle> loop; //�ռ�from�����one-ring����vj
		std::vector<TriMesh::HalfedgeHandle> loop_h;
		TriMesh::HalfedgeHandle h = _hh;
		TriMesh::VertexHandle vj = mesh_.to_vertex_handle(h);// ��ʼʱ��vj = to
		do {
			loop_h.push_back(h);
			loop.push_back(vj);//std::cout << "vj: " << vj << "\n";// for test
			h = mesh_.opposite_halfedge_handle(mesh_.prev_halfedge_handle(h));
			vj = mesh_.to_vertex_handle(h);
		} while(vj != to);
		assert(loop.size() == loop_h.size());

		// map boundary loop to unit circle 2D domain
		unsigned int i, n = loop.size();//�߽�����n������
		TriMesh::Scalar  angle, l, length;
		double len_to_from = (double)(mesh_.point(to) - mesh_.point(from)).norm();
		double len_vn_1_from = (double)(mesh_.point(loop[n-1]) - mesh_.point(from)).norm();

		for (i=0, length=0.0; i<n-1; ++i) //����ܵĳ���
			length += (mesh_.point(loop[i]) - mesh_.point(loop[(i+1)%n])).norm();
		length += (len_to_from + len_vn_1_from);

		for (i=0, l=0.0; i<n; ++i) { //fix the boundary vertices		
			angle = l/length * (2.0*M_PI);
			mesh_.property(vuv_, loop[i]) = OpenMesh::Vec2d(0.5*cos(angle)+0.5, 0.5*sin(angle)+0.5);//�����λԲ����(0.5, 0.5)ΪԲ��, �뾶��0.5.
			l += (mesh_.point(loop[i]) - mesh_.point(loop[(i+1)%n])).norm();
		}
		for (i = 0; i < n; ++i) {
			//std::cout << mesh_.property(vuv_, loop[i]) << "\n";
		}

		// ����from����, ���� 
		double bn_1 = len_to_from / (len_to_from + len_vn_1_from);
		double bt = 1 - bn_1; // = len_vn_1_from / (len_to_from + len_vn_1_from);
		OpenMesh::Vec3d bc_from(bt, 0, bn_1);
		mesh_.property(vuv_, from) = (mesh_.property(vuv_, loop[0]) - mesh_.property(vuv_, loop[n-1])) * bt + mesh_.property(vuv_, loop[n-1]);
		//std::cout << from << ": " << to<< "," << loop[n-2] << "," << loop[n-1] << ": " << bc_from << "\n"; // for test

		// �ռ�from�����One-ring��������ϵĶ���, ��Щ��������֮ǰ�Ĳ����������в��������䵽��Ӧ���������ϵ�.
		std::vector<TriMesh::VertexHandle> free_vertices;
		OpenMesh::VPropHandleT<bool> is_face_located;
		mesh_.add_property(is_face_located);
		for (i = 0; i <= n-2; ++i) { //ѭ��from�����һ������.
			TriMesh::FaceHandle fh = mesh_.face_handle(loop_h[i]);
			TriMesh::VertexHandle v0 = from;//mesh_.from_vertex_handle(loop_h[i]);//������3������, ��һ������from����
			TriMesh::VertexHandle v1 = loop[i];//mesh_.to_vertex_handle(loop_h[i]);//����������from�����һ���򶥵�,
			TriMesh::VertexHandle v2 = loop[i+1];//mesh_.to_vertex_handle(loop_h[i+1]);//����������Ĳ��������������.

			std::vector<TriMesh::VertexHandle> fvset = mesh_.property(fvset_, fh);
			for (std::vector<TriMesh::VertexHandle>::iterator it(fvset.begin()), it_end(fvset.end()); it != it_end; ++it) {
				//��ʼ��ʱ��ÿһ�����ϻ�û�в�������������Ķ���, ��ʱfvsetΪ��.
				if (mesh_.status(*it).deleted() == false) { std::cerr << "Error: \n"; return; } // ��Щ��ǰ�Ѿ����������ĵ�Ӧ�ö���deleted vertices.
				free_vertices.push_back(*it);			

				TriMesh::VertexHandle vf0 = mesh_.property(vf0_, *it); //���*it��������ǰ�Ĳ��������������ڵ������������,
				TriMesh::VertexHandle vf1 = mesh_.property(vf1_, *it); //���������������ɵ���, �պ�Ӧ�þ�������ѭ���е��Ǹ���. 
				TriMesh::VertexHandle vf2 = mesh_.property(vf2_, *it);
				if (vf0 != v0 && vf0 != v1 && vf0 != v2) { std::cerr << "Error: parameterziation vf0.\n"; return; }
				if (vf1 != v0 && vf1 != v1 && vf1 != v2) { std::cerr << "Error: parameterziation vf1.\n"; return; }
				if (vf2 != v0 && vf2 != v1 && vf2 != v2) { std::cerr << "Error: parameterziation vf2.\n"; return; }
				mesh_.property(vuv_, *it) = mesh_.property(vuv_, vf0) * (mesh_.property(vbc_, *it))[0]
				+ mesh_.property(vuv_, vf1) * (mesh_.property(vbc_, *it))[1]
				+ mesh_.property(vuv_, vf2) * (mesh_.property(vbc_, *it))[2];
			}			
			(mesh_.property(fvset_, fh)).clear();//���������Ķ��������, ���������¿�����Щ�����������������.
		}
		for (unsigned int freev_i = 0, freev_size = free_vertices.size(); freev_i < freev_size; ++ freev_i) {
			mesh_.property(is_face_located, free_vertices[freev_i]) = false;//���е�free vertices��ʼΪ��û�ж�λ������.
		}


		for (unsigned int j = 1; j <= n-2; ++j) { //ѭ��n-2����, ��Щ������(from, to)->to֮��ȥ��fl��fr���γɵ�.
			TriMesh::VertexHandle v0 = loop[0];//=to
			TriMesh::VertexHandle vj = loop[j];
			TriMesh::VertexHandle vj1 = loop[j+1];//�����������γ�һ��������, �����ж�free_vertices��������Ķ����Ƿ���������.
			TriMesh::FaceHandle fh = mesh_.face_handle(loop_h[j]);

			for (std::vector<TriMesh::VertexHandle>::iterator freev_it = free_vertices.begin(), freev_end = free_vertices.end(); freev_it != freev_end; ++freev_it) {
				if (mesh_.property(is_face_located, *freev_it) == false) {

					OpenMesh::Vec2d pp = mesh_.property(vuv_, *freev_it);//���free vertex�Ĳ�������
					OpenMesh::Vec3d bc = DGP::calc_barycentric_coordinates(mesh_.property(vuv_,v0), mesh_.property(vuv_, vj), mesh_.property(vuv_, vj1), pp);
					//std::cout << bc << std::endl;//for test
					//(0.5, -1.11022e-016, 0.5)//�����������������������,���е�������С��С��-1.11022e-016��5.55112e-017��Ӧ����Ϊ0��,
					//(0.5, 0.5, 5.55112e-017) //�����Ļ���һ������������ʵ�ǺϺ���һ���������ڲ�(�߽���)���ڲ����Ҫ���.
					if (fabs(bc[0]) < 1.0e-10) bc[0] = 0;
					if (fabs(bc[1]) < 1.0e-10) bc[1] = 0;
					if (fabs(bc[2]) < 1.0e-10) bc[2] = 0;
					if (!(bc[0]<0) && !(bc[1]<0) && !(bc[2]<0)) { //Ҫ�����bc��������������>=0

						// For interior vertex, its' barycentric coordinates is bc according to the triangle(v0, vj, vj1).
						//std::cout << *freev_it << ": " << v0 << "," << vj << "," << vj1 << ": " << bc << "\n"; // for test
						// �������о���ÿһ�β�������Ӧ�û�õĻ��Ǹ��ĵ���Ϣ.
						mesh_.property(vf_, *freev_it) = fh;
						mesh_.property(vf0_, *freev_it) = v0; mesh_.property(vf1_, *freev_it) = vj; mesh_.property(vf2_, *freev_it) = vj1;
						mesh_.property(vbc_, *freev_it) = bc;

						(mesh_.property(fvset_, fh)).push_back(*freev_it);

						mesh_.property(is_face_located, *freev_it) = true;
					}
				} // ���free vertex�����������
			} // һ����free vertex �����������			
		}

		mesh_.property(vbc_, from) = bc_from;
		mesh_.property(vf0_, from) = to; mesh_.property(vf1_, from) = loop[n-2]; mesh_.property(vf2_, from) = loop[n-1];
		mesh_.property(vf_, from) = mesh_.face_handle(loop_h[n-2]);
		(mesh_.property(fvset_, mesh_.face_handle(loop_h[n-2]))).push_back(from);

		for (unsigned int freev_i = 0, freev_size = free_vertices.size(); freev_i < freev_size; ++ freev_i) { //���������жϵ�����,���ƿ���ȥ����.
			if (mesh_.property(is_face_located, free_vertices[freev_i]) == false) { //��ʾ���free vertex�������n-2������
				std::cout << freev_size << ", " << free_vertices[freev_i];
				std::cerr << "Error: parameterization case 3.\n";// û���ҵ�ͶӰ��������һ����, ���Ǵ����.
				return;
			}  
		}
		mesh_.remove_property(is_face_located);

	} else { std::cout << "Error: parameterize case4.\n"; 
	} // end of 3 cases
} // endl of function initial_parameterization /**/

void DecimationModel::decimate(unsigned int _n_vertices) //�û������򻯵����������Ŀ
{
	// build priority queue
	TriMesh::VertexIter  v_it  = mesh_.vertices_begin(), v_end = mesh_.vertices_end();
	queue.clear();
	vhierarchy_.clear();
	for (; v_it!=v_end; ++v_it) {//evaluate the cost of every vertex(pair)'s contraction
		// the minimum cost will have the highest priority, which in the front of the queue
		enqueue_vertex(v_it);//���� enqueue_vertex(v_it.handle());

		// ����һ��ɭ��, ���еĳ�ʼ�ڵ�(Ҳ����Ҷ�ӽڵ�)����ȫ��origianl mesh�е�ԭʼ����
		int new_node_handle = vhierarchy_.add_node();
		mesh_.property(vp_leaf_node_handle_, v_it) = new_node_handle;
		mesh_.property(vp_node_handle_, v_it.handle()) = new_node_handle;
		vhierarchy_.node(new_node_handle).set_vertex_handle(v_it.handle());//�ڵ��Ӧ�ĸ�����
		vhierarchy_.node(new_node_handle).set_active(true);

		assert(new_node_handle == v_it.handle().idx());// ���г�ʼ���㶼��һ��������idx()һ����ŵĽڵ�Node,��ЩNode��ΪҶ��
		//Ҳ����˵ÿһ��original mesh�ϵĶ��㶼��Ӧһ��Ҷ�ӽڵ�, ��Ҷ�ӽڵ����ž��Ƕ�������.
	}
	vhierarchy_.set_n_leaves(mesh_.n_vertices());//Ҷ�Ӷ��ٸ�
	// ��ʼ��ÿһ����߽ṹ, ʹ��ָ������opposite vertex handle.//Ŀ���Ƿ���halfedge collpaseʱ������ҵ�fundamental cut vertex.
	for (TriMesh::HIter he_it(mesh_.halfedges_begin()), he_end(mesh_.halfedges_end()); he_it != he_end; ++he_it) {
		if (mesh_.is_boundary(he_it) == false) {
			mesh_.property(heovh_, he_it) = mesh_.to_vertex_handle(mesh_.next_halfedge_handle(he_it.handle()));
		} else {
			mesh_.property(heovh_, he_it) = TriMesh::VertexHandle(-1);//TriMesh::VertexHandle InvalidVertexHandle;//idx() == -1
		}

		if (mesh_.status(mesh_.edge_handle(he_it)).feature()) { //������ʱ��crease vertex�����ⶨλ
			mesh_.property(hep_heset_, he_it).push_back(he_it.handle());//��ʼ��Ϊ�������Լ�
		}
	}

	unsigned int nv(mesh_.n_vertices());

	TriMesh::HalfedgeHandle hh;// the target halfedge to be collapsed
	TriMesh::VertexHandle   from, to;// the target halfedge's two vertex
	TriMesh::VVIter         vv_it;

	std::vector<TriMesh::VertexHandle>            one_ring;//from������ھ�
	std::vector<TriMesh::VertexHandle>::iterator  or_it, or_end;

	//----һֱ����, �ﵽĿ��ڵ���Ϊֹ. nv�������Ϊmesh_�л�����ô�������.
	while (nv > _n_vertices&& !queue.empty())
	{	// std::cout << nv << "\t";//for test  
		// Decimate using priority queue:
		//   1) take 1st element of queue
		//   2) collapse this halfedge
		//   3) update queue

		// get 1st element Ӧ���Ǵ�����С��,
		TriMesh::VertexHandle vh = *queue.begin();//enqueue_vertex(v_it);ʱ���ǽ�����ߵ����Ž������е�
		queue.erase(queue.begin());//ȡ������ʹ�queue�н���ɾȥ

		hh   = mesh_.property(vp_target_, vh);//��Ҫ�����İ��Ϊhh
		to   = mesh_.to_vertex_handle(hh);
		from = mesh_.from_vertex_handle(hh);//from Ӧ�ú�vh���, Ҳ����from = vh;

		// store one-ring vertex neighbor of the from vertex
		one_ring.clear();
		for (vv_it=TriMesh::VVIter(mesh_, from); vv_it; ++vv_it) {
			//for (vv_it(mesh_, from);            vv_it; ++vv_it) //������һ����
			//for (vv_it = mesh_.vv_iter(from);   vv_it; ++vv_it) //������һ����
			one_ring.push_back(vv_it.handle());
		}

		// perform collapse
		//if (mesh_.is_collapse_ok(hh))
		if (DGP::QEMDecimationT<TriMesh>::is_collapse_legal(mesh_, hh, vp_type_))// is_collapse_legal(hh),�Ѿ�������mesh_.is_collapse_ok(hh)�Ĳ���
		{	//֮���Ի���Ҫ��εļ������Ϊfrom�����һ���������������ܸı���,
			//�ı���ԭ�������һ������ĵ㷢��������.

			//postprocess_collapse(hh);//��仰������mesh_.collapse(hh);֮ǰ,��Ϊ��߱�����֮��topology�ͱ��ı���
			//std::cout << from.idx() << ", " << to.idx() << ": " << mesh_.property(heovh_, hh) << ", " << mesh_.property(heovh_, mesh_.opposite_halfedge_handle(hh)) << ".\n";//for test
			//------------
			// �����������֮����ָ���opposite vertex handle, ��������Ŀ����Ϊ�˸�����ҵ�fundamental cut vertices
			if (mesh_.is_boundary(hh) == false) { // ͼֽ����˵��, ����openmesh��ʹ�õ�halfedge collapse���ص�.
				TriMesh::HHandle n = mesh_.next_halfedge_handle(hh);
				TriMesh::HHandle po = mesh_.opposite_halfedge_handle(mesh_.prev_halfedge_handle(hh));
				mesh_.property(heovh_, n) = mesh_.property(heovh_, po);
			}
			if (mesh_.is_boundary(mesh_.opposite_halfedge_handle(hh)) == false) {
				TriMesh::HHandle o = mesh_.opposite_halfedge_handle(hh);
				TriMesh::HHandle op = mesh_.prev_halfedge_handle(o);
				TriMesh::HHandle ono = mesh_.opposite_halfedge_handle(mesh_.next_halfedge_handle(o));
				mesh_.property(heovh_, op) = mesh_.property(heovh_, ono);
			}
			// �ҳ���Ҫ���������ߵ�����fundamental cut vertices, �����ҵ�fundamental cut vertices����Ӧ��Node�ڵ�
			TriMesh::VHandle fund_cut_vertex_handle0 = mesh_.property(heovh_, hh);//��Щfundamental cut vertex handle�ܿ�����invaild��,��idx() == -1�պ�Ҳ��ʾû��Node�����Ӧ.
			TriMesh::VHandle fund_cut_vertex_handle1 = mesh_.property(heovh_, mesh_.opposite_halfedge_handle(hh));
			int fund_cut_node_handle0 = fund_cut_vertex_handle0.idx();//�������صػ���һ����ʵ: original mesh�ϵ�vertexhandle��idx�պ������Ӧ�Ľڵ�Node(Ҷ��)�����һ��.
			int fund_cut_node_handle1 = fund_cut_vertex_handle1.idx();//����fundamental cut vertex����original mesh�ϵĶ���.���ֵҲ������-1.
			// -------------
			if (do_parameterization_) initial_parameterization(hh);

			//----- to collapse the halfedge
			mesh_.collapse(hh);//halfedge collapse, Ĭ�ϵľ���(from, to)->to
			simplified_mesh_.collapse(hh);
			// (1) ���¶���to��quadric, ����ԭQEM������:
			mesh_.property(vp_quadric_, to) += mesh_.property(vp_quadric_, from); //(from, to)->to, ����û����Vertex Placement Policies// 
			// (2) ���¶���to��һ����ıߺͶ��������, �ⶫ������Ч����������hemisphere.offģ����.
			bool update_type = false;
			if (update_type) {
				int n_fe = 0;
				for (TriMesh::VEIter ve_it(mesh_, to); ve_it; ++ve_it) {
					//if (!mesh_.status(ve_it).deleted()) { // ����ж��������Ƕ�����.
					if (fabs(mesh_.calc_dihedral_angle(ve_it)) >  M_PI*(44.0/180)) { 
						mesh_.status(ve_it).set_feature(true); 
						n_fe ++;
					} else {
						mesh_.status(ve_it).set_feature(false);
					}
					//}
				}
				if (n_fe == 0) { mesh_.property(vp_type_, to) = DGP::SMOOTH_VFT;  } else if (n_fe == 1) { mesh_.property(vp_type_, to) = DGP::DART_VFT; } 
				else if (n_fe == 2) { mesh_.property(vp_type_, to) = DGP::CREASE_VFT;  } else if (n_fe >= 3) { mesh_.property(vp_type_, to) = DGP::CORNER_VFT; }

				for (TriMesh::VVIter vv_it(mesh_, to); vv_it; ++vv_it) {
					n_fe = 0; // ����������ܵ�feature edges�ĸ���
					for (TriMesh::VEIter ve_it = mesh_.ve_iter(vv_it.handle()); ve_it; ++ve_it) {
						if (mesh_.status(ve_it).feature()) n_fe++;
					}
					if (n_fe == 0) { mesh_.property(vp_type_, vv_it) = DGP::SMOOTH_VFT; } else if (n_fe == 1) { mesh_.property(vp_type_, vv_it) = DGP::DART_VFT; } 
					else if (n_fe == 2) { mesh_.property(vp_type_, vv_it) = DGP::CREASE_VFT; } else if (n_fe >= 3) { mesh_.property(vp_type_, vv_it) = DGP::CORNER_VFT; } 
				}
			}
			// (3) ����vertex hierarchy
			int new_node_handle = vhierarchy_.add_node();//ÿ��halfedge collapse֮������һ���µĽڵ�.
			vhierarchy_.node(mesh_.property(vp_node_handle_, from) ).set_parent_node_handle(new_node_handle);
			vhierarchy_.node(mesh_.property(vp_node_handle_, to) ).set_parent_node_handle(new_node_handle);
			vhierarchy_.node(mesh_.property(vp_node_handle_, from) ).set_active(false);
			vhierarchy_.node(mesh_.property(vp_node_handle_, to) ).set_active(false);

			vhierarchy_.node(new_node_handle).set_lchild_node_handle(mesh_.property(vp_node_handle_, from) );
			vhierarchy_.node(new_node_handle).set_rchild_node_handle(mesh_.property(vp_node_handle_, to) );
			vhierarchy_.node(new_node_handle).set_vertex_handle(to);
			vhierarchy_.node(new_node_handle).set_active(true);

			mesh_.property(vp_node_handle_, to) = new_node_handle;//end. to���������ָ���Node.
			assert(vhierarchy_.node(new_node_handle).vertex_handle().idx() == to.idx());

			vhierarchy_.node(new_node_handle).set_fund_cut_node_handle0(fund_cut_node_handle0);
			vhierarchy_.node(new_node_handle).set_fund_cut_node_handle1(fund_cut_node_handle1);

			//std::cout << "fidx2 " << mesh_.status(from).deleted() << ", " << from.is_valid()<< std::endl; //for test
			//���������� fidex2 1, 1, ����ζ�����from�����Ǳ�־Ϊdeleted��, ������������Ч���ڴ洢����������,
			//���Ի�����Ч��, 
			//��Ϊֻ����idx_(���Դ�idx()�������)��������-1,�Ǿ��������ϻ���Ч,����������ʹ��mesh_.point(from)���������.
			if (vertex_optimal_placement == true) {
				double newvp0 = 0.0, newvp1 = 0.0, newvp2 = 0.0;
				if (mesh_.property(vp_quadric_, to).optimal_placement(newvp0, newvp1, newvp2)) {//
					//std::cout << newvp0 << " " <<  newvp1 << " " << newvp2 << "\n"; //for test
					mesh_.set_point(to, TriMesh::Point((float)newvp0, (float)newvp1, (float)newvp2));
				}
			}
			--nv;
			if (nv % 100 == 0) {	
				std::cerr << nv << "\r";//"\n";//
				if (do_parameterization_)  //����һ��ʮ��crazy������, ÿ��һЩ�����ƽ��.
				{	//������mesh_.garbage_collection()֮ǰ,�����޷�ʹ����ǰ��Handle.
					//local_smooth_parameterization_triangledomain();
					//local_smooth_parameterization_circledomain();					
				}
			}

		} else { // end of if(is_collapse_legal(hh))
			std::cout << "Here?\n";
		}

		// update queue, һ�����㱻����֮������one-ring�еĶ���Ҫ�������ǵĴ��ۺ��������Ϊ���Ű��.
		for (or_it=one_ring.begin(), or_end=one_ring.end(); or_it!=or_end; ++or_it)
			enqueue_vertex(*or_it);
	} // end of while(nv > _n_vertices).

	//-------------------����Ѿ����򻯹���mesh_�Ĵ洢�ṹ����ô����
	/*	v_it= mesh_.vertices_begin();
	v_end= mesh_.vertices_end();
	int si(0); 
	for (; v_it != v_end; ++v_it) { //������еĵ�
	if (mesh_.status(v_it).deleted() == false)  {
	//��������, ����ֻ������򻯶��ɵ�base model�еĶ���,
	//����idx()��ʾ��index��������original model��һ��, ���ǳ�ʼ����.off�ļ�ʱ�򶥵������,
	//��Ҫ�ȵ�mesh_.garbage_collection();������֮���������Ż�ı�,��ʱ֮��֮ǰ��handle����Ч��.
	//���Ƿ�˵���˱�ɾ���ĵ�ֻ�Ǳ���־Ϊdeleted(ͨ��mesh_.status(v_it.handle()).set_deleted(true)���), 
	//����ʵ�ⱻɾ���ĵ㻹���ڴ洢�ṹ����, ������Ч��.
	std::cout << "���base model " << v_it.handle().idx() << ": " << mesh_.point(v_it.handle()) << std::endl;
	++si;
	}
	//��������ļ�����Կ�����������дӳ�ʼ����.off�ļ�ʱ��d����Ķ���,������idx()���������
	//�����ж����is_valid()������1��ʾ�ڴ洢�ṹ�ϻ��Ǵ���������Ч, ֻҪidx() != -1��is_valid() == true
	//������deleted()���ڱ�ɾ���Ķ��㷵�ص���1,��Ӧ�򻯺���base model�еĵ�(��λ�ÿ��ܱ����¼���)��0.
	//std::cout << v_it.handle().idx() << ", " 
	//		  << v_it.handle().is_valid() << ", " << mesh_.status(v_it.handle()).deleted() << ", " 
	//		  << mesh_.point(v_it.handle()) << std::endl;
	}	
	std::cerr << "û�б�ɾ���Ķ������: " << si << "\n";*/
	//------------------------------------------------------------
	//write(std::string("test.pm"));

	// ------------------
	// clean up
	queue.clear();
	// �����, vertex hierarchyҲ�������, ��ṹ�͹̶���.
	/*// -----------------------
	for (v_it = mesh_.vertices_begin(), v_end = mesh_.vertices_end(); v_it != v_end; ++v_it) {
	if (mesh_.status(v_it).deleted() == false) { // base mesh
	int node_handle = mesh_.property(vp_node_handle_, v_it);
	//std::cout << v_it.handle() << ", " << mesh_.point(v_it) << ", " << node_handle << "; " << "\n";//for test.
	}
	}*/

	//--------------------------------------------------------------
	// ��mesh_.garbage_collection()ǰ��������, �ο�ModProgMeshT.cc��write����������,������������ӳ��.
	int N = mesh_.n_vertices();
	vpoints_ = std::vector<TriMesh::Point>(N);//���ж������������, �ȴ��base mesh�еĶ���, �ٴ�ű�deleted�Ķ���
	voldhandles_ = std::vector<TriMesh::VHandle> (N);
	int i = 0;//��Ϊ����vpoint_���±�ͼ���
	int n_base_vertices = 0, n_detail_vertices = 0;
	OpenMesh::VPropHandleT<int> idx_vps;//mesh_ÿһ�������Ӧ����vpoints_����һ��Ԫ����? ������idx��Ϊ�±���.
	mesh_.add_property(idx_vps);
	for (v_it = mesh_.vertices_begin(), v_end = mesh_.vertices_end(); v_it != v_end; ++ v_it) {
		if (mesh_.status(v_it).deleted() == false) { 
			vpoints_[i] = mesh_.point(v_it.handle()); 	//�����.
			voldhandles_[i] = v_it.handle();//because: mesh_.point(voldhandles_[i]) is ok.
			mesh_.property(idx_vps, v_it.handle()) = i;//mesh_.property(idx_vps, v_it.handle) = vhandles.size() - 1;
			i ++;

			simplified_mesh_.property(vp_node_handle_sm_, v_it) = mesh_.property(vp_node_handle_, v_it); 
		}
	}
	n_base_vertices = i; //std::cout << "n_base_vertices: " << n_base_vertices << std::endl;// for test//
	for (v_it = mesh_.vertices_begin(), v_end = mesh_.vertices_end(); v_it != v_end; ++ v_it) {
		if (mesh_.status(v_it).deleted() == true) {
			vpoints_[i] = mesh_.point(v_it.handle());	
			voldhandles_[i] = v_it.handle();
			mesh_.property(idx_vps, v_it.handle()) = i;//mesh_.property(idx, v_it.handle) = vhandles.size() - 1;
			i ++;
		}
	}
	n_detail_vertices = N - n_base_vertices;	assert(i == N);

	// ��vertex hierarchy��ÿһ��Node��int point_handle_��������Ϊvpoints_�����еĶ�Ӧ�±�.
	// �Ӷ���vertex hierarchy�еĽڵ�����ڱ���ȫ�����������vpoints_���齨�����ϵ.
	for (int k = 0, vhierarchy_size = vhierarchy_.size(); k < vhierarchy_size; ++k) {
		TriMesh::VertexHandle vh = vhierarchy_.node(k).vertex_handle();
		int point_handle = mesh_.property(idx_vps, vh);//���һ��ʹ����VertexHandle��,��mesh_.garbage_collection();֮�����Ч��.
		vhierarchy_.node(k).set_point_handle(point_handle);//�˺�vhierarchy��ͨ����point_handle����mesh_�ϵĶ�Ӧ����.
		
		vhierarchy_.node(k).set_vertex_handle(TriMesh::VertexHandle(-1));//���ÿһ��Node��vertex_handle��.
	}
	mesh_.remove_property(idx_vps);

	for (TriMesh::FIter f_it(simplified_mesh_.faces_begin()), f_end(simplified_mesh_.faces_end()); f_it != f_end; ++f_it) {
		if (simplified_mesh_.status(f_it.handle()).deleted() == false) {
			simplified_mesh_.property(fp_fh_sm_, f_it) = f_it.handle(); //û�б�ɾ�������¼��garbage_collection()�Ķ�Ӧ��.
		} 
	}
	// -----------------------------------------------------------------------
	//std::cout << mesh_.n_vertices() << ", " << mesh_.n_edges() << ", " << mesh_.n_faces() << "\n";//��������ʾ�㱻����ǰ�����,��ʱ��û��garbage_collection(),
	// now, delete the items marked to be deleted
	//mesh_.garbage_collection();//����Ǳ����,����ģ���б�set_deleted(true)�ĵ㻹����ģ����,���ǲ��Ե�.
	// In the simplification process, the mesh_ was changed, so update the normals of vertices and faces.
	//mesh_.update_normals();//MeshModel::update();
	//std::cout << mesh_.n_vertices() << ", " << mesh_.n_edges() << ", " << mesh_.n_faces() << "\n";
	simplified_mesh_.garbage_collection();//����mesh_���Ǳ���ԭ��������,Ϊ�˺��滹��Ҫ�õ�����������Ϣ.
	simplified_mesh_.update_normals();		// NOTE: ��Ϊ��Ϊmesh_û�и���, ��������Ĵ���Ҫ�ĳ���simplified_mesh_��������.


	/// ----------
	// ���ڼ�֮���ģ��, ����ÿһ����������Ӧ��Node, ����Ӧ��Node����������ָ���vertex handle.
	for (TriMesh::VIter v_it(simplified_mesh_.vertices_begin()), v_end(simplified_mesh_.vertices_end()); v_it != v_end; ++v_it) {
		int node_handle = simplified_mesh_.property(vp_node_handle_sm_, v_it); 
		vhierarchy_.node(node_handle).set_vertex_handle(v_it.handle());
	}

	for (TriMesh::FIter f_it(simplified_mesh_.faces_begin()), f_end(simplified_mesh_.faces_end()); f_it != f_end; ++f_it) {
		mesh_.property(fp_fh_, simplified_mesh_.property(fp_fh_sm_, f_it)) = f_it.handle(); //mesh_��ÿһ��(û�б�ɾ����)�涼֪������simplified_mesh_�Ķ�Ӧ��.
	}

	// ----------- feature and print out 
	int fe = 0;
	int sv = 0, dv = 0, crv = 0, cov = 0, bv = 0;
	/*
	for (TriMesh::EIter e_it(mesh_.edges_begin()), e_end(mesh_.edges_end()); e_it != e_end; ++e_it) {
	if (mesh_.status(e_it).feature()) ++fe;
	}
	for (TriMesh::VIter v_it = mesh_.vertices_begin(), v_end = mesh_.vertices_end(); v_it != v_end; ++ v_it) {
	//	std::cout << v_it.handle() << ", " << mesh_.property(vp_node_handle_, v_it) << ", " << vhierarchy_.node(mesh_.vertex(v_it.handle()).node_handle()).vertex_handle()<< "\n";
	if (mesh_.property(vp_type_, v_it) == DGP::SMOOTH_VFT) sv++;
	if (mesh_.property(vp_type_, v_it) == DGP::DART_VFT) dv++;
	if (mesh_.property(vp_type_, v_it) == DGP::CREASE_VFT) crv++;
	if (mesh_.property(vp_type_, v_it) == DGP::CORNER_VFT) cov++;
	if (mesh_.is_boundary(v_it) == true) bv++;
	//std::cout << mesh_.property(vp_type_, v_it) << ", ";
	}*/
	for (TriMesh::EIter e_it(simplified_mesh_.edges_begin()), e_end(simplified_mesh_.edges_end()); e_it != e_end; ++e_it) {
		if (simplified_mesh_.status(e_it).feature()) ++fe;
	}
	for (TriMesh::VIter v_it = simplified_mesh_.vertices_begin(), v_end = simplified_mesh_.vertices_end(); v_it != v_end; ++ v_it) {
		if (simplified_mesh_.property(vp_type_sm_, v_it) == DGP::SMOOTH_VFT) sv++;
		if (simplified_mesh_.property(vp_type_sm_, v_it) == DGP::DART_VFT) dv++;
		if (simplified_mesh_.property(vp_type_sm_, v_it) == DGP::CREASE_VFT) crv++;
		if (simplified_mesh_.property(vp_type_sm_, v_it) == DGP::CORNER_VFT) cov++;
		if (simplified_mesh_.is_boundary(v_it) == true) bv++; 
	}
	std::cout << "�򻯻�õĻ�������������: " << fe << ".\n";
	std::cout << "�򻯻�õĻ�����Ķ�����: smooth: " << sv << ", dart " << dv << ", crease: " << crv << ", corner: " << cov << ", boundary: " << bv << ".\n";

	/*
	std::cout << "\n";// for test
	for (int k = 0, vhsize = vhierarchy_.size(); k < vhsize; ++k) {
	DGP::VHierarchyNode n = vhierarchy_.node(k);
	std::cout << n.self_node_handle() << ", " << n.is_active() << ", " << n.point_handle() << ", " << n.vertex_handle() << "; " 
	<< n.parent_node_handle() << ", " << n.lchild_node_handle() << ", " << n.rchild_node_handle() << "; " 
	<< n.fund_cut_node_handle0() << ", " << n.fund_cut_node_handle1() << "\n";
	}*/	
	// ������Ϊѹƽ�Ĳ�����֮��Ч������, ����smooth��.
	//local_smooth_parameterization_process();
}

void DecimationModel::copy_backup_mesh() { // ����һ��ԭʼ����Ŀ���, ������Ҫȥ�޸���.
	mesh2_ = mesh_; TriMesh::ConstFaceIter        f_it(mesh2_.faces_sbegin()), f_end(mesh2_.faces_end());
	TriMesh::ConstFaceVertexIter  cfv_it;

	mesh2_indices_for_render_.clear();
	mesh2_indices_for_render_.reserve(mesh2_.n_faces()*3);

	for (; f_it!=f_end; ++f_it)
	{	// �������е���, ��¼��ÿһ����Ķ������
		mesh2_indices_for_render_.push_back((cfv_it=mesh2_.cfv_iter(f_it)).handle().idx());
		mesh2_indices_for_render_.push_back((++cfv_it).handle().idx());
		mesh2_indices_for_render_.push_back((++cfv_it).handle().idx());
	}

	calc_backupmesh_vertex_area();
	calc_backupmesh_halfedge_meanvaluecoor();
}
void DecimationModel::calc_backupmesh_vertex_area() {
	mesh2_.add_property(varea_);
	TriMesh::VertexIter        v_it, v_end(mesh2_.vertices_end());
	TriMesh::VertexFaceIter    vf_it;
	TriMesh::FaceVertexIter    fv_it;
	TriMesh::Scalar area;
	// compute the area associated to a vertex. ����һ�������������, ��һ������ֵ
	// The simplest such area is the barycentric area, 
	// given by 1/3 the area of all incident triangles. ����ʹ�õľ����ⷽ��.
	// A better measure of the area is the Voronoi area associated to the vertex.
	for (v_it=mesh2_.vertices_begin(); v_it!=v_end; ++v_it)
	{
		area = 0.0;

		for (vf_it=mesh2_.vf_iter(v_it); vf_it; ++vf_it)
		{
			fv_it = mesh2_.fv_iter(vf_it);

			const TriMesh::Point& P = mesh2_.point(fv_it);  ++fv_it;
			const TriMesh::Point& Q = mesh2_.point(fv_it);  ++fv_it;
			const TriMesh::Point& R = mesh2_.point(fv_it);

			area += ((Q-P)%(R-P)).norm() * 0.5f * 0.3333f;//���������������֮һ
		}

		mesh2_.property(varea_, v_it) = (fabs(area)>FLT_MIN ? area : 0.0);
		//std::cout << mesh2_.property(varea_, v_it) << " ";
	}

}
void DecimationModel::calc_backupmesh_halfedge_meanvaluecoor() {
	mesh2_.add_property(hmvc_);
	for (TriMesh::HalfedgeIter h_it(mesh2_.halfedges_begin()), h_end(mesh2_.halfedges_end()); h_it != h_end; ++h_it) 
	{
		if (mesh2_.is_boundary(h_it) == false && mesh2_.is_boundary(mesh2_.opposite_halfedge_handle(h_it)) == false) {
			//������h_it.handle()����Ӧ��ϵ��wij, ע������ʹ�õ���mean value coordiante, ������harmonic map�е�cotangent weight
			TriMesh::VertexHandle v0 = mesh2_.from_vertex_handle(h_it);
			TriMesh::VertexHandle v1 = mesh2_.to_vertex_handle(h_it);
			TriMesh::VertexHandle v2 = mesh2_.to_vertex_handle(mesh2_.next_halfedge_handle(h_it.handle()));// ������������ȷ��ǰ���ǷǱ߽����
			TriMesh::VertexHandle v3 = mesh2_.to_vertex_handle(mesh2_.next_halfedge_handle(mesh2_.opposite_halfedge_handle(h_it.handle())));
			OpenMesh::Vec3d v1v0 = mesh2_.point(v1) - mesh2_.point(v0);	// ֮ǰ������ֵ�-1
			OpenMesh::Vec3d v2v0 = mesh2_.point(v2) - mesh2_.point(v0);
			OpenMesh::Vec3d v3v0 = mesh2_.point(v3) - mesh2_.point(v0);
			double v1v0v2_angle = acos(dot(v1v0, v2v0) / (v1v0.norm() * v2v0.norm()));// ʸ��v1v0��v2v0�ļн�.
			double v3v0v1_angle = acos(dot(v3v0, v1v0) / (v3v0.norm() * v1v0.norm()));// ʸ��v3v0��v1v0�ļн�.
			double wij = (tan(v1v0v2_angle / 2) + tan(v3v0v1_angle / 2)) / v1v0.norm();
			if (wij < 0) { std::cout << "Error: calc_backupmesh_halfedge_meanvaluecoor(): wij negative.\n"; return; }//mean value coordinateӦ���Ǵ������.
			else mesh2_.property(hmvc_, h_it.handle()) = wij;
		}		
	}
} // end of function calc_backupmesh_halfedge_meanvaluecoor().
void DecimationModel::decimation_process(unsigned int _num, bool _feature_or_not) {
	std::cout << "\n";
	std::cout << "========��ʼ���򻯹���==================\n";//������ʾ��40�����ȵ�=,ǰ���9��λ�ÿ�ʼ����,��������8��=.
	copy_backup_mesh();
	add_required_properties();
	simplified_mesh_ = mesh_;  //simplified_mesh_��ʼ��mesh_һ��, ������mesh_�ļ򻯶���, ֻ��simplified_mesh_û��garbage_colection().
	simplified_mesh_.add_property(vp_type_sm_); 
	simplified_mesh_.add_property(vp_node_handle_sm_);
	 
	//DGP::identify_all_sharp_features(mesh_, vp_type_, false, 30);	 
	feature_considered_ = _feature_or_not;
	double const gfxDEGTORAD =  0.01745329251994329577;//�Ƕȵ�λ��ת�� pi/180 degree to radian
	unsigned int n_feature_edges = 0;
	if (_feature_or_not ) {
		n_feature_edges = mesh_.find_feature_edges(44 * gfxDEGTORAD); 
	}
	//std::cout << filename_ << ", " << fandisk.size() << ".\n";
	//if (filename_.substr(filename_.rfind('/') == std::string::npos ? filename_.length() : filename_.rfind('/') + 1) 
	//	== std::string("fan_t_6475.off")) { 
	//	for (TriMesh::EIter e_it(mesh_.edges_begin()), e_end(mesh_.edges_end()); e_it != e_end; ++e_it) {
	//		if (std::find(fandisk.begin(), fandisk.end(), e_it.handle().idx()) != fandisk.end()) {
	//			mesh_.status(e_it).set_feature(true);
	//			++n_feature_edges;
	//		}
	//	}
	//}
	std::cout << "Num of feature edges(not include boundaries) " << n_feature_edges << std::endl; //"\t";
	
	int s = 0, d = 0, cr = 0, co = 0, b = 0; // for test 
	for (TriMesh::VIter v_it = mesh_.vertices_begin(), v_end = mesh_.vertices_end(); 
		v_it != v_end; ++v_it) {
		if (mesh_.status(v_it).deleted() == false) { // ���������, ��ʱ���ǵ��������ԭ��.

			if (mesh_.is_boundary(v_it)) {
				b++;
			} 
			unsigned int n_fe = 0; // ����������ܵ�feature edges�ĸ���
			for (TriMesh::VEIter ve_it = mesh_.ve_iter(v_it.handle()); ve_it; ++ve_it) {
				if (mesh_.status(ve_it).feature()) n_fe++;//
				//if (_m.status(ve_it).feature() || _m.is_boundary(ve_it)) n_fe++;//�߽��Ҳ����������.
				//Ϊ�˳���ļ���, boundary edges and vertices ��������.
			}
			if (n_fe == 0) {
				mesh_.property(vp_type_, v_it) = DGP::SMOOTH_VFT;	s++;	
			} else if (n_fe == 1) {
				mesh_.property(vp_type_, v_it) = DGP::DART_VFT;	d++;	
			} else if (n_fe == 2) {
				mesh_.property(vp_type_, v_it) = DGP::CREASE_VFT;	cr++;
			} else if (n_fe >= 3) {
				mesh_.property(vp_type_, v_it) = DGP::CORNER_VFT;	co++;	
			}
		}
	} 
	std::cout << "Num of feature vertex: SMOOTH_VFT " << s << ", DART_VFT " << d 
		<< ", CREASE_VFT " << cr << ", CORNER_VFT " << co << "; ����BOUNDARY " << b << ". \n"; 

	for (TriMesh::VIter v_it = mesh_.vertices_begin(), v_end(mesh_.vertices_end()); v_it != v_end; ++v_it) {
		if (mesh_.status(v_it).deleted() == false) { // ���������
			simplified_mesh_.property(vp_type_sm_, v_it.handle()) = mesh_.property(vp_type_, v_it);
		}
	}
	for (TriMesh::EdgeIter e_it(mesh_.edges_begin()), e_end(mesh_.edges_end()); e_it != e_end; ++e_it) {
		simplified_mesh_.status(e_it).set_feature(mesh_.status(e_it).feature());//copy the mesh_'s feature properties
	}////

	mesh_.add_property(fp_fh_);
	simplified_mesh_.add_property(fp_fh_sm_);

	for (TriMesh::VIter v_it(mesh_.vertices_begin()), v_end(mesh_.vertices_end()); v_it != v_end; ++v_it) 
		mesh_.property(vp_cost_, v_it) = -1.0;
	DGP::QEMDecimationT<TriMesh> qem_decimation;
	qem_decimation.init_calc_quadric(mesh_, vp_quadric_, false);


	do_parameterization_ = true;//false;// //�򻯵�ͬʱ�Ƿ���в�����
	if (do_parameterization_) { mesh_collapsed_ = mesh_; mesh_collapsed_render_or_not_ = false; }
	do_smooth_ = true;//false;//initial parameterization֮���Ƿ����smooth
	n_local_sp_ = 2;//initial parameterization֮�������ô��ε�local smooth parameterization 
	std::cout << "--------����򻯹���--------------------\n";
	std::cout << "DecimationModel: decimate to " << _num << ", begin.\n";
	decimate(_num);
	set_mesh_toberendered(&simplified_mesh_, SIMPLIFIED);
	remove_required_properties();	
	std::cout << "DecimationModel: end. " << simplified_mesh_.n_vertices() << " vertices, " << simplified_mesh_.n_edges() << " edges, " << simplified_mesh_.n_faces()    << " faces\n";
	std::cout << "--------�뿪�򻯹���--------------------\n";

	if (do_parameterization_) {
		// ------------������ʾ�������ĵ�, �����ǽ���Щ��������ļ�����PointViewer.exe�����鿴.
		//std::ofstream file("out.txt", std::ios_base::out);
		//file.clear();
		size_t N = mesh_.n_vertices(); 
		size_t i = 0;
		for (TriMesh::VertexIter v_it(mesh_.vertices_begin()), v_end(mesh_.vertices_end()); v_it != v_end; ++v_it, i++) {
			if (mesh_.status(v_it.handle()).deleted()) { 

				OpenMesh::Vec3d bc = mesh_.property(vbc_, v_it.handle());
				if (!((bc[0] > -1.0e-10) && (bc[1] > -1.0e-10) && (bc[2] > -1.0e-10))) 
				{ std::cerr << "Error: para check, bc negative." << bc << "\n"; } // for test, only assert

				TriMesh::Point vfp0 = mesh_.point(mesh_.property(vf0_, v_it.handle())); 
				if (mesh_.status(mesh_.property(vf0_, v_it.handle())).deleted() == true) { std::cerr << "Error: para check, vf0.\n"; }
				TriMesh::Point vfp1 = mesh_.point(mesh_.property(vf1_, v_it.handle()));
				if (mesh_.status(mesh_.property(vf1_, v_it.handle())).deleted() == true) { std::cerr << "Error: para check, vf1.\n"; }
				TriMesh::Point vfp2 = mesh_.point(mesh_.property(vf2_, v_it.handle()));
				if (mesh_.status(mesh_.property(vf2_, v_it.handle())).deleted() == true) { std::cerr << "Error: para check, vf2.\n"; }

				TriMesh::Point vp = vfp0 * bc[0] + vfp1 * bc[1] + vfp2 * bc[2];//��deletedҲ���Ǳ�parameterized�Ķ����3d����.

				//file << vp << "\n";//������obj��ʽ�Ļ�����<< "v " // �ļ������Щ���������ĵ�����
				//vpoints_for_render_[i] = vp;//��ɾȥ�Ļ��ͻᱻ���������������ĳһ����, ��������ȡ���������

				mesh_collapsed_.set_point(v_it, vp);////
			} else { 
			} 
		}
		//file.close();
	}

	// �������ϵı�, ����ʼ��Ϊ��û���е������.
	mesh_.add_property(empl_);//�ߵ��е��Ƿ��Ѿ�resample��.
	mesh_.add_property(emp_);//���resample,��3D����.
	mesh_.add_property(emp_closest_vh_);
	for (TriMesh::EdgeIter e_it(mesh_.edges_begin()), e_end(mesh_.edges_end()); e_it != e_end; ++ e_it) {
		if (mesh_.status(e_it).deleted() == false) { 
			mesh_.property(empl_, e_it.handle()) = false; // ��ʼ��, ����ߵ��е㻹û�ж�λ.
		}
	}

	std::cout << "========�����򻯹���====================\n";
}

// halfedge collapse������, ��base mesh�ָ���original mesh,
// Ҳ���Ǵ�vertex hierarchy�����ϲ㰴��ʱ����������η��ѵ����²�.
void DecimationModel::sequence_refine()  {
	std::cout << "\n";
	std::cout << "========����Sequence Refinement=========\n";//������ʾ��40�����ȵ�=,ǰ���9��λ�ÿ�ʼ����,��������8��=.
	// ѡ��Ҫsplit�Ķ���
	for (int k = vhierarchy_.size() - 1; k >= vhierarchy_.n_leaves(); --k) {// ecol ���������

		TriMesh::VertexHandle v1 = vhierarchy_.node(k).vertex_handle(); 

		int v1_node_handle = simplified_mesh_.property(vp_node_handle_sm_, v1);
		assert(k == v1_node_handle);
		//std::cout << v1.idx() << ", " << v1_node_handle << "\n";// for test, ע��, ����ʱ��ע����v1_node_handle ������-1��, 
		//-1��ʾ�������vertex handle��û��һ��������Ӧ��Node����.

		// �����������Ӧ��Node����Ҫactive��
		if (vhierarchy_.node(v1_node_handle).is_active()) { 
			int lc_node_handle = vhierarchy_.node(v1_node_handle).lchild_node_handle();
			int rc_node_handle = vhierarchy_.node(v1_node_handle).rchild_node_handle();
			if (lc_node_handle != -1) { // rc_node_handle != -1 //����ڵ㲻��Ҷ�ӽڵ�,��������split
				//std::cout <<"can be splitted: " << lc_node_handle << ", " << rc_node_handle << "\n";//for test
				//vsplit(to)->(from, to), Ҳ����vspit(v1)->(v0, v1).
				TriMesh::Point v_from_point = vpoints_[vhierarchy_.node(lc_node_handle).point_handle()];
				TriMesh::VHandle v_from = simplified_mesh_.add_vertex(v_from_point);//���������from��
				TriMesh::VHandle v_to = v1;		//std::cout << vpoints_[vhierarchy_.node(lc_node_handle).point_handle()] << "\n";//for test
				TriMesh::VHandle v_l = TriMesh::VertexHandle(-1);//��ʾinvalid vertex handle
				TriMesh::VHandle v_r = TriMesh::VertexHandle(-1);
				int fc0 = vhierarchy_.node(v1_node_handle).fund_cut_node_handle0();//std::cout << "fc0 n h: " << fc0 << "\n";//for test, fundamental cut vertex vl^.
				if (fc0 != -1) {
					DGP::VHierarchyNode fund_cut_node = vhierarchy_.node(fc0);
					while (fund_cut_node.is_active() == false) {
						fund_cut_node = vhierarchy_.node(fund_cut_node.parent_node_handle());
					}
					v_l = fund_cut_node.vertex_handle();
				}
				int fc1 = vhierarchy_.node(v1_node_handle).fund_cut_node_handle1();//std::cout << "fc1 n h: " << fc1 << "\n";//for test, fundamental cut vertex vr^.
				if (fc1 != -1) {
					DGP::VHierarchyNode fund_cut_node = vhierarchy_.node(fc1); //std::cout << fund_cut_node.self_node_handle() << "\n";
					while (fund_cut_node.is_active() == false) {
						fund_cut_node = vhierarchy_.node(fund_cut_node.parent_node_handle());
					}  //std::cout << fund_cut_node.self_node_handle() << ", " << fund_cut_node.vertex_handle() << "\n";//for test
					v_r = fund_cut_node.vertex_handle();
				}
				//std::cout << "vsplit " << v_from << ", " << v_to << ", " << v_l << ", " << v_r << ".\n";//for test
				simplified_mesh_.vertex_split(v_from, v_to, v_l, v_r);

				vhierarchy_.node(v1_node_handle).set_active(false);
				vhierarchy_.node(lc_node_handle).set_vertex_handle(v_from); 
				vhierarchy_.node(lc_node_handle).set_active(true);		simplified_mesh_.property(vp_node_handle_sm_, v_from) = lc_node_handle;
				vhierarchy_.node(rc_node_handle).set_vertex_handle(v_to);	
				vhierarchy_.node(rc_node_handle).set_active(true);		simplified_mesh_.property(vp_node_handle_sm_, v_to) = rc_node_handle;

			} // end of if. this node has lc and rc.
		} // end of if. this node is actived.
	}
	set_mesh_toberendered(&simplified_mesh_, ORIGINAL);
	simplified_mesh_.update_normals();//MeshModel::update();
	std::cout << "Now, the vertices num: " << simplified_mesh_.n_vertices() << ".\n";
	std::cout << "========����Sequence Refinement=========\n";//������ʾ��40�����ȵ�=,ǰ���9��λ�ÿ�ʼ����,��������8��=.
}
void DecimationModel::selective_refine() {
	std::cout << "\n";
	std::cout << "========����vsplit����==================\n";//������ʾ��40�����ȵ�=,ǰ���9��λ�ÿ�ʼ����,��������8��=.

	for (TriMesh::VIter v_it = simplified_mesh_.vertices_begin(), v_end = simplified_mesh_.vertices_end(); v_it != v_end; ++v_it) { //ÿ�δ�ɭ�ֵ���Ч�Ķ������·���һ��

		TriMesh::VertexHandle v1 = v_it.handle();// v1�����Ӧ��v1_node_handle����ڵ�Node������.
		int v1_node_handle = simplified_mesh_.property(vp_node_handle_sm_, v1);
		std::cout << v1.idx() << ", " << v1_node_handle << "\n";// for test, ע��, ����ʱ��ע����v1_node_handle ������-1��, -1��ʾ�������vertex handle��ÿһ��������Ӧ��Node����.
		// �����������Ӧ��Node����Ҫactive��
		if (vhierarchy_.node(v1_node_handle).is_active()) { 
			int lc_node_handle = vhierarchy_.node(v1_node_handle).lchild_node_handle();
			int rc_node_handle = vhierarchy_.node(v1_node_handle).rchild_node_handle();
			if (lc_node_handle != -1) { // rc_node_handle != -1 //����ڵ㲻��Ҷ�ӽڵ�,��������split
				//std::cout <<"can be splitted: " << lc_node_handle << ", " << rc_node_handle << "\n";//for test
				//vsplit(to)->(from, to), Ҳ����vspit(v1)->(v0, v1).
				TriMesh::Point v_from_point = vpoints_[vhierarchy_.node(lc_node_handle).point_handle()]; 
				TriMesh::VHandle v_from = simplified_mesh_.add_vertex(v_from_point);//���������from��
				TriMesh::VHandle v_to = v1;											//std::cout << vpoints_[vhierarchy_.node(lc_node_handle).point_handle()] << "\n";//for test
				TriMesh::VHandle v_l = TriMesh::VertexHandle(-1);//��ʾinvalid vertex handle
				TriMesh::VHandle v_r = TriMesh::VertexHandle(-1);
				int fc0 = vhierarchy_.node(v1_node_handle).fund_cut_node_handle0();//std::cout << "fc0 n h: " << fc0 << "\n";//for test, fundamental cut vertex vl^.
				if (fc0 != -1) {
					DGP::VHierarchyNode fund_cut_node = vhierarchy_.node(fc0);
					while (fund_cut_node.is_active() == false) {
						fund_cut_node = vhierarchy_.node(fund_cut_node.parent_node_handle());
					}
					v_l = fund_cut_node.vertex_handle();
				}
				int fc1 = vhierarchy_.node(v1_node_handle).fund_cut_node_handle1();//std::cout << "fc1 n h: " << fc1 << "\n";//for test, fundamental cut vertex vr^.
				if (fc1 != -1) {
					DGP::VHierarchyNode fund_cut_node = vhierarchy_.node(fc1); //std::cout << fund_cut_node.self_node_handle() << "\n";
					while (fund_cut_node.is_active() == false) {
						fund_cut_node = vhierarchy_.node(fund_cut_node.parent_node_handle());
					}  //std::cout << fund_cut_node.self_node_handle() << ", " << fund_cut_node.vertex_handle() << "\n";//for test
					v_r = fund_cut_node.vertex_handle();
				}
				//std::cout << "vsplit " << v_from << ", " << v_to << ", " << v_l << ", " << v_r << ".\n";//for test//
				simplified_mesh_.vertex_split(v_from, v_to, v_l, v_r);//std::cout << "0k.\n";

				vhierarchy_.node(v1_node_handle).set_active(false);
				vhierarchy_.node(lc_node_handle).set_vertex_handle(v_from);	simplified_mesh_.property(vp_node_handle_sm_, v_from) = lc_node_handle;
				vhierarchy_.node(lc_node_handle).set_active(true);
				vhierarchy_.node(rc_node_handle).set_vertex_handle(v_to);	simplified_mesh_.property(vp_node_handle_sm_, v1) = rc_node_handle;
				vhierarchy_.node(rc_node_handle).set_active(true);			

			}

		} // end of if.
	}
	simplified_mesh_.update_normals(); //MeshModel::update();
	set_mesh_toberendered(&simplified_mesh_, ORIGINAL);
	std::cout << "Now, the vertices num: " << simplified_mesh_.n_vertices() << ".\n";
	std::cout << "========����vsplit����==================\n";//������ʾ��40�����ȵ�=,ǰ���9��λ�ÿ�ʼ����,��������8��=.	
}

// ----------------------
// ----------------------
void DecimationModel::local_smooth_parameterization_triangledomain(TriMesh::FaceHandle _fh, int & n_vertices_relocated, int &n_vertices_changeface) //ÿ�����ĸ��ȱ���������Ϊƽ������
{
	TriMesh::FHandle f_it = _fh; //��f_it��������Ϊ�˷��������ǰ�Ĵ���
	std::vector<OpenMesh::Vec2d> vuv(18);
	vuv[0] = OpenMesh::Vec2d(1, 0);		vuv[1] = OpenMesh::Vec2d(1.5, 0.8660254);		vuv[2] = OpenMesh::Vec2d(0.5, 0.8660254);
	vuv[3] = OpenMesh::Vec2d(0, 0);		vuv[4] = OpenMesh::Vec2d(2, 0);					vuv[5] = OpenMesh::Vec2d(1, 1.7320508);
	vuv[6] = OpenMesh::Vec2d(-0.5, 0.8660254);		vuv[7] = OpenMesh::Vec2d(0.5, -0.8660254);
	vuv[8] = OpenMesh::Vec2d(1.5, -0.8660254);		vuv[9] = OpenMesh::Vec2d(2.5, 0.8660254);
	vuv[10] = OpenMesh::Vec2d(2, 1.7320508);		vuv[11] = OpenMesh::Vec2d(0, 1.7320508);
	vuv[12] = OpenMesh::Vec2d(-1, 0);				vuv[13] = OpenMesh::Vec2d(-0.5, -0.8660254);
	vuv[14] = OpenMesh::Vec2d(2.5, -0.8660254);		vuv[15] = OpenMesh::Vec2d(3, 0);
	vuv[16] = OpenMesh::Vec2d(1.5, 2.5980762);		vuv[17] = OpenMesh::Vec2d(0.5, 2.5980762);
	// �����Ƕ���f_it����.
	//////////////////////////////////////////////////////////////////////////
	// (1) Χ�Ƶ�ǰ��f_it���������� //����һ��������, ��4���ȱ����������.
	TriMesh::FHIter fh_it(mesh_, f_it);
	TriMesh::HalfedgeHandle h0 = fh_it.handle(); ++fh_it;
	TriMesh::HalfedgeHandle h1 = fh_it.handle(); ++fh_it;
	TriMesh::HalfedgeHandle h2 = fh_it.handle();
	TriMesh::VHandle v0(mesh_.to_vertex_handle(h0)), v1(mesh_.to_vertex_handle(h1)), v2(mesh_.to_vertex_handle(h2));
	mesh_.property(vuv_, v0) = vuv[0];//�⼸�����ǲ����������������ص���, ���Կ������ھ�ָ��ȥ����.
	mesh_.property(vuv_, v1) = vuv[1];
	mesh_.property(vuv_, v2) = vuv[2];

	TriMesh::FaceHandle f0(-1), f1(-1), f2(-1);
	TriMesh::VHandle v3, v4, v5;//(-1)
	if (mesh_.is_boundary(mesh_.opposite_halfedge_handle(h0)) == false) {
		f0 = mesh_.face_handle(mesh_.opposite_halfedge_handle(h0));
		v3 = mesh_.to_vertex_handle(mesh_.next_halfedge_handle(mesh_.opposite_halfedge_handle(h0)));
		//mesh_.property(vuv_, v3) = vuv[3];//������Ϊv3/v4/v5������ָ����ͬ�ĵ�,�����Ȳ�������������,�����ظ�.
	}
	if (mesh_.is_boundary(mesh_.opposite_halfedge_handle(h1)) == false) {
		f1 = mesh_.face_handle(mesh_.opposite_halfedge_handle(h1));
		v4 = mesh_.to_vertex_handle(mesh_.next_halfedge_handle(mesh_.opposite_halfedge_handle(h1)));				
	}
	if (mesh_.is_boundary(mesh_.opposite_halfedge_handle(h2)) == false) {  //h2�Ķ԰�߲��Ǳ߽�
		f2 = mesh_.face_handle(mesh_.opposite_halfedge_handle(h2));
		v5 = mesh_.to_vertex_handle(mesh_.next_halfedge_handle(mesh_.opposite_halfedge_handle(h2)));			
	} //��֮���f0,f1,f2��Ȼ������.is_valid() == fasle��, ���Ǳ߽���Ľ��.
	TriMesh::FHandle f3, f4, f5, f6, f7, f8;//�������β����������extension
	TriMesh::HHandle h3, h4, h5, h6, h7, h8;
	TriMesh::VHandle v6, v7, v8, v9, v10, v11;
	if (f0.is_valid() == true) {
		TriMesh::HHandle oh0 = mesh_.opposite_halfedge_handle(h0);
		h3 = mesh_.opposite_halfedge_handle(mesh_.next_halfedge_handle(oh0)); h4 = mesh_.opposite_halfedge_handle(mesh_.prev_halfedge_handle(oh0));
		if (mesh_.is_boundary(h3) == false) {
			f3 = mesh_.face_handle(h3);		v6 = mesh_.to_vertex_handle(mesh_.next_halfedge_handle(h3));					
		}
		if (mesh_.is_boundary(h4) == false) {
			f4 = mesh_.face_handle(h4);		v7 = mesh_.to_vertex_handle(mesh_.next_halfedge_handle(h4));
		}							
	}
	if (f1.is_valid() == true) {
		TriMesh::HHandle oh1 = mesh_.opposite_halfedge_handle(h1);
		h5 = mesh_.opposite_halfedge_handle(mesh_.next_halfedge_handle(oh1)); h6 = mesh_.opposite_halfedge_handle(mesh_.prev_halfedge_handle(oh1));
		if (mesh_.is_boundary(h5) == false) {
			f5 = mesh_.face_handle(h5);		v8 = mesh_.to_vertex_handle(mesh_.next_halfedge_handle(h5));								
		}
		if (mesh_.is_boundary(h6) == false) {
			f6 = mesh_.face_handle(h6);		v9 = mesh_.to_vertex_handle(mesh_.next_halfedge_handle(h6));	
		} 
	}
	if (f2.is_valid() == true) {
		TriMesh::HHandle oh2 = mesh_.opposite_halfedge_handle(h2);
		h7 = mesh_.opposite_halfedge_handle(mesh_.next_halfedge_handle(oh2)); h8 = mesh_.opposite_halfedge_handle(mesh_.prev_halfedge_handle(oh2));
		if (mesh_.is_boundary(h7) == false) {
			f7 = mesh_.face_handle(h7);		v10 = mesh_.to_vertex_handle(mesh_.next_halfedge_handle(h7));							
		}
		if (mesh_.is_boundary(h8) == false) {
			f8 = mesh_.face_handle(h8);		v11 = mesh_.to_vertex_handle(mesh_.next_halfedge_handle(h8));	
		} 
	}
	TriMesh::FaceHandle f9, f10, f11, f12, f13, f14;
	TriMesh::HalfedgeHandle h9, h10, h11, h12, h13, h14;
	TriMesh::VertexHandle v12, v13, v14, v15, v16, v17;
	if (f3.is_valid()) {
		h9 = mesh_.opposite_halfedge_handle(mesh_.prev_halfedge_handle(h3));
		if (mesh_.is_boundary(h9) == false) { f9 = mesh_.face_handle(h9); v12 = mesh_.to_vertex_handle(mesh_.next_halfedge_handle(h9)); }
	}
	if (f4.is_valid()) {
		h10 = mesh_.opposite_halfedge_handle(mesh_.next_halfedge_handle(h4));
		if (mesh_.is_boundary(h10) == false) { f10 = mesh_.face_handle(h10); v13 = mesh_.to_vertex_handle(mesh_.next_halfedge_handle(h10)); }
	}
	if (f5.is_valid()) {
		h11 = mesh_.opposite_halfedge_handle(mesh_.prev_halfedge_handle(h5));
		if (mesh_.is_boundary(h11) == false) { f11 = mesh_.face_handle(h11); v14 = mesh_.to_vertex_handle(mesh_.next_halfedge_handle(h11)); }
	}
	if (f6.is_valid()) {
		h12 = mesh_.opposite_halfedge_handle(mesh_.next_halfedge_handle(h6));
		if (mesh_.is_boundary(h12) == false) { f12 = mesh_.face_handle(h12); v15 = mesh_.to_vertex_handle(mesh_.next_halfedge_handle(h12)); }
	}
	if (f7.is_valid()) {
		h13 = mesh_.opposite_halfedge_handle(mesh_.prev_halfedge_handle(h7));
		if (mesh_.is_boundary(h13) == false) { f13 = mesh_.face_handle(h13); v16 = mesh_.to_vertex_handle(mesh_.next_halfedge_handle(h13)); }
	}
	if (f8.is_valid()) {
		h14 = mesh_.opposite_halfedge_handle(mesh_.next_halfedge_handle(h8));
		if (mesh_.is_boundary(h14) == false) { f14 = mesh_.face_handle(h14); v17 = mesh_.to_vertex_handle(mesh_.next_halfedge_handle(h14)); }
	}
	TriMesh::FaceHandle f15, f16, f17, f18, f19, f20;
	TriMesh::HalfedgeHandle h15, h16, h17, h18, h19, h20;
	TriMesh::VertexHandle v18, v19, v20, v21, v22, v23;
	/*if (f4.is_valid()) {
	h15 = mesh_.opposite_halfedge_handle(mesh_.prev_halfedge_handle(h4));
	if (mesh_.is_boundary(h15) == false) { f15 = mesh_.face_handle(h15); v18 = mesh_.to_vertex_handle(mesh_.next_halfedge_handle(h15));}
	}
	if (f5.is_valid()) {
	h16 = mesh_.opposite_halfedge_handle(mesh_.next_halfedge_handle(h5));
	if (mesh_.is_boundary(h16) == false) { f16 = mesh_.face_handle(h16); v19 = mesh_.to_vertex_handle(mesh_.next_halfedge_handle(h16));}
	}
	if (f6.is_valid()) {
	h17 = mesh_.opposite_halfedge_handle(mesh_.prev_halfedge_handle(h6));
	if (mesh_.is_boundary(h17) == false) { f17 = mesh_.face_handle(h17); v20 = mesh_.to_vertex_handle(mesh_.next_halfedge_handle(h17));}
	}
	if (f7.is_valid()) {
	h18 = mesh_.opposite_halfedge_handle(mesh_.next_halfedge_handle(h7));
	if (mesh_.is_boundary(h18) == false) { f18 = mesh_.face_handle(h18); v21 = mesh_.to_vertex_handle(mesh_.next_halfedge_handle(h18));}
	}
	if (f8.is_valid()) {
	h19 = mesh_.opposite_halfedge_handle(mesh_.prev_halfedge_handle(h8));
	if (mesh_.is_boundary(h19) == false) { f19 = mesh_.face_handle(h19); v22 = mesh_.to_vertex_handle(mesh_.next_halfedge_handle(h19));}
	}
	if (f3.is_valid()) {
	h20 = mesh_.opposite_halfedge_handle(mesh_.next_halfedge_handle(h3));
	if (mesh_.is_boundary(h20) == false) { f20 = mesh_.face_handle(h20); v23 = mesh_.to_vertex_handle(mesh_.next_halfedge_handle(h20));}
	}*/

	// ���(1)�Բ�����Ľ���.
	//////////////////////////////////////////////////////////////////////////
	// (2)
	//��������������ϵ����ж���
	std::vector<TriMesh::VertexHandle> fvset = mesh_.property(fvset_, f_it);//ƽ�������������f_it���ϵĲ���������.
	//std::cout << "fvset_'s size " << fvset.size() << ", ";// std::endl; // for test// 
	for (std::vector<TriMesh::VertexHandle>::iterator fvset_it(fvset.begin()), it_end(fvset.end()); fvset_it != it_end; ++fvset_it) {
		//ÿ�� ����һ���������������ϵĶ���*fvset_it.
		//std::cout << "��δ���Ķ�����*fvset_it: " << *fvset_it << ".\n"; //for test
		if (mesh_.property(vf_, *fvset_it) != f_it) { // only to confirm
			std::cout << "Error, para check, faces not match.\n";
		}

		if (mesh2_.is_boundary(*fvset_it)) continue;//����Ǳ߽��ϵĶ���(��������ʱ��Ҳ���������߽��ϵ�)�Ͳ�Ҫ����smooth parameterize.
		if (mesh_.property(vp_type_, *fvset_it) == DGP::CREASE_VFT)  continue;// crease vertexҲͬ������Ҫsmooth					 

		// ���������*fvset_it��original mesh������һ��������Ķ���vj, ���ж�vj�ǲ��Ƕ�������������ε�����(������4��������)����.
		// ����ǵĻ�, ��Ч, �͸������������������.
		bool inside_valid = true; 
		double wij= 0;
		double sum_wij = 0.0;
		OpenMesh::Vec2d sum_wh(0, 0, 0);            
		// �ж�*fvset_it��һ�����vj�Ƿ���Ч
		for (TriMesh::VertexOHalfedgeIter voh_it(mesh2_, *fvset_it); voh_it && inside_valid; ++voh_it) {
			TriMesh::VHandle vj = mesh2_.to_vertex_handle(voh_it.handle());  //���ڶ���*fvset_it�����ڱ���, ���������߶������ڱ߽��.
			//�ھӵ�vjҪô�Ǳ�deleted�˲�����4������, Ҫô��û�б�deleted����������������, ����������Ч��.
			//std::cout << "vj " << vj << ", ";//for test
			if (mesh_.status(vj).deleted()) { //�ڽӵ�vj�Ǳ�deleted��Ҳ�����б���������base complex��ĳһ�����ϵ�.
				TriMesh::FaceHandle pf = mesh_.property(vf_, vj);//�ڵ�vj��initial parameterize����������������
				if (pf == f_it) {														
				} else if (pf == f0) { mesh_.property(vuv_, v3) = vuv[3];
				} else if (pf == f1) { mesh_.property(vuv_, v4) = vuv[4];
				} else if (pf == f2) { mesh_.property(vuv_, v5) = vuv[5]; 
				} else if (pf == f3) { mesh_.property(vuv_, v3) = vuv[3]; mesh_.property(vuv_, v6) = vuv[6];
				} else if (pf == f4) { mesh_.property(vuv_, v3) = vuv[3]; mesh_.property(vuv_, v7) = vuv[7];
				} else if (pf == f5) { mesh_.property(vuv_, v4) = vuv[4]; mesh_.property(vuv_, v8) = vuv[8];
				} else if (pf == f6) { mesh_.property(vuv_, v4) = vuv[4]; mesh_.property(vuv_, v9) = vuv[9];
				} else if (pf == f7) { mesh_.property(vuv_, v5) = vuv[5]; mesh_.property(vuv_, v10) = vuv[10];
				} else if (pf == f8) { mesh_.property(vuv_, v5) = vuv[5]; mesh_.property(vuv_, v11) = vuv[11];
				} else if (pf == f9) { 
					mesh_.property(vuv_, v3) = vuv[3]; mesh_.property(vuv_, v6) = vuv[6]; mesh_.property(vuv_, v12) = vuv[12];
				} else if (pf == f10) { 
					mesh_.property(vuv_, v3) = vuv[3]; mesh_.property(vuv_, v7) = vuv[7]; mesh_.property(vuv_, v13) = vuv[13];
				} else if (pf == f11) { 
					mesh_.property(vuv_, v4) = vuv[4]; mesh_.property(vuv_, v8) = vuv[8]; mesh_.property(vuv_, v14) = vuv[14];
				} else if (pf == f12) { 
					mesh_.property(vuv_, v4) = vuv[4]; mesh_.property(vuv_, v9) = vuv[9]; mesh_.property(vuv_, v15) = vuv[15];
				} else if (pf == f13) { 
					mesh_.property(vuv_, v5) = vuv[5]; mesh_.property(vuv_, v10) = vuv[10]; mesh_.property(vuv_, v16) = vuv[16]; 
				} else if (pf == f14) { 
					mesh_.property(vuv_, v5) = vuv[5]; mesh_.property(vuv_, v11) = vuv[11]; mesh_.property(vuv_, v17) = vuv[17];	 
				} /*else if (pf == f15) { 
				  mesh_.property(vuv_, v7) = vuv[7]; mesh_.property(vuv_, v18) = OpenMesh::Vec2d(1, -0.8660254);
				  } else if (pf == f16) { 
				  mesh_.property(vuv_, v8) = OpenMesh::Vec2d(1.5, -0.8660254); mesh_.property(vuv_, v19) = OpenMesh::Vec2d(1, -0.8660254);	 
				  } else if (pf == f17) { 
				  mesh_.property(vuv_, v9) = OpenMesh::Vec2d(2.5, 0.8660254); mesh_.property(vuv_, v20) = OpenMesh::Vec2d(2.25, 1.2990381);
				  } else if (pf == f18) { 
				  mesh_.property(vuv_, v10) = OpenMesh::Vec2d(2, 1.7320508); mesh_.property(vuv_, v21) = OpenMesh::Vec2d(2.25, 1.2990381);	 
				  } else if (pf == f19) { 
				  mesh_.property(vuv_, v11) = OpenMesh::Vec2d(0, 1.7320508); mesh_.property(vuv_, v22) = OpenMesh::Vec2d(-0.25, 1.2990381);	 
				  } else if (pf == f20) { 
				  mesh_.property(vuv_, v6) = OpenMesh::Vec2d(-0.5, 0.8660254); mesh_.property(vuv_, v23) = OpenMesh::Vec2d(-0.25, 1.2990381);	 
				  } */else inside_valid = false; 

				if (inside_valid == true) { // == falseʱ���û�б�Ҫ���������������.
					mesh_.property(vuv_, vj)
						= mesh_.property(vuv_, mesh_.property(vf0_, vj)) * (mesh_.property(vbc_, vj))[0]
					+ mesh_.property(vuv_, mesh_.property(vf1_, vj)) * (mesh_.property(vbc_, vj))[1]
					+ mesh_.property(vuv_, mesh_.property(vf2_, vj)) * (mesh_.property(vbc_, vj))[2];
				} else { //��֮ǰ�ļ���в�����9����Ч����, ���ﻹҪ��һ�����
					if (mesh_.property(vp_type_, vj) == DGP::CREASE_VFT) {

						if (h3.is_valid()) {
							if (mesh_.property(vp_feature_edge_, vj) == mesh_.edge_handle(mesh_.next_halfedge_handle(h3))) {
								mesh_.property(vuv_, v3) = vuv[3]; mesh_.property(vuv_, v6) = vuv[6]; 
								inside_valid = true;
							}
							if (mesh_.property(vp_feature_edge_, vj) == mesh_.edge_handle(mesh_.prev_halfedge_handle(h3))) { //v2��vuv_������
								//�������Ӧ�ð�����f3 or f9������
							}
						}
						if (h4.is_valid()) {
							if (mesh_.property(vp_feature_edge_, vj) == mesh_.edge_handle(mesh_.next_halfedge_handle(h4))
								|| mesh_.property(vp_feature_edge_, vj) == mesh_.edge_handle(mesh_.prev_halfedge_handle(h4))) { //v0��vuv_������
									mesh_.property(vuv_, v3) = vuv[3]; mesh_.property(vuv_, v7) = vuv[7];
									inside_valid = true;//if ((*fvset_it).idx() == 6424) { std::cout << "tri smooth crease vj 4.\n"; }
							}
						}
						if (h5.is_valid()) {
							if (mesh_.property(vp_feature_edge_, vj) == mesh_.edge_handle(mesh_.next_halfedge_handle(h5))
								|| mesh_.property(vp_feature_edge_, vj) == mesh_.edge_handle(mesh_.prev_halfedge_handle(h5))) { //v0��vuv_������
									mesh_.property(vuv_, v4) = vuv[4]; mesh_.property(vuv_, v8) = vuv[8];
									inside_valid = true;//if ((*fvset_it).idx() == 6424) { std::cout << "tri smooth crease vj 5.\n"; }
							}
						}
						if (h6.is_valid()) {
							if (mesh_.property(vp_feature_edge_, vj) == mesh_.edge_handle(mesh_.next_halfedge_handle(h6))
								|| mesh_.property(vp_feature_edge_, vj) == mesh_.edge_handle(mesh_.prev_halfedge_handle(h6))) { //v1��vuv_������
									mesh_.property(vuv_, v4) = vuv[4]; mesh_.property(vuv_, v9) = vuv[9];
									inside_valid = true;//if ((*fvset_it).idx() == 6424) { std::cout << "tri smooth crease vj 6.\n"; }
							}
						}
						if (h7.is_valid()) {
							if (mesh_.property(vp_feature_edge_, vj) == mesh_.edge_handle(mesh_.next_halfedge_handle(h7))
								|| mesh_.property(vp_feature_edge_, vj) == mesh_.edge_handle(mesh_.prev_halfedge_handle(h7))) { //v1��vuv_������
									mesh_.property(vuv_, v5) = vuv[5]; mesh_.property(vuv_, v10) = vuv[10];
									inside_valid = true;//if ((*fvset_it).idx() == 6424) { std::cout << "tri smooth crease vj 7.\n"; }
							}
						}
						if (h8.is_valid()) {
							if (mesh_.property(vp_feature_edge_, vj) == mesh_.edge_handle(mesh_.next_halfedge_handle(h8))
								|| mesh_.property(vp_feature_edge_, vj) == mesh_.edge_handle(mesh_.prev_halfedge_handle(h8))) { //v2��vuv_������
									mesh_.property(vuv_, v5) = vuv[5]; mesh_.property(vuv_, v11) = vuv[11];
									inside_valid = true;				//if ((*fvset_it).idx() == 6424) { std::cout << 8 << ".\n"; }
							}
						}
						if (h9.is_valid()) {
							if (mesh_.property(vp_feature_edge_, vj) == mesh_.edge_handle(mesh_.next_halfedge_handle(h9))
								|| mesh_.property(vp_feature_edge_, vj) == mesh_.edge_handle(mesh_.prev_halfedge_handle(h9))) {
									mesh_.property(vuv_, v3) = vuv[3]; mesh_.property(vuv_, v6) = vuv[6];
									mesh_.property(vuv_, v12) = vuv[12]; inside_valid = true;
							}
						}
						if (h10.is_valid()) {
							if (mesh_.property(vp_feature_edge_, vj) == mesh_.edge_handle(mesh_.next_halfedge_handle(h10))
								|| mesh_.property(vp_feature_edge_, vj) == mesh_.edge_handle(mesh_.prev_halfedge_handle(h10))) {
									mesh_.property(vuv_, v3) = vuv[3]; mesh_.property(vuv_, v7) = vuv[7];
									mesh_.property(vuv_, v13) = vuv[13]; inside_valid = true;
							}
						}
						if (h11.is_valid()) {
							if (mesh_.property(vp_feature_edge_, vj) == mesh_.edge_handle(mesh_.next_halfedge_handle(h11))
								|| mesh_.property(vp_feature_edge_, vj) == mesh_.edge_handle(mesh_.prev_halfedge_handle(h11))) {
									mesh_.property(vuv_, v4) = vuv[4]; mesh_.property(vuv_, v8) = vuv[8];
									mesh_.property(vuv_, v14) = vuv[14]; inside_valid = true;
							}
						}
						if (h12.is_valid()) {
							if (mesh_.property(vp_feature_edge_, vj) == mesh_.edge_handle(mesh_.next_halfedge_handle(h12))
								|| mesh_.property(vp_feature_edge_, vj) == mesh_.edge_handle(mesh_.prev_halfedge_handle(h12))) {
									mesh_.property(vuv_, v4) = vuv[4]; mesh_.property(vuv_, v9) = vuv[9];
									mesh_.property(vuv_, v15) = vuv[15]; inside_valid = true;
							}
						}
						if (h13.is_valid()) {
							if (mesh_.property(vp_feature_edge_, vj) == mesh_.edge_handle(mesh_.next_halfedge_handle(h13))
								|| mesh_.property(vp_feature_edge_, vj) == mesh_.edge_handle(mesh_.prev_halfedge_handle(h13))) {
									mesh_.property(vuv_, v5) = vuv[5]; mesh_.property(vuv_, v10) = vuv[10];
									mesh_.property(vuv_, v16) = vuv[16]; inside_valid = true;
							}
						}
						if (h14.is_valid()) {
							if (mesh_.property(vp_feature_edge_, vj) == mesh_.edge_handle(mesh_.next_halfedge_handle(h14))
								|| mesh_.property(vp_feature_edge_, vj) == mesh_.edge_handle(mesh_.prev_halfedge_handle(h14))) {
									mesh_.property(vuv_, v5) = vuv[5]; mesh_.property(vuv_, v11) = vuv[11];
									mesh_.property(vuv_, v17) = vuv[17]; inside_valid = true;
							}
						} /*////////////////////////////////////////////
						  if (h15.is_valid()) {
						  if (mesh_.property(vp_feature_edge_, vj) == mesh_.edge_handle(mesh_.next_halfedge_handle(h15))
						  || mesh_.property(vp_feature_edge_, vj) == mesh_.edge_handle(mesh_.prev_halfedge_handle(h15))) {
						  mesh_.property(vuv_, v7) = vuv[7]; mesh_.property(vuv_, v18) = OpenMesh::Vec2d(1, -0.8660254);
						  inside_valid = true;
						  }
						  }
						  if (h16.is_valid()) {
						  if (mesh_.property(vp_feature_edge_, vj) == mesh_.edge_handle(mesh_.next_halfedge_handle(h16))
						  || mesh_.property(vp_feature_edge_, vj) == mesh_.edge_handle(mesh_.prev_halfedge_handle(h16))) {
						  mesh_.property(vuv_, v8) = OpenMesh::Vec2d(1.5, -0.8660254); mesh_.property(vuv_, v19) = OpenMesh::Vec2d(1, -0.8660254);
						  inside_valid = true;
						  }
						  }
						  if (h17.is_valid()) {
						  if (mesh_.property(vp_feature_edge_, vj) == mesh_.edge_handle(mesh_.next_halfedge_handle(h17))
						  || mesh_.property(vp_feature_edge_, vj) == mesh_.edge_handle(mesh_.prev_halfedge_handle(h17))) {
						  mesh_.property(vuv_, v9) = OpenMesh::Vec2d(2.5, 0.8660254); mesh_.property(vuv_, v20) = OpenMesh::Vec2d(2.25, 1.2990381);
						  inside_valid = true;
						  }
						  }
						  if (h18.is_valid()) {
						  if (mesh_.property(vp_feature_edge_, vj) == mesh_.edge_handle(mesh_.next_halfedge_handle(h18))
						  || mesh_.property(vp_feature_edge_, vj) == mesh_.edge_handle(mesh_.prev_halfedge_handle(h18))) {
						  mesh_.property(vuv_, v10) = OpenMesh::Vec2d(2, 1.7320508); mesh_.property(vuv_, v21) = OpenMesh::Vec2d(2.25, 1.2990381);
						  inside_valid = true;
						  }
						  }
						  if (h19.is_valid()) {
						  if (mesh_.property(vp_feature_edge_, vj) == mesh_.edge_handle(mesh_.next_halfedge_handle(h19))
						  || mesh_.property(vp_feature_edge_, vj) == mesh_.edge_handle(mesh_.prev_halfedge_handle(h19))) {
						  mesh_.property(vuv_, v11) = OpenMesh::Vec2d(0, 1.7320508); mesh_.property(vuv_, v22) = OpenMesh::Vec2d(-0.25, 1.2990381);
						  inside_valid = true;
						  }
						  }
						  if (h20.is_valid()) {
						  if (mesh_.property(vp_feature_edge_, vj) == mesh_.edge_handle(mesh_.next_halfedge_handle(h20))
						  || mesh_.property(vp_feature_edge_, vj) == mesh_.edge_handle(mesh_.prev_halfedge_handle(h20))) {
						  mesh_.property(vuv_, v6) = OpenMesh::Vec2d(-0.5, 0.8660254); mesh_.property(vuv_, v23) = OpenMesh::Vec2d(-0.25, 1.2990381);
						  inside_valid = true;
						  }
						  }*/
						if (inside_valid) {/*
										   mesh_.property(vuv_, vj) = mesh_.property(vuv_, mesh_.property(vf0_, vj)) * (mesh_.property(vbc_, vj))[0]
										   + vuv[3]  // vbc_[1] == 0 // ���ﷸ�˸����ʹ���, crease���������crease edge��ֻ����������vbc_, ���Ѿ����Ѿ�����vbc_[1] == 0!
										   + mesh_.property(vuv_, mesh_.property(vf2_, vj)) * (mesh_.property(vbc_, vj))[2];*/
							mesh_.property(vuv_, vj) = mesh_.property(vuv_, mesh_.property(vf0_, vj)) * (mesh_.property(vbc_, vj))[0]
							+ mesh_.property(vuv_, mesh_.property(vf1_, vj)) * (mesh_.property(vbc_, vj))[1]
							+ mesh_.property(vuv_, mesh_.property(vf2_, vj)) * (mesh_.property(vbc_, vj))[2];
							// ����������һ��vbc_[i]������, �պ����Ӧ�Ķ������������Ҳ��������.
						}

					} //else inside_valid = false; �������Ѿ�==false  // end of if-else vj is crease 
				} // end of if-else. ��deleted�Ļ��Ƿ���Ч
			} else { //�ڽӵ�vjû�б�deleted��Ҳ����������base complex�ϵ�ĳһ������.
				if ((vj == v0) || (vj == v1) || (vj == v2)) { //mesh_.property(vuv_, vj);��vj�͵�������������Ļ���ô�����������Ѿ����������õ���.
				} else if (vj == v3) { mesh_.property(vuv_, v3) = vuv[3];
				} else if (vj == v4) { mesh_.property(vuv_, v4) = vuv[4];
				} else if (vj == v5) { mesh_.property(vuv_, v5) = vuv[5]; 
				} else if (vj == v6) { mesh_.property(vuv_, v3) = vuv[3]; mesh_.property(vuv_, v6) = vuv[6]; 
				} else if (vj == v7) { mesh_.property(vuv_, v3) = vuv[3]; mesh_.property(vuv_, v7) = vuv[7]; 			
				} else if (vj == v8) { mesh_.property(vuv_, v4) = vuv[4]; mesh_.property(vuv_, v8) = vuv[8]; 
				} else if (vj == v9) { mesh_.property(vuv_, v4) = vuv[4]; mesh_.property(vuv_, v9) = vuv[9];  
				} else if (vj == v10) { mesh_.property(vuv_, v5) = vuv[5]; mesh_.property(vuv_, v10) = vuv[10]; 
				} else if (vj == v11) { mesh_.property(vuv_, v5) = vuv[5]; mesh_.property(vuv_, v11) = vuv[11]; 
				} else if (vj == v12) { mesh_.property(vuv_, v12) = vuv[12]; 
				} else if (vj == v13) { mesh_.property(vuv_, v13) = vuv[13]; 
				} else if (vj == v14) { mesh_.property(vuv_, v14) = vuv[14]; 
				} else if (vj == v15) { mesh_.property(vuv_, v15) = vuv[15]; 
				} else if (vj == v16) { mesh_.property(vuv_, v16) = vuv[16]; 
				} else if (vj == v17) { mesh_.property(vuv_, v17) = vuv[17]; 
				} else { inside_valid = false; }
			}

			if (inside_valid) { //vi������һ��vj��Ч//
				//std::cout << "vj is valid, " << mesh_.property(vuv_, vj) << ";\t";//for test
				//������voh_it.handle()����Ӧ��ϵ��wij, ע������ʹ�õ���mean value coordiante, ������harmonic map�е�cotangent weight
				sum_wij += mesh2_.property(hmvc_, voh_it);
				sum_wh += mesh2_.property(hmvc_, voh_it) * mesh_.property(vuv_, vj);							
			} else { // �����Ӧ�ò���Ŷ�, �Ͼ�չ�������Ѿ�������.
				std::cout << "������ƽ��In face " << f_it << "չ��������: " << f0 << " " << f1 << " " << f2 << " " // for test
					<< f3 << " " << f4 << " " << f5 << " " << f6 << " " << f7 << " " << f8 << " "
					<< f9 << " " << f10 << " " << f11 << " " << f12 << " " << f13 << " " << f14  << ".\n"; 
				std::cout << "(vi, vj)=(" << *fvset_it << ", " << vj << "), vj's vf_" << mesh_.property(vf_, vj) << ", vj's type" << mesh_.property(vp_type_, vj) << ".\n";  
				TriMesh::FaceFaceIter fit = mesh_.ff_iter(mesh_.property(vf_, vj)); std::cout << fit.handle() << " ";
				++ fit; std::cout << fit.handle() << " ";++ fit; std::cout << fit.handle() << ":vj's one-ring faces.\n";
				/**/
			}
		} // end of for.ѭ������˶���*fvset_it��һ���򶥵�vj, ���ж϶���*fvset_it�Ƿ���Ч,Ҳ����˵�ǲ��Ƕ���һ��������.

		if (inside_valid == true) { 
			// ����vi������vj����Ч�Ļ������vi(Ҳ����*fvset_it)������²���������, ���ҳ��������ĸ���������.
			OpenMesh::Vec2d old_uv = mesh_.property(vuv_, *fvset_it); 
			mesh_.property(vuv_, *fvset_it) = sum_wh * (1.0/sum_wij);
			// std::cout << "vi " << mesh_.property(vuv_, *fvset_it) << std::endl; // for test
			// �����ж϶���*fvset_it�ǲ��������ĸ���(������չ��10����)����һ����������.
			std::vector<TriMesh::FHandle> fvector(22);
			fvector[0] = f_it; fvector[1] = f0; fvector[2] = f1; fvector[3] = f2;
			fvector[4] = f3; fvector[5] = f4; fvector[6] = f5; fvector[7] = f6;fvector[8] = f7; fvector[9] = f8; 
			fvector[10] = f9; fvector[11] = f10; fvector[12] = f11; fvector[13] = f12;fvector[14] = f13; fvector[15] = f14; 
			fvector[16] = f15; fvector[17] = f16; fvector[18] = f17; fvector[19] = f18;fvector[20] = f19; fvector[21] = f20; 
			bool is_face_located = false;
			for (size_t i = 0; i < 16 && is_face_located == false; i++) { //ѭ���ĸ����ܵ��� //22
				if ((fvector[i]).is_valid() == false ) continue;//f0, f1, f2�������涼���ܲ�һ�����ڵ�.
				TriMesh::VHandle vx, vy, vz;
				if (i == 0) {// fvector[0]����Ӧ����������
					vx = v0; vy = v1; vz = v2; 
				} else if (i == 1) { // fvector[1]����Ӧ����������
					vx = v0; vy = v2; vz = v3; mesh_.property(vuv_, v3) = vuv[3];
				} else if (i == 2) { // fvector[2]����Ӧ����������
					vx = v0; vy = v4; vz = v1; mesh_.property(vuv_, v4) = vuv[4];
				} else if (i == 3) { // fvector[3]����Ӧ����������
					vx = v1; vy = v5; vz = v2; mesh_.property(vuv_, v5) = vuv[5]; 
				} else if (i == 4) { vx = v3; vy = v2; vz = v6; // f3
				mesh_.property(vuv_, v3) = vuv[3]; mesh_.property(vuv_, v6) = vuv[6];
				} else if (i == 5) { vx = v0; vy = v3; vz = v7; // f4
				mesh_.property(vuv_, v3) = vuv[3]; mesh_.property(vuv_, v7) = vuv[7];
				} else if (i == 6) { vx = v4; vy = v0; vz = v8; // f5
				mesh_.property(vuv_, v4) = vuv[4]; mesh_.property(vuv_, v8) = vuv[8];
				} else if (i == 7) { vx = v1; vy = v4; vz = v9; // f6
				mesh_.property(vuv_, v4) = vuv[4]; mesh_.property(vuv_, v9) = vuv[9];  
				} else if (i == 8) { vx = v5; vy = v1; vz = v10; //f7
				mesh_.property(vuv_, v5) = vuv[5]; mesh_.property(vuv_, v10) = vuv[10];
				} else if (i == 9) { vx = v2; vy = v5; vz = v11; //f8
				mesh_.property(vuv_, v5) = vuv[5]; mesh_.property(vuv_, v11) = vuv[11];
				} else if (i == 10) { vx = v3; vy = v6; vz = v12; //f9
				mesh_.property(vuv_, v3) = vuv[3]; mesh_.property(vuv_, v6) = vuv[6]; mesh_.property(vuv_, v12) = vuv[12];
				} else if (i == 11) { vx = v7; vy = v3; vz = v13; //f10
				mesh_.property(vuv_, v3) = vuv[3]; mesh_.property(vuv_, v7) = vuv[7]; mesh_.property(vuv_, v13) = vuv[13]; 
				} else if (i == 12) { vx = v4; vy = v8; vz = v14; //f11
				mesh_.property(vuv_, v4) = vuv[4]; mesh_.property(vuv_, v8) = vuv[8]; mesh_.property(vuv_, v14) = vuv[14];
				} else if (i == 13) { vx = v9; vy = v4; vz = v15; //f12
				mesh_.property(vuv_, v4) = vuv[4]; mesh_.property(vuv_, v9) = vuv[9]; mesh_.property(vuv_, v15) = vuv[15]; 
				} else if (i == 14) { vx = v5; vy = v10; vz = v16; //f13
				mesh_.property(vuv_, v5) = vuv[5]; mesh_.property(vuv_, v10) = vuv[10]; mesh_.property(vuv_, v16) = vuv[16];
				} else if (i == 15) { vx = v11; vy = v5; vz = v17; //f14
				mesh_.property(vuv_, v5) = vuv[5]; mesh_.property(vuv_, v11) = vuv[11]; mesh_.property(vuv_, v17) = vuv[17]; 
				} else if (i == 16) { vx = v0; vy = v7; vz = v18; // 
				mesh_.property(vuv_, v7) = vuv[7]; mesh_.property(vuv_, v18) = OpenMesh::Vec2d(1, -0.8660254);
				} else if (i == 17) { vx = v8; vy = v0; vz = v19; // 
				mesh_.property(vuv_, v8) = vuv[8]; mesh_.property(vuv_, v19) = OpenMesh::Vec2d(1, -0.8660254);
				} else if (i == 18) { vx = v1; vy = v9; vz = v20; // 
				mesh_.property(vuv_, v9) = vuv[9]; mesh_.property(vuv_, v20) = OpenMesh::Vec2d(2.25, 1.2990381);
				} else if (i == 19) { vx = v10; vy = v1; vz = v21; // 
				mesh_.property(vuv_, v10) = vuv[10]; mesh_.property(vuv_, v21) = OpenMesh::Vec2d(2.25, 1.2990381);
				} else if (i == 20) { vx = v2; vy = v11; vz = v22; // 
				mesh_.property(vuv_, v11) = vuv[11]; mesh_.property(vuv_, v22) = OpenMesh::Vec2d(-0.25, 1.2990381);
				} else if (i == 21) { vx = v6; vy = v2; vz = v23; // 
				mesh_.property(vuv_, v6) = vuv[6]; mesh_.property(vuv_, v23) = OpenMesh::Vec2d(-0.25, 1.2990381);
				} 

				OpenMesh::Vec3d bc 
					= DGP::calc_barycentric_coordinates(mesh_.property(vuv_, vx), mesh_.property(vuv_, vy), mesh_.property(vuv_, vz), mesh_.property(vuv_, *fvset_it));
				//std::cout << "BC: " << bc << " of face " << i << std::endl; //for test
				if (fabs(bc[0]) < 1.0e-10) bc[0] = 0; if (fabs(bc[1]) < 1.0e-10) bc[1] = 0; if (fabs(bc[2]) < 1.0e-10) bc[2] = 0;
				if (!(bc[0]<0) && !(bc[1]<0) && !(bc[2]<0)) { //Ҫ�����bc��������������>=0
					//std::cout << "BC: " << bc << std::endl; // for test
					// �������о���ÿһ�β�������Ӧ�û�õĻ��Ǹ��ĵ���Ϣ.
					mesh_.property(vf0_, *fvset_it) = vx; mesh_.property(vf1_, *fvset_it) = vy; mesh_.property(vf2_, *fvset_it) = vz;
					mesh_.property(vbc_, *fvset_it) = bc;
					if (f_it != fvector[i]) { //mesh_.property(vf_, *fvset_it) �͵���f_it
						n_vertices_changeface++;//ǰ�����Ƿ�һ��

						// ɾȥ����ԭ���Ǹ���(f_it)�ϵļ�¼
						std::vector<TriMesh::VHandle>::iterator to_be_earsed 
							= remove((mesh_.property(fvset_, f_it)).begin(), (mesh_.property(fvset_, f_it)).end(), TriMesh::VHandle((*fvset_it).idx()));
						(mesh_.property(fvset_, f_it)).erase(to_be_earsed, (mesh_.property(fvset_, f_it)).end());

						(mesh_.property(fvset_, fvector[i])).push_back(*fvset_it);
						mesh_.property(vf_, *fvset_it) = fvector[i];
					}
					is_face_located = true; 
					n_vertices_relocated++;//���smooth para�����¶�λ�Ķ�����.
				}
			} // ѭ����10�����ܵ���֮��û�к��ʵĻ��ͳ�����.//����Ҳû�й�ϵ, ��Ϊ�������ǵ����������������ʱ��Ŵ�ԭ������ɾȥ��¼, 
			if (is_face_located == false) { 
				std::cerr << "Error:tri smooth para face relocated," << old_uv << ", " << mesh_.property(vuv_, *fvset_it) << "; " << *fvset_it << ".\n"; 
				//mesh_.property(vuv_, *fvset_it) = old_uv;
				//test_vertex_ = *fvset_it;
			}// ����Ҳû��ɾȥ����ԭ�����ϵļ�¼.
		} // end of if (inside_valid == true). ������������������f_it�ϵ�һ�����relocate.

	} // end of for. ������������������f_it�ϵ�ÿһ����*fvset_it��ƽ����relocate.
}

void DecimationModel::local_smooth_parameterization_triangledomain() {
	std::cout << "���뺯�� local_smooth_parameterization_triangledomain(), before smooth parameter.\n";
	int n_vertices_base_mesh = 0;
	int n_vertices_deleted = 0, n_vertices_paraed = 0;
	for (TriMesh::VertexIter v_it(mesh_.vertices_begin()), v_end(mesh_.vertices_end()); v_it != v_end; ++v_it) {
		if (mesh_.status(v_it.handle()).deleted()) { 
			n_vertices_deleted++;
		} else n_vertices_base_mesh++;
		// std::cout << mesh_.property(vbc_, v_it.handle()) << "\t"; //�ܳɹ����, ��ʾ���custom properties����.
	}
	for (TriMesh::FaceIter f_it(mesh_.faces_begin()), f_end(mesh_.faces_end()); f_it != f_end; ++f_it) {
		if (mesh_.status(f_it.handle()).deleted() == false) {
			n_vertices_paraed += (mesh_.property(fvset_, f_it.handle())).size();//��������������ϵĶ�����.
		}
	}// ����������Ǳ�deleted�Ķ������ͱ��������Ķ�����Ŀһ��.
	std::cout << "para check: V " << mesh_.n_vertices() << "= " << n_vertices_base_mesh << " + " << n_vertices_deleted << "(=" << n_vertices_paraed << ").\n";

	// -------------
	int n_edges_nondeleted(0), n_edges_deleted(0), count_crosspatches(0), n_edges_original(0);
	for (TriMesh::EdgeIter e_it(mesh_.edges_begin()), e_end(mesh_.edges_end()); e_it != e_end; ++ e_it) {
		if (mesh_.status(e_it).deleted() == true) {
			n_edges_deleted++;
			TriMesh::VHandle v0 = mesh2_.to_vertex_handle(mesh2_.halfedge_handle(e_it, 0));
			TriMesh::VHandle v1 = mesh2_.to_vertex_handle(mesh2_.halfedge_handle(e_it, 1));
			if (mesh_.property(vf_, v0) != mesh_.property(vf_, v1)) {
				count_crosspatches++;//���ٸ��߲��ǲ�������ͬһ�����ϵ�.
			}
		} else { //�����û�б�ɾȥ,���ǿ����������˵㶼�ı���,Ҳ����˵����ڶ����֮ǰ�Ǹ�ԭʼ����, ���ڼ򻯹��������ɵ�.
			n_edges_nondeleted++;
			TriMesh::VHandle v0 = mesh_.to_vertex_handle(mesh_.halfedge_handle(e_it, 0));
			TriMesh::VHandle v1 = mesh_.to_vertex_handle(mesh_.halfedge_handle(e_it, 1));
			if (((v0 == mesh2_.to_vertex_handle(mesh2_.halfedge_handle(e_it, 0))) && (v1 == mesh2_.to_vertex_handle(mesh2_.halfedge_handle(e_it, 1))))
				|| ((v1 == mesh2_.to_vertex_handle(mesh2_.halfedge_handle(e_it, 0))) && (v0 == mesh2_.to_vertex_handle(mesh2_.halfedge_handle(e_it, 1))))) {
					// ������ڼ򻯵Ĺ��̵���û�з����仯,�������˵㻹�Ǽ�֮���������.
					n_edges_original ++;//�ж��ٸ���û�б仯
			}
		}
	}
	std::cout << "para check: E " << mesh_.n_edges() << "= " << n_edges_nondeleted << " + " << n_edges_deleted 
		<< "; crop " << count_crosspatches << ", ori " << n_edges_original << ".\n";

	// ---------//�������ƽ�������� ---------------------
	std::cout << "���濪ʼһ��" << n_local_sp_ << "�εľֲ�ƽ��.\n";
	for ( int i_iterance = 0; i_iterance < n_local_sp_; ++i_iterance) { // begin to locally smooth parameterize, n_local_sp_ times
		std::cout << "Begin the " << i_iterance << " time to locally smooth parameterization.\n";
		int n_vertices_relocated = 0; // ��¼�����smooth parameterization���ж��ٸ����㱻���¶�λ��.
		int n_vertices_changeface = 0;// �ж��ٸ������ھֲ�ƽ��ʱ���Ƶ������ƽ����
		for (TriMesh::FaceIter f_it(mesh_.faces_begin()), f_end(mesh_.faces_end()); f_it != f_end; ++f_it) {
			if (mesh_.status(f_it.handle()).deleted() == false) { //ѭ��ÿһ��û�б�deleted����,Ҳ����base complex�ϵ���
				local_smooth_parameterization_triangledomain(f_it.handle(), n_vertices_relocated, n_vertices_changeface);
			} // end of if.�ж�f_it.handle()�Ƿ�deleted, Ҳ���Ƿ���base complex�ϵ�һ����			
		} // end of for.��������ѭ��һ��, ѭ��ʱ������ϵ����в������ĵ����һ��smooth
		std::cout << "End the " << i_iterance << " time to locally smooth parameterization. " 
			<< (n_vertices_relocated * 100.0/ n_vertices_paraed) << "% vertices relocated, " << n_vertices_changeface*100/n_vertices_paraed << "% change face.\n";

	} // end of n_local_sp_ times local smooth parameterization.
	n_vertices_paraed = 0;
	for (TriMesh::FaceIter f_it(mesh_.faces_begin()), f_end(mesh_.faces_end()); f_it != f_end; ++f_it) {
		if (mesh_.status(f_it.handle()).deleted() == false) {
			n_vertices_paraed += (mesh_.property(fvset_, f_it.handle())).size();//��������������ϵĶ�����.
		}
	}
	std::cout << "para check: n_vertices_parameterized " << n_vertices_paraed << "\n";
	std::cout << "�뿪���� local_smooth_parameterization_triangledomain().\n";
}
void DecimationModel::local_smooth_parameterization_circledomain(TriMesh::VertexHandle _vh) {
	TriMesh::VHandle vh = _vh;
	if (mesh_.status(vh).deleted() == false && mesh_.is_boundary(vh) == false) 
	{
		std::vector<TriMesh::HHandle> loop_h; 
		std::vector<TriMesh::VHandle> loop; 
		std::vector<TriMesh::FHandle> loop_f;
		std::vector<TriMesh::EHandle> loop_nh;
		for (TriMesh::VOHIter voh_it(mesh_, vh); voh_it; ++voh_it) {
			loop_h.push_back(voh_it.handle());
			loop_nh.push_back(mesh_.edge_handle(mesh_.next_halfedge_handle(voh_it)));
			loop.push_back(mesh_.to_vertex_handle(voh_it)); 
			loop_f.push_back(mesh_.face_handle(voh_it.handle()));
		}

		unsigned int n = loop_h.size();//�ж��ٸ��߽綥��
		unsigned int i = 0;  
		double length = 0.0, rou_angle = 0.0;  
		std::vector<double> vec_angle;
		for (i=0 ; i<n; ++i) {//����ܵĳ���, �ܳ� 
			length += (OpenMesh::vector_cast<OpenMesh::Vec3d>(mesh_.point(loop[i])) - OpenMesh::vector_cast<OpenMesh::Vec3d>(mesh_.point(loop[(i+1)%n]))).norm();

			OpenMesh::Vec3d d1 = (OpenMesh::vector_cast<OpenMesh::Vec3d>(mesh_.point(loop[i])) - OpenMesh::vector_cast<OpenMesh::Vec3d>(mesh_.point(vh))).normalize();
			OpenMesh::Vec3d d2 = (OpenMesh::vector_cast<OpenMesh::Vec3d>(mesh_.point(loop[(i+1)%n])) - OpenMesh::vector_cast<OpenMesh::Vec3d>(mesh_.point(vh))).normalize();
			rou_angle += acos(dot(d1, d2));
			vec_angle.push_back(acos(dot(d1, d2)));
		}  //std::cout << "angle " << rou_angle << ".\n";//std::cout << "len: " << length << std::endl; //for test	
		if (vec_angle.size() != n) { std::cout << "Error: ������һ��.\n"; }

		/* // ����һ, �õ�λԲ��������
		double l = 0.0, angle = 0.0;
		double sum_wij = 0.0;
		OpenMesh::Vec2d sum_wh(0, 0, 0);
		for (i = 0; i < n; ++i) { //����vi��one-ring�ϵ���n���ڽӵ�// fix the boundary/one-ring vertices, 				
		angle = l / length * (2.0*M_PI); //std::cout << "ang: " << angle << std::endl;
		mesh_.property(vuv_, loop[i]) = OpenMesh::Vec2d(0.5*cos(angle)+0.5, 0.5*sin(angle)+0.5);
		l += (mesh_.point(loop[i]) - mesh_.point(loop[(i+1)%n])).norm();

		// ���vi��vj�����ϵ�ϵ��, ע������ʹ�õ���mean value coordiante, ������harmonic map�е�cotangent weight
		TriMesh::VertexHandle v0 = vh; // = mesh_.from_vertex_handle(loop_h[i]);
		TriMesh::VertexHandle v1 = loop[i];// =mesh_.to_vertex_handle(loop_h[i]);
		TriMesh::VertexHandle v2 = loop[(i+1)%n];//������������ȷ��ǰ���ǷǱ߽����
		TriMesh::VertexHandle v3 = loop[(i+n-1)%n];
		OpenMesh::Vec3d v1v0 = mesh_.point(v1) - mesh_.point(v0);	// ֮ǰ������ֵ�-1
		OpenMesh::Vec3d v2v0 = mesh_.point(v2) - mesh_.point(v0);
		OpenMesh::Vec3d v3v0 = mesh_.point(v3) - mesh_.point(v0);
		double v1v0v2_angle = acos(dot(v1v0, v2v0) / (v1v0.norm() * v2v0.norm()));//ʸ��v1v0��v2v0�ļн�.
		double v3v0v1_angle = acos(dot(v3v0, v1v0) / (v3v0.norm() * v1v0.norm()));//ʸ��v3v0��v1v0�ļн�.
		double wij = (tan(v1v0v2_angle / 2) + tan(v3v0v1_angle / 2)) / v1v0.norm();
		if (wij < 0) { std::cout << "Error: wij is negative.\n"; return; }

		sum_wij += wij;
		sum_wh += wij * mesh_.property(vuv_, loop[i]);
		}
		mesh_.property(vuv_, vh) = sum_wh * (1.0/sum_wij);//�����vi�����ڵ�λԲ�ϵĲ���������.
		*///std::cout << mesh_.property(vuv_, vh) << ", " << n << ":" ;// for test.
		// ���Զ�, ѹƽ, �Ƕ�����, ���Ǻ���Ľ���Ч�������԰�, ���뿴��ȥ��һ��
		double angle_scale_ratio = 2 * M_PI / rou_angle; //���ű���
		double temp_sum_angle = 0.0, l = 0;
		for (i = 0; i < n; ++i) {
			temp_sum_angle += (vec_angle[i] * angle_scale_ratio);
			l = (OpenMesh::vector_cast<OpenMesh::Vec3d>(mesh_.point(loop[(i+1)%n])) - OpenMesh::vector_cast<OpenMesh::Vec3d>(mesh_.point(vh))).norm(); 
			l *= angle_scale_ratio;
			mesh_.property(vuv_, loop[(i+1)%n]) = OpenMesh::Vec2d(l*cos(temp_sum_angle), l*sin(temp_sum_angle));
		}
		mesh_.property(vuv_, vh) = OpenMesh::Vec2d(0, 0);

		for (i = 0; i < n; ++i)  // �Զ���vh��һ���������ÿһ����loop_f[i].
		{	//std::cout << "In face " << loop_f[i] << ", ";
			std::vector<TriMesh::VHandle> fvset = mesh_.property(fvset_, loop_f[i]);
			for (std::vector<TriMesh::VHandle>::iterator fvset_it = fvset.begin(), fvset_end = fvset.end(); 
				fvset_it != fvset_end; ++fvset_it) //  �Բ���������loop_f[i]�ϵ�ÿһ������.
			{ // TriMesh::VertexHandle *fvset_it�ǲ����Ķ���.
				if (mesh2_.is_boundary(*fvset_it) == true) continue; //���������������㲻����ƽ��.
				if (mesh_.property(vp_type_, *fvset_it) == DGP::CREASE_VFT) continue; 

				bool one_ring_vertex_is_valid = true; //�ȼ��趥��*fvset_it��һ���򶥵㶼���������Բ��, �������Ч.
				double sum_wij = 0.0; OpenMesh::Vec2d sum_wh(0, 0, 0);
				for (TriMesh::VOHIter voh_it(mesh2_, *fvset_it); voh_it; ++voh_it) 
				{ 
					TriMesh::VHandle vj = mesh2_.to_vertex_handle(voh_it.handle());
					std::vector<TriMesh::FHandle>::iterator result = find(loop_f.begin(), loop_f.end(), mesh_.property(vf_, vj)); 
					if (result != loop_f.end()) { //�������vj�Ƿ����������Բ��?
						mesh_.property(vuv_, vj) = mesh_.property(vuv_, mesh_.property(vf0_, vj)) * (mesh_.property(vbc_, vj))[0]
						+ mesh_.property(vuv_, mesh_.property(vf1_, vj)) * (mesh_.property(vbc_, vj))[1] 
						+ mesh_.property(vuv_, mesh_.property(vf2_, vj)) * (mesh_.property(vbc_, vj))[2];

						sum_wij += mesh2_.property(hmvc_, voh_it.handle());
						sum_wh += mesh2_.property(hmvc_, voh_it.handle()) * mesh_.property(vuv_, vj);
					} else {/*
							if (mesh_.property(vp_type_, vj) == DGP::CREASE_VFT) {
							std::vector<TriMesh::EHandle>::iterator result = find(loop_nh.begin(), loop_nh.end(), mesh_.property(vp_feature_edge_, vj));
							if (result != loop_nh.end()) {
							mesh_.property(vuv_, vj) = mesh_.property(vuv_, mesh_.property(vf0_, vj)) * (mesh_.property(vbc_, vj))[0]
							+ mesh_.property(vuv_, mesh_.property(vf2_, vj)) * (mesh_.property(vbc_, vj))[1]//+ OpenMesh::Vec2d(0, 0)
							+ mesh_.property(vuv_, mesh_.property(vf2_, vj)) * (mesh_.property(vbc_, vj))[2];

							sum_wij += mesh2_.property(hmvc_, voh_it.handle());
							sum_wh += mesh2_.property(hmvc_, voh_it.handle()) * mesh_.property(vuv_, vj);
							} else { one_ring_vertex_is_valid = false; } 
							} else */one_ring_vertex_is_valid = false; 
					}
				} // end of for ����*fvset_it��һ���򶥵�.
				//std::cout << one_ring_vertex_is_valid << "; ";
				if (one_ring_vertex_is_valid == true) {
					mesh_.property(vuv_, *fvset_it) = sum_wh * (1.0 / sum_wij); //����㱻ƽ��֮����²�������.

					// ������*fvset_it�ҳ�������������һ������.
					bool relocated = false;
					for (int ii = 0; ii < n && relocated == false; ++ii) { // loop_f[ii]
						TriMesh::VHandle vx = vh, vy = loop[ii], vz = mesh_.to_vertex_handle(mesh_.next_halfedge_handle(loop_h[ii])); 
						OpenMesh::Vec3d bc 
							= DGP::calc_barycentric_coordinates(mesh_.property(vuv_, vx), mesh_.property(vuv_, vy), mesh_.property(vuv_, vz), mesh_.property(vuv_, *fvset_it));

						if (fabs(bc[0]) < 1.0e-10) bc[0] = 0; if (fabs(bc[1]) < 1.0e-10) bc[1] = 0; if (fabs(bc[2]) < 1.0e-10) bc[2] = 0;
						if (!(bc[0]<0) && !(bc[1]<0) && !(bc[2]<0)) { //Ҫ�����bc��������������>=0
							//std::cout << "BC: " << bc << std::endl; // for test
							// �������о���ÿһ�β�������Ӧ�û�õĻ��Ǹ��ĵ���Ϣ.
							mesh_.property(vf0_, *fvset_it) = vx; mesh_.property(vf1_, *fvset_it) = vy; mesh_.property(vf2_, *fvset_it) = vz;
							mesh_.property(vbc_, *fvset_it) = bc;//std::cout << "a ";
							if (loop_f[i] != loop_f[ii]) { //mesh_.property(vf_, *fvset_it) �͵���loop_f[i].
								// ɾȥ����ԭ���Ǹ���(loop_f[i])�ϵļ�¼
								std::vector<TriMesh::VHandle>::iterator to_be_earsed 
									= remove((mesh_.property(fvset_, loop_f[i])).begin(), (mesh_.property(fvset_, loop_f[i])).end(), TriMesh::VHandle((*fvset_it).idx()));
								(mesh_.property(fvset_, loop_f[i])).erase(to_be_earsed, (mesh_.property(fvset_, loop_f[i])).end());
								//std::cout << "b ";
								(mesh_.property(fvset_, loop_f[ii])).push_back(*fvset_it); //std::cout << "bb ";
								mesh_.property(vf_, *fvset_it) = loop_f[ii];//std::cout << "c ";
							}
							relocated = true; 
							//n_vertices_relocated++;//���smooth para�����¶�λ�Ķ�����.
						} //std::cout << "d "; 
					} //std::cout << " ok. ";
					if (relocated == false ) { std::cout << "Error: relocated, it seems don't matter.\n"; }
				} // end of if. �������*fvset_it��һ��������Ч��Բ����.

			} // end of for.�Բ���������loop_f[i]�ϵ�ÿһ������*fvset_it.
		} // end of for.�Զ���vh��һ���������ÿһ����loop_f[i].
	} // end of if ����vh�����ڱ߽��� 
}
void DecimationModel::local_smooth_parameterization_circledomain() {
	std::cout << "���뺯�� local_smooth_parameterization_circledomain.\n";
	for (TriMesh::VIter v_it(mesh_.vertices_begin()), v_end(mesh_.vertices_end()); v_it != v_end; ++v_it) 
	{ //��v_itΪԲ�Ľ�������
		local_smooth_parameterization_circledomain(v_it.handle());		
	} // end of for.�Ի������ϵ�ÿһ������.
	std::cout << "�뿪���� local_smooth_parameterization_circledomain.\n";
} // end of function local_smooth_parameterization_circledomain().
void DecimationModel::local_smooth_parameterization_process() {
	if (do_parameterization_)  {	
		if (do_smooth_) {//////////ƽ�����е������������mesh_.garbage_collection()֮ǰ,�����޷�ʹ����ǰ��Handle.
			//local_smooth_parameterization_triangledomain();
			local_smooth_parameterization_circledomain();

			for (TriMesh::VertexIter v_it(mesh_.vertices_begin()), v_end(mesh_.vertices_end()); v_it != v_end; ++v_it) {
				if (mesh_.status(v_it.handle()).deleted()) { 

					OpenMesh::Vec3d bc = mesh_.property(vbc_, v_it.handle());
					if (!((bc[0] > -1.0e-10) && (bc[1] > -1.0e-10) && (bc[2] > -1.0e-10))) { 
						std::cerr << "Error: para check, bc negative." << bc << "\n"; 
					} // for test, only assert

					TriMesh::Point vfp0 = mesh_.point(mesh_.property(vf0_, v_it.handle())); 
					if (mesh_.status(mesh_.property(vf0_, v_it.handle())).deleted() == true) { std::cerr << "Error: para check, vf0.\n"; }
					TriMesh::Point vfp1 = mesh_.point(mesh_.property(vf1_, v_it.handle()));
					if (mesh_.status(mesh_.property(vf1_, v_it.handle())).deleted() == true) { std::cerr << "Error: para check, vf1.\n"; }
					TriMesh::Point vfp2 = mesh_.point(mesh_.property(vf2_, v_it.handle()));
					if (mesh_.status(mesh_.property(vf2_, v_it.handle())).deleted() == true) { std::cerr << "Error: para check, vf2.\n"; }

					TriMesh::Point vp = vfp0 * bc[0] + vfp1 * bc[1] + vfp2 * bc[2];//��deletedҲ���Ǳ�parameterized�Ķ����3d����.

					mesh_collapsed_.set_point(v_it, vp);////
				} else {
					mesh_collapsed_.set_point(v_it, mesh_.point(v_it.handle()));//û�б�ɾȥ�Ķ�����ǻ�����Ķ���,����ֱ��ȡ������Ϳ�����.
				}
			} // end of for
		} // end of if (do_smooth_)
	} // end of if ()
}