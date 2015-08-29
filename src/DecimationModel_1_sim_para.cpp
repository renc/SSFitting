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
	// 这里只添加mesh_.add_property().
	// used by the whole 
	mesh_.add_property(vp_type_);//顶点的类型

	// only used by the simplification
	mesh_.add_property(vp_quadric_);//顶点的Quadric
	mesh_.add_property(vp_cost_);//顶点在队列中的优先级
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
	mesh_.remove_property(vp_quadric_);//顶点的Quadric
	mesh_.remove_property(vp_cost_);//顶点在队列中的优先级
	mesh_.remove_property(vp_target_);

	mesh_.remove_property(heovh_);
}

//-----------------------------------------------------------------------------
//被enqueue_vertex函数调用
// calculate and return priority: the smaller the better. 决定顶点次序的.
double DecimationModel::priority_qem(TriMesh::HalfedgeHandle _heh)
{	// use quadrics to estimate approximation error, approximation cost.
	// (from, to) -> to
	TriMesh::VertexHandle from(mesh_.from_vertex_handle(_heh));
	TriMesh::VertexHandle to(mesh_.to_vertex_handle(_heh));

	DGP::Quadricd q = mesh_.property(vp_quadric_, from);  
	q += mesh_.property(vp_quadric_, to);//q现在是这个半边对应的点对的Quadric之和

	if (vertex_optimal_placement == true) { //收缩到哪一个目标点就用哪一个点来计算代价.
		double newvp0 = 0.0, newvp1 = 0.0, newvp2 = 0.0;
		if (q.optimal_placement(newvp0, newvp1, newvp2)) {//
			return q(TriMesh::Point(newvp0, newvp1, newvp2));
		}
	}

	// 尝试考虑生成之后的三角形的quality.
	// collect faces 半边的对应的面
	std::vector<TriMesh::HalfedgeHandle> loop_h;
	TriMesh::HHandle hh = _heh;
	do {
		loop_h.push_back(hh);
		hh = mesh_.opposite_halfedge_handle(mesh_.prev_halfedge_handle(hh));
	} while(hh != _heh);
	unsigned int n = loop_h.size();
	TriMesh::FaceHandle fl = mesh_.face_handle(_heh);
	TriMesh::FaceHandle fr = mesh_.face_handle(mesh_.opposite_halfedge_handle(_heh));
	// backup point positions 这两个变量的值不被修改	
	TriMesh::Point p0 = mesh_.point(from);
	TriMesh::Point p1 = mesh_.point(to);
	mesh_.set_point(from, p1);

	double aspect_ratio(0); int n_count(0);// (1) 三角形形状最长和最短边的比率
	const double av_ang = 60*M_PI/180;
	double angle_deviation = 60*M_PI/180;//(2)平均角度与60读的差距
	double tri_area = 0.0; // (3)trianle area
	unsigned int fvset_size = 1;//(4)参数化到这个面上的顶点数量,注意不要初始化为0,否则下面可能乘以0, 这个参数好像并不现实的.
	double area_ratio = 1;//(5) 三角形面积比率
	double valence_deviation = 0.0;//(6)与阶为6的差距
	for (unsigned int i = 1; i <= n-2; ++i) { // 不算fl和fr这两个三角形的.
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
	//else valence_deviation = 6 / valence_deviation;//上面对太大的valence有一个放大代价的作用
	if (valence_deviation != 6) valence_deviation *= 100;
	//std::cout << aspect_ratio << ", " << q(mesh_.point(to)) << ", ";
	// undo simulation changes
	mesh_.set_point(from, p0);

	double weight = 1//aspect_ratio;//效果不好

		// fvset_size * (angle_deviation + aspect_ratio);
		////* valence_deviation 
		//* area_ratio//
		//* aspect_ratio 
		//* angle_deviation //一直都除去的.
		;
	//fvset_size * angle_deviation * aspect_ratio;
	//(angle_deviation + aspect_ratio);//fandisk.off有些顶点出错
	//angle_deviation;//fandisk.off效果不好

	return weight  * q(mesh_.point(to));
	//if (fabs(q(mesh_.point(to))) < 1.e-6) { return aspect_ratio; }
	//else return aspect_ratio * q(mesh_.point(to)); //fandisk.off的效果也不好.
	//return 0.5*aspect_ratio + 0.5*q(mesh_.point(to));	// +好像效果不好.
	//return q(mesh_.point(to));// evaluate quadric Q at vector v: v*Q*v, 求出error cost
	// 这个好像需要是点对收缩而成的新点来作为衡量代价的v的,而非from或to, 则q(NewPoint(from and to)), 对, 本该是这样的.
	// 但是这里使用的是halfedge collapse, 和QEM中原始使用的vertex pair collpase有点区别.
}

//-----------------------------------------------------------------------------
//被decimate函数调用
//求出这个顶点(的某一条出半边)收缩时候的代价, 这作为在queque中的排序依据,
//将这个顶点最小的代价的那半边存起来以备后用.
void DecimationModel::enqueue_vertex(TriMesh::VertexHandle _vh)
{
	float	 min_prio(FLT_MAX);//float值越大也就是代价越大, 优先级越小
	TriMesh::HalfedgeHandle  min_hh;//与这个顶点相关的最小代价的半边, 

	// find best out-going halfedge, _vh是作为四周的出半边的起点的
	for (TriMesh::VOHIter vh_it(mesh_, _vh); vh_it; ++vh_it) // 这个顶点的四周的出半边
	{
		if (DGP::QEMDecimationT<TriMesh>::is_collapse_legal(mesh_, vh_it, vp_type_))
		{	//这个出半边, 也就是这个出半边的两个顶点(一点对)可以收缩
			float prio = priority_qem(vh_it);//求一个点对(由这点对组成的半边表示)收缩时候的代价
			if (prio != -1.0 && prio < min_prio)
			{
				min_prio = prio;
				min_hh   = vh_it.handle();
			}
		}
	} //获得与这个顶点相关的所有点对中哪个点对的代价最小(优先级最大)记为min_prio, 这里点对用半边表示
	//其实在求点对的代价,也就是 float priority(HalfedgeHandle )函数中只是求这个半边的to vertex(target vertex)的代价

	// update queue
	if (mesh_.property(vp_cost_, _vh) != -1.0) 
	{
		queue.erase(_vh);
		mesh_.property(vp_cost_, _vh) = -1.0;
	}

	if (min_hh.is_valid()) 
	{
		mesh_.property(vp_cost_, _vh) = min_prio;//这个东西用于在queue中排序的.
		mesh_.property(vp_target_, _vh)   = min_hh;
		queue.insert(_vh);//插入这个顶点_vh, 需要使用到mesh_.property(vp_cost_, _vh)来决定其在队列中的位置
		//这样的话从队列里面获得最小代价的点_vh,在mesh_.property(vp_target_, _vh)中获得其对应的min_hh
	}
}

void DecimationModel::initial_parameterization(TriMesh::HalfedgeHandle _hh) {
	// (from, to)->to
	TriMesh::VertexHandle from = mesh_.from_vertex_handle(_hh);
	TriMesh::VertexHandle to = mesh_.to_vertex_handle(_hh);
	TriMesh::HalfedgeHandle o = mesh_.opposite_halfedge_handle(_hh);
	//std::cout << "Para (" << from << ", " << to <<")\n";

	// 下面分三种情况来讨论, 第一种半边_hh及其对半边都不是边界边,也就是说from顶点是interior vertex,
	// 第二种情况是_hh是边界边, from顶点是边界顶点, 第三种情况是_hh的对半边是边界边, 这时from顶点也是边界点.
	if (!mesh_.is_boundary(_hh) && !mesh_.is_boundary(o)) { //非边界情况, from是interior vertex
		// 第一种case
		//std::cout << "c1\n"; // for test 
		// (from, to)->to, 将from和其一邻域上的点参数化到一个单位圆中, 做Mean Value Mapping. 
		// 将from顶点看做是vi顶点, 而from顶点的one-ring邻域都是vj顶点

		std::vector<TriMesh::VertexHandle> loop; //用于循序地记录vi顶点(这里也就是from顶点)的one-ring邻域里面的顶点
		std::vector<TriMesh::HalfedgeHandle> loop_h;//用于记录vi顶点(这里也就是from顶点)的one-ring邻域里面的出半边
		std::vector<TriMesh::FHandle> loop_f;
		TriMesh::HalfedgeHandle h = _hh;
		do {	// 这里不使用VertexOHalfedgeIter和VertexVertexIter的原因是我要保证loop[j]的顺序.loop[0]=to
			loop_h.push_back(h);
			loop_f.push_back(mesh_.face_handle(h));
			loop.push_back(mesh_.to_vertex_handle(h));//std::cout << mesh_.to_vertex_handle(h) << "\t";	 ////	for test	
			h = mesh_.opposite_halfedge_handle(mesh_.prev_halfedge_handle(h));
		} while(h != _hh);
		if (loop.size() != loop_h.size()) { std::cerr << "Error: collect one-ring vertices and halfedges.\n"; };

		unsigned int n = loop_h.size();//有多少个边界顶点
		unsigned int i = 0;
		double length = 0.0, rou_angle = 0.0;//
		std::vector<double> vec_angle;
		for (i=0 ; i<n; ++i) {//求出总的长度(周长), 以及角度之和. 
			length += (OpenMesh::vector_cast<OpenMesh::Vec3d>(mesh_.point(loop[i])) - OpenMesh::vector_cast<OpenMesh::Vec3d>(mesh_.point(loop[(i+1)%n]))).norm();

			OpenMesh::Vec3d d1 = (OpenMesh::vector_cast<OpenMesh::Vec3d>(mesh_.point(loop[i])) - OpenMesh::vector_cast<OpenMesh::Vec3d>(mesh_.point(from))).normalize();
			OpenMesh::Vec3d d2 = (OpenMesh::vector_cast<OpenMesh::Vec3d>(mesh_.point(loop[(i+1)%n])) - OpenMesh::vector_cast<OpenMesh::Vec3d>(mesh_.point(from))).normalize();
			double angle = atan2(cross(d1, d2).norm(), dot(d1, d2)); //acos(dot(d1, d2));
			rou_angle += angle;
			vec_angle.push_back(angle);
		}  //std::cout << "angle " << rou_angle << ".\n";//std::cout << "len: " << length << std::endl; //for test	
		if (vec_angle.size() != n) { std::cout << "Error: 容量不一致.\n"; } // 1000

		bool is_crease_edge = false; //先判断要收缩的_hh是不是crease edge, 也就是from顶点是不是crease vertex, 
		if (mesh_.property(vp_type_, from) == DGP::CREASE_VFT) { // if (mesh_.status(mesh_.edge_handle(_hh)).feature() == true)是一样的.
			is_crease_edge = true;
		}
		int feature_vertex_index = -1;// 当is_crease_edge == true;标记loop[]数组中哪一个是特征点.
		if (is_crease_edge == true) {//找出feature_vertex_index等于多少.
			for (i = 2; i < n-1; ++i) //另一个crease edge只允许在loop_h[2]到loop_h[n-2]这范围之内
				if (mesh_.status(mesh_.edge_handle(loop_h[i])).feature() == true) 
					feature_vertex_index = i;
			if (feature_vertex_index == -1) std::cerr << "Error: initial_parameterization(): feature vertex index.\n";
		} // 只要is_crease_edge == true, 那么
		double bcv, bt;//当is_crease_edge == true时候才有意义, 是重心坐标的两小部分, 实际上是比例.

		//---------------------------------------------------------------------------------------------------
		/*// 法一 单位圆
		double l = 0.0, angle = 0.0;
		// to decide the coordiante of vertex from.
		if (is_crease_edge == false) // using the mean value coordiante
		{
		//std::vector<OpenMesh::Vec2d > h_vj(loop_h.size());//vj顶点的参数化坐标
		std::vector<double> wij(loop_h.size());//这是顶点vi和vj的系数
		double sum_wij = 0.0;
		OpenMesh::Vec2d sum_wh(0, 0, 0);
		for (i = 0; i < n; ++i) { //顶点vi在one-ring上的有n个邻接点
		// fix the boundary/one-ring vertices, 
		angle = l / length * (2.0*M_PI); //std::cout << "ang: " << angle << std::endl;
		mesh_.property(vuv_, loop[i]) = OpenMesh::Vec2d(0.5*cos(angle)+0.5, 0.5*sin(angle)+0.5);
		//这个单位圆是以(0.5, 0.5)为圆心, 半径是0.5. 圆心要是定在(0, 0)的话, 例如本来是(0, 0.5)的坐标会机器误差为(e, 0.5),e是一个很小很小的值.
		l += (mesh_.point(loop[i]) - mesh_.point(loop[(i+1)%n])).norm();

		// 求出vi和vj这半边上的系数, 注意这里使用的是mean value coordiante, 而不是harmonic map中的cotangent weight
		TriMesh::VertexHandle v0 = from; // = mesh_.from_vertex_handle(loop_h[i]);
		TriMesh::VertexHandle v1 = loop[i];// =mesh_.to_vertex_handle(loop_h[i]);
		TriMesh::VertexHandle v2 = loop[(i+1)%n];//= mesh_.to_vertex_handle(mesh_.next_halfedge_handle(loop_h[i]));//这以下两句正确的前提是非边界情况
		TriMesh::VertexHandle v3 = loop[(i+n-1)%n];//mesh_.to_vertex_handle(mesh_.next_halfedge_handle(mesh_.opposite_halfedge_handle(loop_h[i])));
		OpenMesh::Vec3d v1v0 = mesh_.point(v1) - mesh_.point(v0);	// 之前这里出现的-1
		OpenMesh::Vec3d v2v0 = mesh_.point(v2) - mesh_.point(v0);
		OpenMesh::Vec3d v3v0 = mesh_.point(v3) - mesh_.point(v0);
		double v1v0v2_angle = acos(dot(v1v0, v2v0) / (v1v0.norm() * v2v0.norm()));//矢量v1v0和v2v0的夹角.
		double v3v0v1_angle = acos(dot(v3v0, v1v0) / (v3v0.norm() * v1v0.norm()));//矢量v3v0和v1v0的夹角.
		wij[i] = (tan(v1v0v2_angle / 2) + tan(v3v0v1_angle / 2)) / v1v0.norm();
		if (wij[i] < 0) { std::cout << "Error: wij is negative.\n"; return; }

		sum_wij += wij[i];
		sum_wh += wij[i] * mesh_.property(vuv_, loop[i]);
		}
		mesh_.property(vuv_, from) = sum_wh * (1.0/sum_wij);//这就是vi顶点在单位圆上的参数化坐标.
		// 至此求出了from顶点在其一邻域内的参数化uv坐标.
		} else {
		for (i = 0; i < n; ++i) { //顶点vi在one-ring上的有n个邻接点
		// fix the boundary/one-ring vertices, 
		angle = l / length * (2.0*M_PI); //std::cout << "ang: " << angle << std::endl;
		mesh_.property(vuv_, loop[i]) = OpenMesh::Vec2d(0.5*cos(angle)+0.5, 0.5*sin(angle)+0.5);
		l += (mesh_.point(loop[i]) - mesh_.point(loop[(i+1)%n])).norm();
		} 
		double len_to_from = (double)(mesh_.point(to) - mesh_.point(from)).norm();
		double len_creasevertex_from = (double)(mesh_.point(loop[feature_vertex_index]) - mesh_.point(from)).norm();
		bcv = len_to_from / (len_to_from + len_creasevertex_from);//barycentric coordinate中相对于顶点crease vertex的分量.
		bt = 1- bcv;
		mesh_.property(vuv_, from) = (mesh_.property(vuv_, to) - mesh_.property(vuv_, loop[feature_vertex_index])) * bt + mesh_.property(vuv_, loop[feature_vertex_index]);
		}*/
		// 法二 压平
		double angle_scale_ratio = 2 * M_PI / rou_angle; //缩放比率
		double temp_sum_angle = 0.0, len = 0;
		for (int i = 0; i < n; ++i) { //一邻域点的参数化坐标.
			temp_sum_angle += (vec_angle[i] * angle_scale_ratio);
			len = (OpenMesh::vector_cast<OpenMesh::Vec3d>(mesh_.point(loop[(i+1)%n])) - OpenMesh::vector_cast<OpenMesh::Vec3d>(mesh_.point(from))).norm(); 
			len *= angle_scale_ratio;
			mesh_.property(vuv_, loop[(i+1)%n]) = OpenMesh::Vec2d(len*cos(temp_sum_angle), len*sin(temp_sum_angle));
		} 		
		if (feature_vertex_index != -1) {
			double len_to_from = (double)(mesh_.point(to) - mesh_.point(from)).norm();
			double len_creasevertex_from = (double)(mesh_.point(loop[feature_vertex_index]) - mesh_.point(from)).norm();
			bcv = len_to_from / (len_to_from + len_creasevertex_from);//barycentric coordinate中相对于顶点crease vertex的分量.
			bt = 1- bcv; //from这特征点的参数坐标直接给出, 而参数化到哪个面就后面统一决定.
			mesh_.property(vuv_, from) = mesh_.property(vuv_, to)
				+ (mesh_.property(vuv_, loop[feature_vertex_index]) - mesh_.property(vuv_, to)) * bcv;
		} else { // not crease 
			mesh_.property(vuv_, from) = OpenMesh::Vec2d(0, 0); //完成一邻域参数化 
		} 

		//---------------------------------------------------------------------------------------------------		

		//收集所有的free vertices, 这些free vertices包括了两个部分.
		//分别是from顶点, 以及from顶点一邻域的n个面上所包含的在之前的参数化过程中包含的顶点.
		// 做法(1)无论from顶点是否crease vertex都加入到free_vertics数组中来计算其落到哪一个面上,
		// 做法(2)from顶点不是crease vertex时候才加入到free_vertics数组, 因为已经可以直接给出其参数坐标信息.
		// 这里使用的是做法(2).
		std::vector<TriMesh::VertexHandle> free_vertices;		
		if (is_crease_edge == false) {
			free_vertices.push_back(from);
		}////
		for (i = 0; i < n; ++i) //循环n个面
		{	
			TriMesh::FaceHandle fh = loop_f[i];//= mesh_.face_handle(loop_h[i]);
			TriMesh::VertexHandle v0 = from;//=mesh_.from_vertex_handle(loop_h[i]);//这个面的3个顶点, 这个就是from顶点
			TriMesh::VertexHandle v1 = loop[i];//=mesh_.to_vertex_handle(loop_h[i]);
			TriMesh::VertexHandle v2 = loop[(i+1)%n];//=mesh_.to_vertex_handle(mesh_.next_halfedge_handle(loop_h[i]));

			//取出这个面所包含的在之前的参数化过程当中参数化到这个面上的顶点, 求出它们在新压平的one-ring中的参数化uv坐标
			//这些点其中有一部分是crease edge顶点.
			std::vector<TriMesh::VertexHandle> fvset = mesh_.property(fvset_, fh);
			for (std::vector<TriMesh::VertexHandle>::iterator it(fvset.begin()), it_end(fvset.end()); it != it_end; ++it) {
				//初始的时候每一个面上还没有参数化落在上面的顶点, 此时fvset为空.
				if (mesh_.status(*it).deleted() == false) { std::cerr << "Error: \n"; return; } // 这些以前已经参数话过的点应该都是deleted vertices.
				free_vertices.push_back(*it);				

				// 同时求出这些free vertices在这个一邻域里面的参数坐标.
				TriMesh::VertexHandle vf0 = mesh_.property(vf0_, *it); //这个*it顶点在以前的参数化过程中落在的面的三个顶点,
				TriMesh::VertexHandle vf1 = mesh_.property(vf1_, *it); //这三个顶点所构成的面, 刚好应该就是上面循环中的那个面. 
				TriMesh::VertexHandle vf2 = mesh_.property(vf2_, *it); //而且三个点之中有一个就是from顶点.
				if (vf0 != v0 && vf0 != v1 && vf0 != v2) { std::cerr << "Error: not match vf0.\n"; return; }// 
				if (vf1 != v0 && vf1 != v1 && vf1 != v2) { std::cerr << "Error: not match vf1.\n"; return; }// 
				if (vf2 != v0 && vf2 != v1 && vf2 != v2) { std::cerr << "Error: not match vf2.\n"; return; }//
				mesh_.property(vuv_, *it) = mesh_.property(vuv_, vf0) * (mesh_.property(vbc_, *it))[0]
				+ mesh_.property(vuv_, vf1) * (mesh_.property(vbc_, *it))[1]
				+ mesh_.property(vuv_, vf2) * (mesh_.property(vbc_, *it))[2];
			}

			(mesh_.property(fvset_, fh)).clear();//这个面上面的顶点先清空, 下面再重新看有哪些顶点参数化落在这里.
		}

		// 在initial parameterization中加入了一点平滑,但是有个问题是平滑后可能有些顶点被拉到from顶点的一邻域面之外了.
		for (std::vector<TriMesh::VertexHandle>::iterator freev_it = free_vertices.begin(), freev_end = free_vertices.end(); freev_it != freev_end; ++freev_it) 
		{ //freev_it(from顶点),不用做平滑得. 
			if (*freev_it == from) continue;
			if (mesh2_.is_boundary(*freev_it)) continue;//如果是边界上的顶点(参数化的时候也参数化到边界上的)就不要进行smooth parameterize.
			if (mesh_.property(vp_type_, *freev_it) == DGP::CREASE_VFT)  continue;					 

			bool one_ring_vertex_is_valid = true; //先假设顶点*freev_it的一邻域顶点都在这个参数圆内, 这就是有效.
			double sum_wij = 0.0; OpenMesh::Vec2d sum_wh(0, 0, 0);
			for (TriMesh::VOHIter voh_it(mesh2_, *freev_it); voh_it; ++voh_it) 
			{ 
				TriMesh::VHandle vj = mesh2_.to_vertex_handle(voh_it.handle());
				std::vector<TriMesh::FHandle>::iterator result = find(loop_f.begin(), loop_f.end(), mesh_.property(vf_, vj)); 
				if (result != loop_f.end()) { //这个顶点vj是否在这个参数圆里?
					sum_wij += mesh2_.property(hmvc_, voh_it.handle());
					sum_wh += mesh2_.property(hmvc_, voh_it.handle()) * mesh_.property(vuv_, vj);
				} else { one_ring_vertex_is_valid = false; }
			} // end of for
			if (one_ring_vertex_is_valid == true) {
				mesh_.property(vuv_, *freev_it) = sum_wh * (1.0 / sum_wij);
			}
		}/**//**/

		// 下面是定位这些free vertices是落在半边_hh收缩后形成的哪一个三角形之内, 半边_hh收缩之后一邻域的面从n减为n-2个.
		// 所有的free vertices初始化为还没有定位到三角形.
		OpenMesh::VPropHandleT<bool> is_face_located;
		mesh_.add_property(is_face_located);
		for (std::vector<TriMesh::VHandle>::size_type freev_i = 0, freev_size = free_vertices.size(); freev_i < freev_size; ++ freev_i) {
			mesh_.property(is_face_located, free_vertices[freev_i]) = false;
		}

		for (i = 1 ;i <= n-2; ++i) { // 循环n-2个面, 这n-2个面是假设(from, to)->to之后形成的
			// 判断那所有需要重新定位的free vertices是否在这个triangle(loop[0], loop[i], loop[i+1])之内
			for (std::vector<TriMesh::VertexHandle>::iterator freev_it = free_vertices.begin(), freev_end = free_vertices.end(); freev_it != freev_end; ++freev_it) {
				if (mesh_.property(is_face_located, *freev_it) == false) {

					OpenMesh::Vec2d pp = mesh_.property(vuv_, *freev_it);//这个free vertex的参数坐标
					OpenMesh::Vec3d bc = DGP::calc_barycentric_coordinates(mesh_.property(vuv_, loop[0]), mesh_.property(vuv_, loop[i]), mesh_.property(vuv_, loop[i+1]), pp);
					//std::cout << bc << std::endl;//for test

					//(0.5, -1.11022e-016, 0.5)//这里是两个特殊的重心坐标,其中的两个很小的小数-1.11022e-016和5.55112e-017都应该判为0的,
					//(0.5, 0.5, 5.55112e-017) //这样的话第一个重心坐标其实是合乎在一个三角形内部(边界上)的内部点的要求的.
					if (fabs(bc[0]) < 1.0e-10) bc[0] = 0;//原来是1.0e-6,现在改为1.0e-5, 否则像bc(-1.37927e-006 0.19194 0.808062)
					if (fabs(bc[1]) < 1.0e-10) bc[1] = 0;//就被认为是无效的. 而事实上第一项该为0, 后两项和为1.
					if (fabs(bc[2]) < 1.0e-10) bc[2] = 0;//还有(0.303709 -0.000569765 0.696861)

					//if (!(bc[0]<0) && !(bc[1]<0) && !(bc[2]<0)) { //要求这个bc的三个分量都是>=0
					if ((bc[0] >= 0) && (bc[1] >= 0) && (bc[2] >= 0)) { //要求这个bc的三个分量都是>=0

						mesh_.property(is_face_located, *freev_it) = true;

						// For interior vertex, its' barycentric coordinates is bc according to the triangle(0, i, i+1).
						//std::cout << *freev_it << ": " << loop[0] << "," << loop[i] << "," << loop[i+1] << ": " << bc << "\n"; // for test

						TriMesh::FaceHandle fh = mesh_.face_handle(mesh_.find_halfedge(loop[i], loop[i+1]));
						if (fh.is_valid() == false) std::cerr << "Error: face invalid." << std::endl;//Only for assert

						// 下面四行就是每一次参数化所应该获得的或是更改的信息.
						mesh_.property(vf_,  *freev_it) = fh;
						mesh_.property(vf0_, *freev_it) = loop[0]; mesh_.property(vf1_, *freev_it) = loop[i]; mesh_.property(vf2_, *freev_it) = loop[i+1];
						mesh_.property(vbc_, *freev_it) = bc;

						(mesh_.property(fvset_, fh)).push_back(*freev_it); 
					}
				} // 这个free vertex落在这个面上
			} // 一部分free vertex 落在这个面上			
		} // 所有的free vertices 落在这n-2个面上
		for (std::vector<TriMesh::VHandle>::size_type freev_i = 0, freev_size = free_vertices.size(); freev_i < freev_size; ++ freev_i) { 
			//这里起了判断的作用,估计可以去掉的.
			if (mesh_.property(is_face_located, free_vertices[freev_i]) == false) { //表示这个free vertex这上面的n-2个面中没有找到投影到其中哪一个上, 这是错误的.
				std::cerr << "Error: free vertex reloacted 1. " << free_vertices[freev_i] << "; " << from << ", " << to << "\n";// 
				return;
			}  
		}
		mesh_.remove_property(is_face_located);

		if (is_crease_edge == true) {
			// 首先, crease vertex类型的from顶点没有加入到free_vertcies数组中, from顶点参数化信息还没有完整.
			mesh_.property(vf_,  from) = loop_f[feature_vertex_index - 1]; // = mesh_.face_handle(loop_h[feature_vertex_index - 1]);左右两个面中默认选左的那个
			mesh_.property(vf0_, from) = loop[0]; mesh_.property(vf1_, from) = loop[feature_vertex_index - 1]; mesh_.property(vf2_, from) = loop[feature_vertex_index];
			mesh_.property(vbc_, from) = OpenMesh::Vec3d(bt, 0, bcv);

			(mesh_.property(fvset_, loop_f[feature_vertex_index-1])).push_back(from);

			mesh_.property(vp_feature_edge_, from) = mesh_.edge_handle(loop_h[feature_vertex_index]);//这个from点被参数化到那个特征边上.

			// 特征半边记录其有什么半边顺序组成, 下面的for语句应该由copy算法函数代替
			// 正面
			for (std::vector<TriMesh::HalfedgeHandle>::iterator it(mesh_.property(hep_heset_, _hh).begin()), end(mesh_.property(hep_heset_, _hh).end()); it != end; ++it)
			{ //将被收缩的半边_hh上的之前被参数化的特征边的记录都移到新的特征边上.
				mesh_.property(hep_heset_, mesh_.opposite_halfedge_handle(loop_h[feature_vertex_index])).push_back(*it);

				if (mesh2_.to_vertex_handle(*it) != to) { //将被参数化到_hh边上的点, 都更新到新的边上
					mesh_.property(vp_feature_edge_, mesh2_.to_vertex_handle(*it)) = mesh_.edge_handle(loop_h[feature_vertex_index]);
				}
			}
			// 反面. 将被参数化到_hh的反半边上的点, 这里没有更新到新的边上是因为上面其实已经达到效果了, 点是一样的.
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
		// 第二种case
		//std::cout << "c2\n"; // for test // means in case 2
		std::vector<TriMesh::VertexHandle> loop; //收集from顶点的one-ring顶点vj// 初始时候vj = to, loop[0] = to 
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
		unsigned int i, n = loop.size();//边界上有n个顶点
		TriMesh::Scalar  angle, l, length;
		double len_to_from = (double)(mesh_.point(to) - mesh_.point(from)).norm();
		double len_vn_1_from = (double)(mesh_.point(loop[n-1]) - mesh_.point(from)).norm();

		for (i = 0, length = 0.0; i < n-1; ++i) //求出总的长度
			length += (mesh_.point(loop[i]) - mesh_.point(loop[(i+1)%n])).norm();
		length += (len_to_from + len_vn_1_from);

		for (i=0, l=0.0; i<n; ++i) { //fix the boundary vertices		
			angle = l/length * (2.0*M_PI);
			mesh_.property(vuv_, loop[i]) = OpenMesh::Vec2d(0.5*cos(angle)+0.5, 0.5*sin(angle)+0.5);//这个单位圆是以(0.5, 0.5)为圆心, 半径是0.5.
			l += (mesh_.point(loop[i]) - mesh_.point(loop[(i+1)%n])).norm();
		} // 至此, from顶点的一邻域顶点已经参数化完毕.
		for (i = 0; i < n; ++i) {
			//std::cout << mesh_.property(vuv_, loop[i]) << "\n";
		}

		// 对于from顶点, 特殊 
		double bn_1 = len_to_from / (len_to_from + len_vn_1_from);
		double bt = 1 - bn_1; // = len_vn_1_from / (len_to_from + len_vn_1_from);
		OpenMesh::Vec3d bc_from(bt, 0, bn_1);
		mesh_.property(vuv_, from) = (mesh_.property(vuv_, loop[0]) - mesh_.property(vuv_, loop[n-1])) * bt + mesh_.property(vuv_, loop[n-1]);
		//std::cout << from << ": " << to<< "," << loop[n-2] << "," << loop[n-1] << ": " << bc_from << "\n"; // for test

		// 收集from顶点的One-ring邻域的面上的顶点, 这些顶点是在之前的参数化过程中参数化到落到相应的三角面上的.
		std::vector<TriMesh::VertexHandle> free_vertices;
		OpenMesh::VPropHandleT<bool> is_face_located;
		mesh_.add_property(is_face_located);
		for (i = 1; i <= n-1; ++i) { //循环from顶点的一邻域面, 共n-1个.
			TriMesh::FaceHandle fh = mesh_.face_handle(loop_h[i]);
			TriMesh::VertexHandle v0 = from;//=mesh_.from_vertex_handle(loop_h[i]);//这个面的3个顶点, 有一个就是from顶点
			TriMesh::VertexHandle v1 = loop[i-1];////另外两个是from顶点的一邻域顶点,
			TriMesh::VertexHandle v2 = loop[i];//mesh_.to_vertex_handle(loop_h[i]);//这三个顶点的参数坐标已求出了.

			std::vector<TriMesh::VertexHandle> fvset = mesh_.property(fvset_, fh);
			for (std::vector<TriMesh::VertexHandle>::iterator it(fvset.begin()), it_end(fvset.end()); it != it_end; ++it) {
				//初始的时候每一个面上还没有参数化落在上面的顶点, 此时fvset为空.
				if (mesh_.status(*it).deleted() == false) { std::cerr << "Error: \n"; return; } // 这些以前已经参数话过的点应该都是deleted vertices.
				free_vertices.push_back(*it);

				TriMesh::VertexHandle vf0 = mesh_.property(vf0_, *it); //这个*it顶点在以前的参数化过程中落在的面的三个顶点,
				TriMesh::VertexHandle vf1 = mesh_.property(vf1_, *it); //这三个顶点所构成的面, 刚好应该就是上面循环中的那个面. 
				TriMesh::VertexHandle vf2 = mesh_.property(vf2_, *it);
				if (vf0 != v0 && vf0 != v1 && vf0 != v2) { std::cerr << "Error: parameterziation vf0.\n"; return; }
				if (vf1 != v0 && vf1 != v1 && vf1 != v2) { std::cerr << "Error: parameterziation vf1.\n"; return; }
				if (vf2 != v0 && vf2 != v1 && vf2 != v2) { std::cerr << "Error: parameterziation vf2.\n"; return; }
				mesh_.property(vuv_, *it) = mesh_.property(vuv_, vf0) * (mesh_.property(vbc_, *it))[0]
				+ mesh_.property(vuv_, vf1) * (mesh_.property(vbc_, *it))[1]
				+ mesh_.property(vuv_, vf2) * (mesh_.property(vbc_, *it))[2];//求出这个自由顶点的新的参数化坐标.
			}
			(mesh_.property(fvset_, fh)).clear();//这个面上面的顶点先清空, 下面再重新看有哪些顶点参数化落在这里.
		}
		for (unsigned int freev_i = 0, freev_size = free_vertices.size(); freev_i < freev_size; ++ freev_i) {
			mesh_.property(is_face_located, free_vertices[freev_i]) = false;//所有的free vertices初始为还没有定位三角形.
		}

		for (unsigned int j = 1; j <= n-2; ++j) { //循环n-2个面, 这些面是在(from, to)->to之后去除fl和fr后形成的.
			TriMesh::VertexHandle v0 = loop[0];//=to
			TriMesh::VertexHandle vj = loop[j];
			TriMesh::VertexHandle vj1 = loop[j+1];//这三个顶点形成一个三角形, 下面判断free_vertices数组里面的顶点是否落在其中.
			TriMesh::FaceHandle fh = mesh_.face_handle(loop_h[j+1]);

			for (std::vector<TriMesh::VertexHandle>::iterator freev_it = free_vertices.begin(), freev_end = free_vertices.end(); freev_it != freev_end; ++freev_it) {
				if (mesh_.property(is_face_located, *freev_it) ==false) {

					OpenMesh::Vec2d pp = mesh_.property(vuv_, *freev_it);//这个free vertex的参数坐标
					OpenMesh::Vec3d bc = DGP::calc_barycentric_coordinates(mesh_.property(vuv_,loop[0]), mesh_.property(vuv_, loop[j]), mesh_.property(vuv_, loop[j+1]), pp);
					//std::cout << bc << std::endl;//for test
					//(0.5, -1.11022e-016, 0.5)//这里是两个特殊的重心坐标,其中的两个很小的小数-1.11022e-016和5.55112e-017都应该判为0的,
					//(0.5, 0.5, 5.55112e-017) //这样的话第一个重心坐标其实是合乎在一个三角形内部(边界上)的内部点的要求的.
					if (fabs(bc[0]) < 1.0e-10) bc[0] = 0;
					if (fabs(bc[1]) < 1.0e-10) bc[1] = 0;
					if (fabs(bc[2]) < 1.0e-10) bc[2] = 0;
					if (!(bc[0]<0) && !(bc[1]<0) && !(bc[2]<0)) { //要求这个bc的三个分量都是>=0

						// For interior vertex, its' barycentric coordinates is bc according to the triangle(v0, vj, vj1).
						//std::cout << *freev_it << ": " << loop[0] << "," << loop[j] << "," << loop[j+1] << ": " << bc << "\n"; // for test
						// 下面四行就是每一次参数化所应该获得的或是更改的信息.
						mesh_.property(vf_, *freev_it) = fh;
						mesh_.property(vf0_, *freev_it) = loop[0]; mesh_.property(vf1_, *freev_it) = loop[j]; mesh_.property(vf2_, *freev_it) = loop[j+1];
						mesh_.property(vbc_, *freev_it) = bc;

						(mesh_.property(fvset_, fh)).push_back(*freev_it);

						mesh_.property(is_face_located, *freev_it) = true;
					}
				} // 这个free vertex落在这个面上
			} // 一部分free vertex 落在这个面上			
		}
		mesh_.property(vbc_, from) = bc_from;
		mesh_.property(vf0_, from) = to; mesh_.property(vf1_, from) = loop[n-2]; mesh_.property(vf2_, from) = loop[n-1];
		mesh_.property(vf_, from) = mesh_.face_handle(loop_h[n-1]);
		(mesh_.property(fvset_, mesh_.face_handle(loop_h[n-1]))).push_back(from);

		for (unsigned int freev_i = 0, freev_size = free_vertices.size(); freev_i < freev_size; ++ freev_i) { //这里起了判断的作用,估计可以去掉的.
			if (mesh_.property(is_face_located, free_vertices[freev_i]) == false) { //表示这个free vertex这上面的n-2个面中
				std::cerr << "Error: parameterization case 2.\n";// 没有找到投影到其中哪一个上, 这是错误的.
				return;
			}  
		}
		mesh_.remove_property(is_face_located);

	} else if (!mesh_.is_boundary(_hh) && mesh_.is_boundary(o)) {
		// 第三种case
		//std::cout << "c3\n"; // for test 
		std::vector<TriMesh::VertexHandle> loop; //收集from顶点的one-ring顶点vj
		std::vector<TriMesh::HalfedgeHandle> loop_h;
		TriMesh::HalfedgeHandle h = _hh;
		TriMesh::VertexHandle vj = mesh_.to_vertex_handle(h);// 初始时候vj = to
		do {
			loop_h.push_back(h);
			loop.push_back(vj);//std::cout << "vj: " << vj << "\n";// for test
			h = mesh_.opposite_halfedge_handle(mesh_.prev_halfedge_handle(h));
			vj = mesh_.to_vertex_handle(h);
		} while(vj != to);
		assert(loop.size() == loop_h.size());

		// map boundary loop to unit circle 2D domain
		unsigned int i, n = loop.size();//边界上有n个顶点
		TriMesh::Scalar  angle, l, length;
		double len_to_from = (double)(mesh_.point(to) - mesh_.point(from)).norm();
		double len_vn_1_from = (double)(mesh_.point(loop[n-1]) - mesh_.point(from)).norm();

		for (i=0, length=0.0; i<n-1; ++i) //求出总的长度
			length += (mesh_.point(loop[i]) - mesh_.point(loop[(i+1)%n])).norm();
		length += (len_to_from + len_vn_1_from);

		for (i=0, l=0.0; i<n; ++i) { //fix the boundary vertices		
			angle = l/length * (2.0*M_PI);
			mesh_.property(vuv_, loop[i]) = OpenMesh::Vec2d(0.5*cos(angle)+0.5, 0.5*sin(angle)+0.5);//这个单位圆是以(0.5, 0.5)为圆心, 半径是0.5.
			l += (mesh_.point(loop[i]) - mesh_.point(loop[(i+1)%n])).norm();
		}
		for (i = 0; i < n; ++i) {
			//std::cout << mesh_.property(vuv_, loop[i]) << "\n";
		}

		// 对于from顶点, 特殊 
		double bn_1 = len_to_from / (len_to_from + len_vn_1_from);
		double bt = 1 - bn_1; // = len_vn_1_from / (len_to_from + len_vn_1_from);
		OpenMesh::Vec3d bc_from(bt, 0, bn_1);
		mesh_.property(vuv_, from) = (mesh_.property(vuv_, loop[0]) - mesh_.property(vuv_, loop[n-1])) * bt + mesh_.property(vuv_, loop[n-1]);
		//std::cout << from << ": " << to<< "," << loop[n-2] << "," << loop[n-1] << ": " << bc_from << "\n"; // for test

		// 收集from顶点的One-ring邻域的面上的顶点, 这些顶点是在之前的参数化过程中参数化到落到相应的三角面上的.
		std::vector<TriMesh::VertexHandle> free_vertices;
		OpenMesh::VPropHandleT<bool> is_face_located;
		mesh_.add_property(is_face_located);
		for (i = 0; i <= n-2; ++i) { //循环from顶点的一邻域面.
			TriMesh::FaceHandle fh = mesh_.face_handle(loop_h[i]);
			TriMesh::VertexHandle v0 = from;//mesh_.from_vertex_handle(loop_h[i]);//这个面的3个顶点, 有一个就是from顶点
			TriMesh::VertexHandle v1 = loop[i];//mesh_.to_vertex_handle(loop_h[i]);//另外两个是from顶点的一邻域顶点,
			TriMesh::VertexHandle v2 = loop[i+1];//mesh_.to_vertex_handle(loop_h[i+1]);//这三个顶点的参数坐标已求出了.

			std::vector<TriMesh::VertexHandle> fvset = mesh_.property(fvset_, fh);
			for (std::vector<TriMesh::VertexHandle>::iterator it(fvset.begin()), it_end(fvset.end()); it != it_end; ++it) {
				//初始的时候每一个面上还没有参数化落在上面的顶点, 此时fvset为空.
				if (mesh_.status(*it).deleted() == false) { std::cerr << "Error: \n"; return; } // 这些以前已经参数话过的点应该都是deleted vertices.
				free_vertices.push_back(*it);			

				TriMesh::VertexHandle vf0 = mesh_.property(vf0_, *it); //这个*it顶点在以前的参数化过程中落在的面的三个顶点,
				TriMesh::VertexHandle vf1 = mesh_.property(vf1_, *it); //这三个顶点所构成的面, 刚好应该就是上面循环中的那个面. 
				TriMesh::VertexHandle vf2 = mesh_.property(vf2_, *it);
				if (vf0 != v0 && vf0 != v1 && vf0 != v2) { std::cerr << "Error: parameterziation vf0.\n"; return; }
				if (vf1 != v0 && vf1 != v1 && vf1 != v2) { std::cerr << "Error: parameterziation vf1.\n"; return; }
				if (vf2 != v0 && vf2 != v1 && vf2 != v2) { std::cerr << "Error: parameterziation vf2.\n"; return; }
				mesh_.property(vuv_, *it) = mesh_.property(vuv_, vf0) * (mesh_.property(vbc_, *it))[0]
				+ mesh_.property(vuv_, vf1) * (mesh_.property(vbc_, *it))[1]
				+ mesh_.property(vuv_, vf2) * (mesh_.property(vbc_, *it))[2];
			}			
			(mesh_.property(fvset_, fh)).clear();//这个面上面的顶点先清空, 下面再重新看有哪些顶点参数化落在这里.
		}
		for (unsigned int freev_i = 0, freev_size = free_vertices.size(); freev_i < freev_size; ++ freev_i) {
			mesh_.property(is_face_located, free_vertices[freev_i]) = false;//所有的free vertices初始为还没有定位三角形.
		}


		for (unsigned int j = 1; j <= n-2; ++j) { //循环n-2个面, 这些面是在(from, to)->to之后去除fl和fr后形成的.
			TriMesh::VertexHandle v0 = loop[0];//=to
			TriMesh::VertexHandle vj = loop[j];
			TriMesh::VertexHandle vj1 = loop[j+1];//这三个顶点形成一个三角形, 下面判断free_vertices数组里面的顶点是否落在其中.
			TriMesh::FaceHandle fh = mesh_.face_handle(loop_h[j]);

			for (std::vector<TriMesh::VertexHandle>::iterator freev_it = free_vertices.begin(), freev_end = free_vertices.end(); freev_it != freev_end; ++freev_it) {
				if (mesh_.property(is_face_located, *freev_it) == false) {

					OpenMesh::Vec2d pp = mesh_.property(vuv_, *freev_it);//这个free vertex的参数坐标
					OpenMesh::Vec3d bc = DGP::calc_barycentric_coordinates(mesh_.property(vuv_,v0), mesh_.property(vuv_, vj), mesh_.property(vuv_, vj1), pp);
					//std::cout << bc << std::endl;//for test
					//(0.5, -1.11022e-016, 0.5)//这里是两个特殊的重心坐标,其中的两个很小的小数-1.11022e-016和5.55112e-017都应该判为0的,
					//(0.5, 0.5, 5.55112e-017) //这样的话第一个重心坐标其实是合乎在一个三角形内部(边界上)的内部点的要求的.
					if (fabs(bc[0]) < 1.0e-10) bc[0] = 0;
					if (fabs(bc[1]) < 1.0e-10) bc[1] = 0;
					if (fabs(bc[2]) < 1.0e-10) bc[2] = 0;
					if (!(bc[0]<0) && !(bc[1]<0) && !(bc[2]<0)) { //要求这个bc的三个分量都是>=0

						// For interior vertex, its' barycentric coordinates is bc according to the triangle(v0, vj, vj1).
						//std::cout << *freev_it << ": " << v0 << "," << vj << "," << vj1 << ": " << bc << "\n"; // for test
						// 下面四行就是每一次参数化所应该获得的或是更改的信息.
						mesh_.property(vf_, *freev_it) = fh;
						mesh_.property(vf0_, *freev_it) = v0; mesh_.property(vf1_, *freev_it) = vj; mesh_.property(vf2_, *freev_it) = vj1;
						mesh_.property(vbc_, *freev_it) = bc;

						(mesh_.property(fvset_, fh)).push_back(*freev_it);

						mesh_.property(is_face_located, *freev_it) = true;
					}
				} // 这个free vertex落在这个面上
			} // 一部分free vertex 落在这个面上			
		}

		mesh_.property(vbc_, from) = bc_from;
		mesh_.property(vf0_, from) = to; mesh_.property(vf1_, from) = loop[n-2]; mesh_.property(vf2_, from) = loop[n-1];
		mesh_.property(vf_, from) = mesh_.face_handle(loop_h[n-2]);
		(mesh_.property(fvset_, mesh_.face_handle(loop_h[n-2]))).push_back(from);

		for (unsigned int freev_i = 0, freev_size = free_vertices.size(); freev_i < freev_size; ++ freev_i) { //这里起了判断的作用,估计可以去掉的.
			if (mesh_.property(is_face_located, free_vertices[freev_i]) == false) { //表示这个free vertex这上面的n-2个面中
				std::cout << freev_size << ", " << free_vertices[freev_i];
				std::cerr << "Error: parameterization case 3.\n";// 没有找到投影到其中哪一个上, 这是错误的.
				return;
			}  
		}
		mesh_.remove_property(is_face_located);

	} else { std::cout << "Error: parameterize case4.\n"; 
	} // end of 3 cases
} // endl of function initial_parameterization /**/

void DecimationModel::decimate(unsigned int _n_vertices) //用户期望简化到这个顶点数目
{
	// build priority queue
	TriMesh::VertexIter  v_it  = mesh_.vertices_begin(), v_end = mesh_.vertices_end();
	queue.clear();
	vhierarchy_.clear();
	for (; v_it!=v_end; ++v_it) {//evaluate the cost of every vertex(pair)'s contraction
		// the minimum cost will have the highest priority, which in the front of the queue
		enqueue_vertex(v_it);//或者 enqueue_vertex(v_it.handle());

		// 建立一个森林, 所有的初始节点(也就是叶子节点)就是全部origianl mesh中的原始顶点
		int new_node_handle = vhierarchy_.add_node();
		mesh_.property(vp_leaf_node_handle_, v_it) = new_node_handle;
		mesh_.property(vp_node_handle_, v_it.handle()) = new_node_handle;
		vhierarchy_.node(new_node_handle).set_vertex_handle(v_it.handle());//节点对应哪个顶点
		vhierarchy_.node(new_node_handle).set_active(true);

		assert(new_node_handle == v_it.handle().idx());// 所有初始顶点都有一个和它的idx()一样序号的节点Node,这些Node作为叶子
		//也就是说每一个original mesh上的顶点都对应一个叶子节点, 那叶子节点的序号就是顶点的序号.
	}
	vhierarchy_.set_n_leaves(mesh_.n_vertices());//叶子多少个
	// 初始化每一个半边结构, 使其指向它的opposite vertex handle.//目的是方便halfedge collpase时候快速找到fundamental cut vertex.
	for (TriMesh::HIter he_it(mesh_.halfedges_begin()), he_end(mesh_.halfedges_end()); he_it != he_end; ++he_it) {
		if (mesh_.is_boundary(he_it) == false) {
			mesh_.property(heovh_, he_it) = mesh_.to_vertex_handle(mesh_.next_halfedge_handle(he_it.handle()));
		} else {
			mesh_.property(heovh_, he_it) = TriMesh::VertexHandle(-1);//TriMesh::VertexHandle InvalidVertexHandle;//idx() == -1
		}

		if (mesh_.status(mesh_.edge_handle(he_it)).feature()) { //参数化时候crease vertex的特殊定位
			mesh_.property(hep_heset_, he_it).push_back(he_it.handle());//初始化为包含了自己
		}
	}

	unsigned int nv(mesh_.n_vertices());

	TriMesh::HalfedgeHandle hh;// the target halfedge to be collapsed
	TriMesh::VertexHandle   from, to;// the target halfedge's two vertex
	TriMesh::VVIter         vv_it;

	std::vector<TriMesh::VertexHandle>            one_ring;//from顶点的邻居
	std::vector<TriMesh::VertexHandle>::iterator  or_it, or_end;

	//----一直进行, 达到目标节点数为止. nv可以理解为mesh_中还有这么多个顶点.
	while (nv > _n_vertices&& !queue.empty())
	{	// std::cout << nv << "\t";//for test  
		// Decimate using priority queue:
		//   1) take 1st element of queue
		//   2) collapse this halfedge
		//   3) update queue

		// get 1st element 应该是代价最小的,
		TriMesh::VertexHandle vh = *queue.begin();//enqueue_vertex(v_it);时候是将出半边的起点放进队列中的
		queue.erase(queue.begin());//取出来后就从queue中将其删去

		hh   = mesh_.property(vp_target_, vh);//需要收缩的半边为hh
		to   = mesh_.to_vertex_handle(hh);
		from = mesh_.from_vertex_handle(hh);//from 应该和vh相等, 也就是from = vh;

		// store one-ring vertex neighbor of the from vertex
		one_ring.clear();
		for (vv_it=TriMesh::VVIter(mesh_, from); vv_it; ++vv_it) {
			//for (vv_it(mesh_, from);            vv_it; ++vv_it) //和上面一样的
			//for (vv_it = mesh_.vv_iter(from);   vv_it; ++vv_it) //和上面一样的
			one_ring.push_back(vv_it.handle());
		}

		// perform collapse
		//if (mesh_.is_collapse_ok(hh))
		if (DGP::QEMDecimationT<TriMesh>::is_collapse_legal(mesh_, hh, vp_type_))// is_collapse_legal(hh),已经包括了mesh_.is_collapse_ok(hh)的测试
		{	//之所以还需要这次的检查是因为from顶点的一邻域里面的情况可能改变了,
			//改变是原因可能是一邻域里的点发生了收缩.

			//postprocess_collapse(hh);//这句话必须在mesh_.collapse(hh);之前,因为半边被收缩之后topology就被改变了
			//std::cout << from.idx() << ", " << to.idx() << ": " << mesh_.property(heovh_, hh) << ", " << mesh_.property(heovh_, mesh_.opposite_halfedge_handle(hh)) << ".\n";//for test
			//------------
			// 修正半边收缩之后所指向的opposite vertex handle, 这样做的目的是为了更快地找到fundamental cut vertices
			if (mesh_.is_boundary(hh) == false) { // 图纸上有说明, 这是openmesh中使用的halfedge collapse的特点.
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
			// 找出将要收缩的这半边的左右fundamental cut vertices, 进而找到fundamental cut vertices所对应的Node节点
			TriMesh::VHandle fund_cut_vertex_handle0 = mesh_.property(heovh_, hh);//这些fundamental cut vertex handle很可能是invaild的,其idx() == -1刚好也表示没有Node和其对应.
			TriMesh::VHandle fund_cut_vertex_handle1 = mesh_.property(heovh_, mesh_.opposite_halfedge_handle(hh));
			int fund_cut_node_handle0 = fund_cut_vertex_handle0.idx();//这是严重地基于一个事实: original mesh上的vertexhandle的idx刚好于其对应的节点Node(叶子)的序号一样.
			int fund_cut_node_handle1 = fund_cut_vertex_handle1.idx();//而且fundamental cut vertex就是original mesh上的顶点.这个值也可能是-1.
			// -------------
			if (do_parameterization_) initial_parameterization(hh);

			//----- to collapse the halfedge
			mesh_.collapse(hh);//halfedge collapse, 默认的就是(from, to)->to
			simplified_mesh_.collapse(hh);
			// (1) 更新顶点to的quadric, 按照原QEM的做法:
			mesh_.property(vp_quadric_, to) += mesh_.property(vp_quadric_, from); //(from, to)->to, 但是没有用Vertex Placement Policies// 
			// (2) 更新顶点to的一邻域的边和顶点的类型, 这东西好像效果更差了在hemisphere.off模型中.
			bool update_type = false;
			if (update_type) {
				int n_fe = 0;
				for (TriMesh::VEIter ve_it(mesh_, to); ve_it; ++ve_it) {
					//if (!mesh_.status(ve_it).deleted()) { // 这个判断语句可能是多余了.
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
					n_fe = 0; // 这个顶点四周的feature edges的个数
					for (TriMesh::VEIter ve_it = mesh_.ve_iter(vv_it.handle()); ve_it; ++ve_it) {
						if (mesh_.status(ve_it).feature()) n_fe++;
					}
					if (n_fe == 0) { mesh_.property(vp_type_, vv_it) = DGP::SMOOTH_VFT; } else if (n_fe == 1) { mesh_.property(vp_type_, vv_it) = DGP::DART_VFT; } 
					else if (n_fe == 2) { mesh_.property(vp_type_, vv_it) = DGP::CREASE_VFT; } else if (n_fe >= 3) { mesh_.property(vp_type_, vv_it) = DGP::CORNER_VFT; } 
				}
			}
			// (3) 更新vertex hierarchy
			int new_node_handle = vhierarchy_.add_node();//每次halfedge collapse之后都生成一个新的节点.
			vhierarchy_.node(mesh_.property(vp_node_handle_, from) ).set_parent_node_handle(new_node_handle);
			vhierarchy_.node(mesh_.property(vp_node_handle_, to) ).set_parent_node_handle(new_node_handle);
			vhierarchy_.node(mesh_.property(vp_node_handle_, from) ).set_active(false);
			vhierarchy_.node(mesh_.property(vp_node_handle_, to) ).set_active(false);

			vhierarchy_.node(new_node_handle).set_lchild_node_handle(mesh_.property(vp_node_handle_, from) );
			vhierarchy_.node(new_node_handle).set_rchild_node_handle(mesh_.property(vp_node_handle_, to) );
			vhierarchy_.node(new_node_handle).set_vertex_handle(to);
			vhierarchy_.node(new_node_handle).set_active(true);

			mesh_.property(vp_node_handle_, to) = new_node_handle;//end. to顶点更新其指向的Node.
			assert(vhierarchy_.node(new_node_handle).vertex_handle().idx() == to.idx());

			vhierarchy_.node(new_node_handle).set_fund_cut_node_handle0(fund_cut_node_handle0);
			vhierarchy_.node(new_node_handle).set_fund_cut_node_handle1(fund_cut_node_handle1);

			//std::cout << "fidx2 " << mesh_.status(from).deleted() << ", " << from.is_valid()<< std::endl; //for test
			//上面的输出是 fidex2 1, 1, 这意味这这个from顶点是标志为deleted了, 但是它本身还有效还在存储的数组里面,
			//所以还是有效的, 
			//因为只有其idx_(可以从idx()函数获得)还不等于-1,那就是物理上还有效,甚至还可以使用mesh_.point(from)获得其坐标.
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
				if (do_parameterization_)  //这是一种十分crazy的做法, 每简化一些点就做平滑.
				{	//必须在mesh_.garbage_collection()之前,否则无法使用以前的Handle.
					//local_smooth_parameterization_triangledomain();
					//local_smooth_parameterization_circledomain();					
				}
			}

		} else { // end of if(is_collapse_legal(hh))
			std::cout << "Here?\n";
		}

		// update queue, 一个顶点被收缩之后它的one-ring中的都需要升级它们的代价和哪条半边为最优半边.
		for (or_it=one_ring.begin(), or_end=one_ring.end(); or_it!=or_end; ++or_it)
			enqueue_vertex(*or_it);
	} // end of while(nv > _n_vertices).

	//-------------------检查已经被简化过的mesh_的存储结构是怎么样的
	/*	v_it= mesh_.vertices_begin();
	v_end= mesh_.vertices_end();
	int si(0); 
	for (; v_it != v_end; ++v_it) { //检查所有的点
	if (mesh_.status(v_it).deleted() == false)  {
	//经过检验, 这里只输出被简化而成的base model中的顶点,
	//但是idx()显示其index还是与在original model的一样, 就是初始读入.off文件时候顶点的索引,
	//这要等到mesh_.garbage_collection();被调用之后其索引才会改变,那时之后之前的handle就无效了.
	//这是否说明了被删除的点只是被标志为deleted(通过mesh_.status(v_it.handle()).set_deleted(true)完成), 
	//但其实这被删除的点还是在存储结构里面, 还是有效的.
	std::cout << "检查base model " << v_it.handle().idx() << ": " << mesh_.point(v_it.handle()) << std::endl;
	++si;
	}
	//经过下面的检验可以看出输出了所有从初始读入.off文件时候d读入的顶点,其索引idx()就是那序号
	//那所有顶点的is_valid()都返回1表示在存储结构上还是存在物理有效, 只要idx() != -1就is_valid() == true
	//但是那deleted()对于被删除的顶点返回的是1,对应简化后还在base model中的点(其位置可能被重新计算)是0.
	//std::cout << v_it.handle().idx() << ", " 
	//		  << v_it.handle().is_valid() << ", " << mesh_.status(v_it.handle()).deleted() << ", " 
	//		  << mesh_.point(v_it.handle()) << std::endl;
	}	
	std::cerr << "没有被删除的顶点个数: " << si << "\n";*/
	//------------------------------------------------------------
	//write(std::string("test.pm"));

	// ------------------
	// clean up
	queue.clear();
	// 简化完毕, vertex hierarchy也建立完毕, 其结构就固定了.
	/*// -----------------------
	for (v_it = mesh_.vertices_begin(), v_end = mesh_.vertices_end(); v_it != v_end; ++v_it) {
	if (mesh_.status(v_it).deleted() == false) { // base mesh
	int node_handle = mesh_.property(vp_node_handle_, v_it);
	//std::cout << v_it.handle() << ", " << mesh_.point(v_it) << ", " << node_handle << "; " << "\n";//for test.
	}
	}*/

	//--------------------------------------------------------------
	// 在mesh_.garbage_collection()前重组网格, 参考ModProgMeshT.cc中write函数的做法,用数组来重新映射.
	int N = mesh_.n_vertices();
	vpoints_ = std::vector<TriMesh::Point>(N);//所有顶点的物理坐标, 先存放base mesh中的顶点, 再存放被deleted的顶点
	voldhandles_ = std::vector<TriMesh::VHandle> (N);
	int i = 0;//作为数组vpoint_的下标和计算
	int n_base_vertices = 0, n_detail_vertices = 0;
	OpenMesh::VPropHandleT<int> idx_vps;//mesh_每一个顶点对应数组vpoints_中哪一个元素呢? 就是用idx作为下标标记.
	mesh_.add_property(idx_vps);
	for (v_it = mesh_.vertices_begin(), v_end = mesh_.vertices_end(); v_it != v_end; ++ v_it) {
		if (mesh_.status(v_it).deleted() == false) { 
			vpoints_[i] = mesh_.point(v_it.handle()); 	//多余的.
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

	// 将vertex hierarchy中每一个Node的int point_handle_属性设置为vpoints_数组中的对应下标.
	// 从而将vertex hierarchy中的节点和用于保存全部顶点坐标的vpoints_数组建立起关系.
	for (int k = 0, vhierarchy_size = vhierarchy_.size(); k < vhierarchy_size; ++k) {
		TriMesh::VertexHandle vh = vhierarchy_.node(k).vertex_handle();
		int point_handle = mesh_.property(idx_vps, vh);//最后一次使用这VertexHandle了,在mesh_.garbage_collection();之后就无效了.
		vhierarchy_.node(k).set_point_handle(point_handle);//此后vhierarchy就通过这point_handle来找mesh_上的对应顶点.
		
		vhierarchy_.node(k).set_vertex_handle(TriMesh::VertexHandle(-1));//清空每一个Node的vertex_handle项.
	}
	mesh_.remove_property(idx_vps);

	for (TriMesh::FIter f_it(simplified_mesh_.faces_begin()), f_end(simplified_mesh_.faces_end()); f_it != f_end; ++f_it) {
		if (simplified_mesh_.status(f_it.handle()).deleted() == false) {
			simplified_mesh_.property(fp_fh_sm_, f_it) = f_it.handle(); //没有被删除的面记录其garbage_collection()的对应面.
		} 
	}
	// -----------------------------------------------------------------------
	//std::cout << mesh_.n_vertices() << ", " << mesh_.n_edges() << ", " << mesh_.n_faces() << "\n";//这两个显示点被收缩前的情况,这时候还没有garbage_collection(),
	// now, delete the items marked to be deleted
	//mesh_.garbage_collection();//这个是必须的,否则模型中被set_deleted(true)的点还留在模型中,这是不对的.
	// In the simplification process, the mesh_ was changed, so update the normals of vertices and faces.
	//mesh_.update_normals();//MeshModel::update();
	//std::cout << mesh_.n_vertices() << ", " << mesh_.n_edges() << ", " << mesh_.n_faces() << "\n";
	simplified_mesh_.garbage_collection();//现在mesh_还是保存原来的索引,为了后面还需要用到参数化的信息.
	simplified_mesh_.update_normals();		// NOTE: 因为改为mesh_没有更新, 所以下面的代码要改成用simplified_mesh_来操作了.


	/// ----------
	// 对于简化之后的模型, 根据每一个顶点所对应的Node, 给相应的Node重新设置其指向的vertex handle.
	for (TriMesh::VIter v_it(simplified_mesh_.vertices_begin()), v_end(simplified_mesh_.vertices_end()); v_it != v_end; ++v_it) {
		int node_handle = simplified_mesh_.property(vp_node_handle_sm_, v_it); 
		vhierarchy_.node(node_handle).set_vertex_handle(v_it.handle());
	}

	for (TriMesh::FIter f_it(simplified_mesh_.faces_begin()), f_end(simplified_mesh_.faces_end()); f_it != f_end; ++f_it) {
		mesh_.property(fp_fh_, simplified_mesh_.property(fp_fh_sm_, f_it)) = f_it.handle(); //mesh_的每一个(没有被删除的)面都知道其在simplified_mesh_的对应面.
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
	std::cout << "简化获得的基网格中特征边: " << fe << ".\n";
	std::cout << "简化获得的基网格的顶点中: smooth: " << sv << ", dart " << dv << ", crease: " << crv << ", corner: " << cov << ", boundary: " << bv << ".\n";

	/*
	std::cout << "\n";// for test
	for (int k = 0, vhsize = vhierarchy_.size(); k < vhsize; ++k) {
	DGP::VHierarchyNode n = vhierarchy_.node(k);
	std::cout << n.self_node_handle() << ", " << n.is_active() << ", " << n.point_handle() << ", " << n.vertex_handle() << "; " 
	<< n.parent_node_handle() << ", " << n.lchild_node_handle() << ", " << n.rchild_node_handle() << "; " 
	<< n.fund_cut_node_handle0() << ", " << n.fund_cut_node_handle1() << "\n";
	}*/	
	// 后来改为压平的参数化之后效果不错, 不用smooth了.
	//local_smooth_parameterization_process();
}

void DecimationModel::copy_backup_mesh() { // 复制一份原始网格的拷贝, 尽量不要去修改它.
	mesh2_ = mesh_; TriMesh::ConstFaceIter        f_it(mesh2_.faces_sbegin()), f_end(mesh2_.faces_end());
	TriMesh::ConstFaceVertexIter  cfv_it;

	mesh2_indices_for_render_.clear();
	mesh2_indices_for_render_.reserve(mesh2_.n_faces()*3);

	for (; f_it!=f_end; ++f_it)
	{	// 遍历所有的面, 记录下每一个面的顶点序号
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
	// compute the area associated to a vertex. 计算一个顶点的相关面积, 是一个标量值
	// The simplest such area is the barycentric area, 
	// given by 1/3 the area of all incident triangles. 这里使用的就是这方法.
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

			area += ((Q-P)%(R-P)).norm() * 0.5f * 0.3333f;//这个面的面积的三分之一
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
			//求出半边h_it.handle()所对应的系数wij, 注意这里使用的是mean value coordiante, 而不是harmonic map中的cotangent weight
			TriMesh::VertexHandle v0 = mesh2_.from_vertex_handle(h_it);
			TriMesh::VertexHandle v1 = mesh2_.to_vertex_handle(h_it);
			TriMesh::VertexHandle v2 = mesh2_.to_vertex_handle(mesh2_.next_halfedge_handle(h_it.handle()));// 这以下两句正确的前提是非边界情况
			TriMesh::VertexHandle v3 = mesh2_.to_vertex_handle(mesh2_.next_halfedge_handle(mesh2_.opposite_halfedge_handle(h_it.handle())));
			OpenMesh::Vec3d v1v0 = mesh2_.point(v1) - mesh2_.point(v0);	// 之前这里出现的-1
			OpenMesh::Vec3d v2v0 = mesh2_.point(v2) - mesh2_.point(v0);
			OpenMesh::Vec3d v3v0 = mesh2_.point(v3) - mesh2_.point(v0);
			double v1v0v2_angle = acos(dot(v1v0, v2v0) / (v1v0.norm() * v2v0.norm()));// 矢量v1v0和v2v0的夹角.
			double v3v0v1_angle = acos(dot(v3v0, v1v0) / (v3v0.norm() * v1v0.norm()));// 矢量v3v0和v1v0的夹角.
			double wij = (tan(v1v0v2_angle / 2) + tan(v3v0v1_angle / 2)) / v1v0.norm();
			if (wij < 0) { std::cout << "Error: calc_backupmesh_halfedge_meanvaluecoor(): wij negative.\n"; return; }//mean value coordinate应该是大于零的.
			else mesh2_.property(hmvc_, h_it.handle()) = wij;
		}		
	}
} // end of function calc_backupmesh_halfedge_meanvaluecoor().
void DecimationModel::decimation_process(unsigned int _num, bool _feature_or_not) {
	std::cout << "\n";
	std::cout << "========初始化简化过程==================\n";//这种提示用40个长度的=,前面第9个位置开始文字,即最先有8个=.
	copy_backup_mesh();
	add_required_properties();
	simplified_mesh_ = mesh_;  //simplified_mesh_初始和mesh_一致, 并随着mesh_的简化而简化, 只是simplified_mesh_没有garbage_colection().
	simplified_mesh_.add_property(vp_type_sm_); 
	simplified_mesh_.add_property(vp_node_handle_sm_);
	 
	//DGP::identify_all_sharp_features(mesh_, vp_type_, false, 30);	 
	feature_considered_ = _feature_or_not;
	double const gfxDEGTORAD =  0.01745329251994329577;//角度单位的转换 pi/180 degree to radian
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
		if (mesh_.status(v_it).deleted() == false) { // 后来加入的, 当时考虑到简化网格的原因.

			if (mesh_.is_boundary(v_it)) {
				b++;
			} 
			unsigned int n_fe = 0; // 这个顶点四周的feature edges的个数
			for (TriMesh::VEIter ve_it = mesh_.ve_iter(v_it.handle()); ve_it; ++ve_it) {
				if (mesh_.status(ve_it).feature()) n_fe++;//
				//if (_m.status(ve_it).feature() || _m.is_boundary(ve_it)) n_fe++;//边界边也当成特征边.
				//为了程序的简明, boundary edges and vertices 单独处理.
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
		<< ", CREASE_VFT " << cr << ", CORNER_VFT " << co << "; 其中BOUNDARY " << b << ". \n"; 

	for (TriMesh::VIter v_it = mesh_.vertices_begin(), v_end(mesh_.vertices_end()); v_it != v_end; ++v_it) {
		if (mesh_.status(v_it).deleted() == false) { // 后来加入的
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


	do_parameterization_ = true;//false;// //简化的同时是否进行参数化
	if (do_parameterization_) { mesh_collapsed_ = mesh_; mesh_collapsed_render_or_not_ = false; }
	do_smooth_ = true;//false;//initial parameterization之后是否进行smooth
	n_local_sp_ = 2;//initial parameterization之后进行这么多次的local smooth parameterization 
	std::cout << "--------进入简化过程--------------------\n";
	std::cout << "DecimationModel: decimate to " << _num << ", begin.\n";
	decimate(_num);
	set_mesh_toberendered(&simplified_mesh_, SIMPLIFIED);
	remove_required_properties();	
	std::cout << "DecimationModel: end. " << simplified_mesh_.n_vertices() << " vertices, " << simplified_mesh_.n_edges() << " edges, " << simplified_mesh_.n_faces()    << " faces\n";
	std::cout << "--------离开简化过程--------------------\n";

	if (do_parameterization_) {
		// ------------关于显示参数化的点, 或者是将这些点输出到文件在用PointViewer.exe独立查看.
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

				TriMesh::Point vp = vfp0 * bc[0] + vfp1 * bc[1] + vfp2 * bc[2];//被deleted也就是被parameterized的顶点的3d坐标.

				//file << vp << "\n";//换成是obj格式的话加上<< "v " // 文件输出那些被参数化的点坐标
				//vpoints_for_render_[i] = vp;//被删去的话就会被参数化到基网格的某一个面, 所以这里取其参数坐标

				mesh_collapsed_.set_point(v_it, vp);////
			} else { 
			} 
		}
		//file.close();
	}

	// 基网格上的边, 都初始化为还没有中点采样的.
	mesh_.add_property(empl_);//边的中点是否已经resample了.
	mesh_.add_property(emp_);//如果resample,其3D坐标.
	mesh_.add_property(emp_closest_vh_);
	for (TriMesh::EdgeIter e_it(mesh_.edges_begin()), e_end(mesh_.edges_end()); e_it != e_end; ++ e_it) {
		if (mesh_.status(e_it).deleted() == false) { 
			mesh_.property(empl_, e_it.handle()) = false; // 初始化, 这个边的中点还没有定位.
		}
	}

	std::cout << "========结束简化过程====================\n";
}

// halfedge collapse的逆序, 从base mesh恢复到original mesh,
// 也就是从vertex hierarchy的最上层按简化时候的逆序依次分裂到最下层.
void DecimationModel::sequence_refine()  {
	std::cout << "\n";
	std::cout << "========进入Sequence Refinement=========\n";//这种提示用40个长度的=,前面第9个位置开始文字,即最先有8个=.
	// 选择要split的顶点
	for (int k = vhierarchy_.size() - 1; k >= vhierarchy_.n_leaves(); --k) {// ecol 的逆序分裂

		TriMesh::VertexHandle v1 = vhierarchy_.node(k).vertex_handle(); 

		int v1_node_handle = simplified_mesh_.property(vp_node_handle_sm_, v1);
		assert(k == v1_node_handle);
		//std::cout << v1.idx() << ", " << v1_node_handle << "\n";// for test, 注意, 检查的时候注意那v1_node_handle 不能是-1的, 
		//-1表示这个顶点vertex handle还没有一个和它对应的Node相连.

		// 这个顶点所对应的Node是需要active的
		if (vhierarchy_.node(v1_node_handle).is_active()) { 
			int lc_node_handle = vhierarchy_.node(v1_node_handle).lchild_node_handle();
			int rc_node_handle = vhierarchy_.node(v1_node_handle).rchild_node_handle();
			if (lc_node_handle != -1) { // rc_node_handle != -1 //这个节点不是叶子节点,可以往下split
				//std::cout <<"can be splitted: " << lc_node_handle << ", " << rc_node_handle << "\n";//for test
				//vsplit(to)->(from, to), 也就是vspit(v1)->(v0, v1).
				TriMesh::Point v_from_point = vpoints_[vhierarchy_.node(lc_node_handle).point_handle()];
				TriMesh::VHandle v_from = simplified_mesh_.add_vertex(v_from_point);//这个新增的from点
				TriMesh::VHandle v_to = v1;		//std::cout << vpoints_[vhierarchy_.node(lc_node_handle).point_handle()] << "\n";//for test
				TriMesh::VHandle v_l = TriMesh::VertexHandle(-1);//表示invalid vertex handle
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
	std::cout << "========结束Sequence Refinement=========\n";//这种提示用40个长度的=,前面第9个位置开始文字,即最先有8个=.
}
void DecimationModel::selective_refine() {
	std::cout << "\n";
	std::cout << "========进入vsplit过程==================\n";//这种提示用40个长度的=,前面第9个位置开始文字,即最先有8个=.

	for (TriMesh::VIter v_it = simplified_mesh_.vertices_begin(), v_end = simplified_mesh_.vertices_end(); v_it != v_end; ++v_it) { //每次从森林的有效的顶层往下分裂一层

		TriMesh::VertexHandle v1 = v_it.handle();// v1顶点对应的v1_node_handle这个节点Node做分裂.
		int v1_node_handle = simplified_mesh_.property(vp_node_handle_sm_, v1);
		std::cout << v1.idx() << ", " << v1_node_handle << "\n";// for test, 注意, 检查的时候注意那v1_node_handle 不能是-1的, -1表示这个顶点vertex handle还每一个和它对应的Node相连.
		// 这个顶点所对应的Node是需要active的
		if (vhierarchy_.node(v1_node_handle).is_active()) { 
			int lc_node_handle = vhierarchy_.node(v1_node_handle).lchild_node_handle();
			int rc_node_handle = vhierarchy_.node(v1_node_handle).rchild_node_handle();
			if (lc_node_handle != -1) { // rc_node_handle != -1 //这个节点不是叶子节点,可以往下split
				//std::cout <<"can be splitted: " << lc_node_handle << ", " << rc_node_handle << "\n";//for test
				//vsplit(to)->(from, to), 也就是vspit(v1)->(v0, v1).
				TriMesh::Point v_from_point = vpoints_[vhierarchy_.node(lc_node_handle).point_handle()]; 
				TriMesh::VHandle v_from = simplified_mesh_.add_vertex(v_from_point);//这个新增的from点
				TriMesh::VHandle v_to = v1;											//std::cout << vpoints_[vhierarchy_.node(lc_node_handle).point_handle()] << "\n";//for test
				TriMesh::VHandle v_l = TriMesh::VertexHandle(-1);//表示invalid vertex handle
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
	std::cout << "========结束vsplit过程==================\n";//这种提示用40个长度的=,前面第9个位置开始文字,即最先有8个=.	
}

// ----------------------
// ----------------------
void DecimationModel::local_smooth_parameterization_triangledomain(TriMesh::FaceHandle _fh, int & n_vertices_relocated, int &n_vertices_changeface) //每次以四个等边三角形作为平滑区域
{
	TriMesh::FHandle f_it = _fh; //用f_it这名字上为了方便配合以前的代码
	std::vector<OpenMesh::Vec2d> vuv(18);
	vuv[0] = OpenMesh::Vec2d(1, 0);		vuv[1] = OpenMesh::Vec2d(1.5, 0.8660254);		vuv[2] = OpenMesh::Vec2d(0.5, 0.8660254);
	vuv[3] = OpenMesh::Vec2d(0, 0);		vuv[4] = OpenMesh::Vec2d(2, 0);					vuv[5] = OpenMesh::Vec2d(1, 1.7320508);
	vuv[6] = OpenMesh::Vec2d(-0.5, 0.8660254);		vuv[7] = OpenMesh::Vec2d(0.5, -0.8660254);
	vuv[8] = OpenMesh::Vec2d(1.5, -0.8660254);		vuv[9] = OpenMesh::Vec2d(2.5, 0.8660254);
	vuv[10] = OpenMesh::Vec2d(2, 1.7320508);		vuv[11] = OpenMesh::Vec2d(0, 1.7320508);
	vuv[12] = OpenMesh::Vec2d(-1, 0);				vuv[13] = OpenMesh::Vec2d(-0.5, -0.8660254);
	vuv[14] = OpenMesh::Vec2d(2.5, -0.8660254);		vuv[15] = OpenMesh::Vec2d(3, 0);
	vuv[16] = OpenMesh::Vec2d(1.5, 2.5980762);		vuv[17] = OpenMesh::Vec2d(0.5, 2.5980762);
	// 现在是对面f_it操作.
	//////////////////////////////////////////////////////////////////////////
	// (1) 围绕当前面f_it建立参数域 //建立一个参数面, 由4个等边三角形组成.
	TriMesh::FHIter fh_it(mesh_, f_it);
	TriMesh::HalfedgeHandle h0 = fh_it.handle(); ++fh_it;
	TriMesh::HalfedgeHandle h1 = fh_it.handle(); ++fh_it;
	TriMesh::HalfedgeHandle h2 = fh_it.handle();
	TriMesh::VHandle v0(mesh_.to_vertex_handle(h0)), v1(mesh_.to_vertex_handle(h1)), v2(mesh_.to_vertex_handle(h2));
	mesh_.property(vuv_, v0) = vuv[0];//这几个点是不可能与其它顶点重叠的, 所以可以现在就指定去坐标.
	mesh_.property(vuv_, v1) = vuv[1];
	mesh_.property(vuv_, v2) = vuv[2];

	TriMesh::FaceHandle f0(-1), f1(-1), f2(-1);
	TriMesh::VHandle v3, v4, v5;//(-1)
	if (mesh_.is_boundary(mesh_.opposite_halfedge_handle(h0)) == false) {
		f0 = mesh_.face_handle(mesh_.opposite_halfedge_handle(h0));
		v3 = mesh_.to_vertex_handle(mesh_.next_halfedge_handle(mesh_.opposite_halfedge_handle(h0)));
		//mesh_.property(vuv_, v3) = vuv[3];//这里因为v3/v4/v5可能是指向相同的点,所以先不给出参数坐标,以免重复.
	}
	if (mesh_.is_boundary(mesh_.opposite_halfedge_handle(h1)) == false) {
		f1 = mesh_.face_handle(mesh_.opposite_halfedge_handle(h1));
		v4 = mesh_.to_vertex_handle(mesh_.next_halfedge_handle(mesh_.opposite_halfedge_handle(h1)));				
	}
	if (mesh_.is_boundary(mesh_.opposite_halfedge_handle(h2)) == false) {  //h2的对半边不是边界
		f2 = mesh_.face_handle(mesh_.opposite_halfedge_handle(h2));
		v5 = mesh_.to_vertex_handle(mesh_.next_halfedge_handle(mesh_.opposite_halfedge_handle(h2)));			
	} //这之后的f0,f1,f2仍然可能是.is_valid() == fasle的, 这是边界面的结果.
	TriMesh::FHandle f3, f4, f5, f6, f7, f8;//正三角形参数化区域的extension
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

	// 完成(1)对参数域的建立.
	//////////////////////////////////////////////////////////////////////////
	// (2)
	//参数化到这个面上的所有顶点
	std::vector<TriMesh::VertexHandle> fvset = mesh_.property(fvset_, f_it);//平滑参数化都这个f_it面上的参数化顶点.
	//std::cout << "fvset_'s size " << fvset.size() << ", ";// std::endl; // for test// 
	for (std::vector<TriMesh::VertexHandle>::iterator fvset_it(fvset.begin()), it_end(fvset.end()); fvset_it != it_end; ++fvset_it) {
		//每次 处理一个参数化到此面上的顶点*fvset_it.
		//std::cout << "这次处理的顶点是*fvset_it: " << *fvset_it << ".\n"; //for test
		if (mesh_.property(vf_, *fvset_it) != f_it) { // only to confirm
			std::cout << "Error, para check, faces not match.\n";
		}

		if (mesh2_.is_boundary(*fvset_it)) continue;//如果是边界上的顶点(参数化的时候也参数化到边界上的)就不要进行smooth parameterize.
		if (mesh_.property(vp_type_, *fvset_it) == DGP::CREASE_VFT)  continue;// crease vertex也同样不需要smooth					 

		// 用这个顶点*fvset_it到original mesh中找其一邻域里面的顶点vj, 再判断vj是不是都是在这个三角形的邻域(最多就是4个三角形)里面.
		// 如果是的话, 有效, 就更新其参数化重心坐标.
		bool inside_valid = true; 
		double wij= 0;
		double sum_wij = 0.0;
		OpenMesh::Vec2d sum_wh(0, 0, 0);            
		// 判断*fvset_it的一邻域点vj是否有效
		for (TriMesh::VertexOHalfedgeIter voh_it(mesh2_, *fvset_it); voh_it && inside_valid; ++voh_it) {
			TriMesh::VHandle vj = mesh2_.to_vertex_handle(voh_it.handle());  //由于顶点*fvset_it不是在边上, 所以其出半边都不是在边界的.
			//邻居点vj要么是被deleted了并在那4个面上, 要么是没有被deleted并在那六个顶点上, 否则都算是无效的.
			//std::cout << "vj " << vj << ", ";//for test
			if (mesh_.status(vj).deleted()) { //邻接点vj是被deleted的也就是有被参数化到base complex的某一个面上的.
				TriMesh::FaceHandle pf = mesh_.property(vf_, vj);//邻点vj在initial parameterize是所参数化到的面
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

				if (inside_valid == true) { // == false时候就没有必要求这个参数坐标了.
					mesh_.property(vuv_, vj)
						= mesh_.property(vuv_, mesh_.property(vf0_, vj)) * (mesh_.property(vbc_, vj))[0]
					+ mesh_.property(vuv_, mesh_.property(vf1_, vj)) * (mesh_.property(vbc_, vj))[1]
					+ mesh_.property(vuv_, mesh_.property(vf2_, vj)) * (mesh_.property(vbc_, vj))[2];
				} else { //在之前的检查中不在这9个有效面内, 这里还要进一步检查
					if (mesh_.property(vp_type_, vj) == DGP::CREASE_VFT) {

						if (h3.is_valid()) {
							if (mesh_.property(vp_feature_edge_, vj) == mesh_.edge_handle(mesh_.next_halfedge_handle(h3))) {
								mesh_.property(vuv_, v3) = vuv[3]; mesh_.property(vuv_, v6) = vuv[6]; 
								inside_valid = true;
							}
							if (mesh_.property(vp_feature_edge_, vj) == mesh_.edge_handle(mesh_.prev_halfedge_handle(h3))) { //v2的vuv_不用设
								//这种情况应该包含在f3 or f9里面了
							}
						}
						if (h4.is_valid()) {
							if (mesh_.property(vp_feature_edge_, vj) == mesh_.edge_handle(mesh_.next_halfedge_handle(h4))
								|| mesh_.property(vp_feature_edge_, vj) == mesh_.edge_handle(mesh_.prev_halfedge_handle(h4))) { //v0的vuv_不用设
									mesh_.property(vuv_, v3) = vuv[3]; mesh_.property(vuv_, v7) = vuv[7];
									inside_valid = true;//if ((*fvset_it).idx() == 6424) { std::cout << "tri smooth crease vj 4.\n"; }
							}
						}
						if (h5.is_valid()) {
							if (mesh_.property(vp_feature_edge_, vj) == mesh_.edge_handle(mesh_.next_halfedge_handle(h5))
								|| mesh_.property(vp_feature_edge_, vj) == mesh_.edge_handle(mesh_.prev_halfedge_handle(h5))) { //v0的vuv_不用设
									mesh_.property(vuv_, v4) = vuv[4]; mesh_.property(vuv_, v8) = vuv[8];
									inside_valid = true;//if ((*fvset_it).idx() == 6424) { std::cout << "tri smooth crease vj 5.\n"; }
							}
						}
						if (h6.is_valid()) {
							if (mesh_.property(vp_feature_edge_, vj) == mesh_.edge_handle(mesh_.next_halfedge_handle(h6))
								|| mesh_.property(vp_feature_edge_, vj) == mesh_.edge_handle(mesh_.prev_halfedge_handle(h6))) { //v1的vuv_不用设
									mesh_.property(vuv_, v4) = vuv[4]; mesh_.property(vuv_, v9) = vuv[9];
									inside_valid = true;//if ((*fvset_it).idx() == 6424) { std::cout << "tri smooth crease vj 6.\n"; }
							}
						}
						if (h7.is_valid()) {
							if (mesh_.property(vp_feature_edge_, vj) == mesh_.edge_handle(mesh_.next_halfedge_handle(h7))
								|| mesh_.property(vp_feature_edge_, vj) == mesh_.edge_handle(mesh_.prev_halfedge_handle(h7))) { //v1的vuv_不用设
									mesh_.property(vuv_, v5) = vuv[5]; mesh_.property(vuv_, v10) = vuv[10];
									inside_valid = true;//if ((*fvset_it).idx() == 6424) { std::cout << "tri smooth crease vj 7.\n"; }
							}
						}
						if (h8.is_valid()) {
							if (mesh_.property(vp_feature_edge_, vj) == mesh_.edge_handle(mesh_.next_halfedge_handle(h8))
								|| mesh_.property(vp_feature_edge_, vj) == mesh_.edge_handle(mesh_.prev_halfedge_handle(h8))) { //v2的vuv_不用设
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
										   + vuv[3]  // vbc_[1] == 0 // 这里犯了个典型错误啊, crease点参数化到crease edge上只有两个邮箱vbc_, 但已经不已经就是vbc_[1] == 0!
										   + mesh_.property(vuv_, mesh_.property(vf2_, vj)) * (mesh_.property(vbc_, vj))[2];*/
							mesh_.property(vuv_, vj) = mesh_.property(vuv_, mesh_.property(vf0_, vj)) * (mesh_.property(vbc_, vj))[0]
							+ mesh_.property(vuv_, mesh_.property(vf1_, vj)) * (mesh_.property(vbc_, vj))[1]
							+ mesh_.property(vuv_, mesh_.property(vf2_, vj)) * (mesh_.property(vbc_, vj))[2];
							// 上面那里有一项vbc_[i]等于零, 刚好其对应的顶点参数化坐标也不用设置.
						}

					} //else inside_valid = false; 本来就已经==false  // end of if-else vj is crease 
				} // end of if-else. 被deleted的话是否有效
			} else { //邻接点vj没有被deleted的也就是它还是base complex上的某一个顶点.
				if ((vj == v0) || (vj == v1) || (vj == v2)) { //mesh_.property(vuv_, vj);当vj就等于那六个顶点的话那么其参数坐标就已经在上面设置的了.
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

			if (inside_valid) { //vi的其中一个vj有效//
				//std::cout << "vj is valid, " << mesh_.property(vuv_, vj) << ";\t";//for test
				//求出半边voh_it.handle()所对应的系数wij, 注意这里使用的是mean value coordiante, 而不是harmonic map中的cotangent weight
				sum_wij += mesh2_.property(hmvc_, voh_it);
				sum_wh += mesh2_.property(hmvc_, voh_it) * mesh_.property(vuv_, vj);							
			} else { // 这情况应该不多才对, 毕竟展开的面已经不少了.
				std::cout << "三角形平滑In face " << f_it << "展开的面有: " << f0 << " " << f1 << " " << f2 << " " // for test
					<< f3 << " " << f4 << " " << f5 << " " << f6 << " " << f7 << " " << f8 << " "
					<< f9 << " " << f10 << " " << f11 << " " << f12 << " " << f13 << " " << f14  << ".\n"; 
				std::cout << "(vi, vj)=(" << *fvset_it << ", " << vj << "), vj's vf_" << mesh_.property(vf_, vj) << ", vj's type" << mesh_.property(vp_type_, vj) << ".\n";  
				TriMesh::FaceFaceIter fit = mesh_.ff_iter(mesh_.property(vf_, vj)); std::cout << fit.handle() << " ";
				++ fit; std::cout << fit.handle() << " ";++ fit; std::cout << fit.handle() << ":vj's one-ring faces.\n";
				/**/
			}
		} // end of for.循环完毕了顶点*fvset_it的一邻域顶点vj, 再判断顶点*fvset_it是否有效,也就是说是不是都在一邻域里面.

		if (inside_valid == true) { 
			// 如若vi的所有vj都有效的话就求出vi(也就是*fvset_it)顶点的新参数化坐标, 并找出其落在哪个三角形中.
			OpenMesh::Vec2d old_uv = mesh_.property(vuv_, *fvset_it); 
			mesh_.property(vuv_, *fvset_it) = sum_wh * (1.0/sum_wij);
			// std::cout << "vi " << mesh_.property(vuv_, *fvset_it) << std::endl; // for test
			// 重新判断顶点*fvset_it是参数到这四个面(后来扩展到10个面)中哪一个三角形上.
			std::vector<TriMesh::FHandle> fvector(22);
			fvector[0] = f_it; fvector[1] = f0; fvector[2] = f1; fvector[3] = f2;
			fvector[4] = f3; fvector[5] = f4; fvector[6] = f5; fvector[7] = f6;fvector[8] = f7; fvector[9] = f8; 
			fvector[10] = f9; fvector[11] = f10; fvector[12] = f11; fvector[13] = f12;fvector[14] = f13; fvector[15] = f14; 
			fvector[16] = f15; fvector[17] = f16; fvector[18] = f17; fvector[19] = f18;fvector[20] = f19; fvector[21] = f20; 
			bool is_face_located = false;
			for (size_t i = 0; i < 16 && is_face_located == false; i++) { //循环四个可能的面 //22
				if ((fvector[i]).is_valid() == false ) continue;//f0, f1, f2这三个面都可能不一定存在的.
				TriMesh::VHandle vx, vy, vz;
				if (i == 0) {// fvector[0]所对应的三个顶点
					vx = v0; vy = v1; vz = v2; 
				} else if (i == 1) { // fvector[1]所对应的三个顶点
					vx = v0; vy = v2; vz = v3; mesh_.property(vuv_, v3) = vuv[3];
				} else if (i == 2) { // fvector[2]所对应的三个顶点
					vx = v0; vy = v4; vz = v1; mesh_.property(vuv_, v4) = vuv[4];
				} else if (i == 3) { // fvector[3]所对应的三个顶点
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
				if (!(bc[0]<0) && !(bc[1]<0) && !(bc[2]<0)) { //要求这个bc的三个分量都是>=0
					//std::cout << "BC: " << bc << std::endl; // for test
					// 下面四行就是每一次参数化所应该获得的或是更改的信息.
					mesh_.property(vf0_, *fvset_it) = vx; mesh_.property(vf1_, *fvset_it) = vy; mesh_.property(vf2_, *fvset_it) = vz;
					mesh_.property(vbc_, *fvset_it) = bc;
					if (f_it != fvector[i]) { //mesh_.property(vf_, *fvset_it) 就等于f_it
						n_vertices_changeface++;//前后面是否一样

						// 删去其在原来那个面(f_it)上的记录
						std::vector<TriMesh::VHandle>::iterator to_be_earsed 
							= remove((mesh_.property(fvset_, f_it)).begin(), (mesh_.property(fvset_, f_it)).end(), TriMesh::VHandle((*fvset_it).idx()));
						(mesh_.property(fvset_, f_it)).erase(to_be_earsed, (mesh_.property(fvset_, f_it)).end());

						(mesh_.property(fvset_, fvector[i])).push_back(*fvset_it);
						mesh_.property(vf_, *fvset_it) = fvector[i];
					}
					is_face_located = true; 
					n_vertices_relocated++;//这次smooth para中重新定位的顶点数.
				}
			} // 循环完10个可能的面之后还没有合适的话就出错了.//但是也没有关系, 因为我上面是当参数化到别的面上时候才从原来面上删去记录, 
			if (is_face_located == false) { 
				std::cerr << "Error:tri smooth para face relocated," << old_uv << ", " << mesh_.property(vuv_, *fvset_it) << "; " << *fvset_it << ".\n"; 
				//mesh_.property(vuv_, *fvset_it) = old_uv;
				//test_vertex_ = *fvset_it;
			}// 现在也没有删去其在原来面上的记录.
		} // end of if (inside_valid == true). 处理完参数化到这个面f_it上的一个点的relocate.

	} // end of for. 处理完参数化到这个面f_it上的每一个点*fvset_it的平滑和relocate.
}

void DecimationModel::local_smooth_parameterization_triangledomain() {
	std::cout << "进入函数 local_smooth_parameterization_triangledomain(), before smooth parameter.\n";
	int n_vertices_base_mesh = 0;
	int n_vertices_deleted = 0, n_vertices_paraed = 0;
	for (TriMesh::VertexIter v_it(mesh_.vertices_begin()), v_end(mesh_.vertices_end()); v_it != v_end; ++v_it) {
		if (mesh_.status(v_it.handle()).deleted()) { 
			n_vertices_deleted++;
		} else n_vertices_base_mesh++;
		// std::cout << mesh_.property(vbc_, v_it.handle()) << "\t"; //能成功输出, 表示这个custom properties可用.
	}
	for (TriMesh::FaceIter f_it(mesh_.faces_begin()), f_end(mesh_.faces_end()); f_it != f_end; ++f_it) {
		if (mesh_.status(f_it.handle()).deleted() == false) {
			n_vertices_paraed += (mesh_.property(fvset_, f_it.handle())).size();//参数化到这个面上的顶点数.
		}
	}// 正常的情况是被deleted的顶点数和被参数化的顶点数目一致.
	std::cout << "para check: V " << mesh_.n_vertices() << "= " << n_vertices_base_mesh << " + " << n_vertices_deleted << "(=" << n_vertices_paraed << ").\n";

	// -------------
	int n_edges_nondeleted(0), n_edges_deleted(0), count_crosspatches(0), n_edges_original(0);
	for (TriMesh::EdgeIter e_it(mesh_.edges_begin()), e_end(mesh_.edges_end()); e_it != e_end; ++ e_it) {
		if (mesh_.status(e_it).deleted() == true) {
			n_edges_deleted++;
			TriMesh::VHandle v0 = mesh2_.to_vertex_handle(mesh2_.halfedge_handle(e_it, 0));
			TriMesh::VHandle v1 = mesh2_.to_vertex_handle(mesh2_.halfedge_handle(e_it, 1));
			if (mesh_.property(vf_, v0) != mesh_.property(vf_, v1)) {
				count_crosspatches++;//多少个边不是参数化到同一个面上的.
			}
		} else { //这个边没有被删去,但是可能那两个端点都改变了,也就是说虽存在都亦非之前那个原始边了, 是在简化过程中生成的.
			n_edges_nondeleted++;
			TriMesh::VHandle v0 = mesh_.to_vertex_handle(mesh_.halfedge_handle(e_it, 0));
			TriMesh::VHandle v1 = mesh_.to_vertex_handle(mesh_.halfedge_handle(e_it, 1));
			if (((v0 == mesh2_.to_vertex_handle(mesh2_.halfedge_handle(e_it, 0))) && (v1 == mesh2_.to_vertex_handle(mesh2_.halfedge_handle(e_it, 1))))
				|| ((v1 == mesh2_.to_vertex_handle(mesh2_.halfedge_handle(e_it, 0))) && (v0 == mesh2_.to_vertex_handle(mesh2_.halfedge_handle(e_it, 1))))) {
					// 这个边在简化的过程当中没有发生变化,其两个端点还是简化之间的那两个.
					n_edges_original ++;//有多少个边没有变化
			}
		}
	}
	std::cout << "para check: E " << mesh_.n_edges() << "= " << n_edges_nondeleted << " + " << n_edges_deleted 
		<< "; crop " << count_crosspatches << ", ori " << n_edges_original << ".\n";

	// ---------//下面进行平滑参数化 ---------------------
	std::cout << "下面开始一共" << n_local_sp_ << "次的局部平滑.\n";
	for ( int i_iterance = 0; i_iterance < n_local_sp_; ++i_iterance) { // begin to locally smooth parameterize, n_local_sp_ times
		std::cout << "Begin the " << i_iterance << " time to locally smooth parameterization.\n";
		int n_vertices_relocated = 0; // 记录在这次smooth parameterization中有多少个顶点被重新定位了.
		int n_vertices_changeface = 0;// 有多少个顶点在局部平滑时候移到另外的平面了
		for (TriMesh::FaceIter f_it(mesh_.faces_begin()), f_end(mesh_.faces_end()); f_it != f_end; ++f_it) {
			if (mesh_.status(f_it.handle()).deleted() == false) { //循环每一个没有被deleted的面,也就是base complex上的面
				local_smooth_parameterization_triangledomain(f_it.handle(), n_vertices_relocated, n_vertices_changeface);
			} // end of if.判断f_it.handle()是否被deleted, 也就是否是base complex上的一个面			
		} // end of for.对所有面循环一次, 循环时候对面上的所有参数化的点进行一个smooth
		std::cout << "End the " << i_iterance << " time to locally smooth parameterization. " 
			<< (n_vertices_relocated * 100.0/ n_vertices_paraed) << "% vertices relocated, " << n_vertices_changeface*100/n_vertices_paraed << "% change face.\n";

	} // end of n_local_sp_ times local smooth parameterization.
	n_vertices_paraed = 0;
	for (TriMesh::FaceIter f_it(mesh_.faces_begin()), f_end(mesh_.faces_end()); f_it != f_end; ++f_it) {
		if (mesh_.status(f_it.handle()).deleted() == false) {
			n_vertices_paraed += (mesh_.property(fvset_, f_it.handle())).size();//参数化到这个面上的顶点数.
		}
	}
	std::cout << "para check: n_vertices_parameterized " << n_vertices_paraed << "\n";
	std::cout << "离开函数 local_smooth_parameterization_triangledomain().\n";
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

		unsigned int n = loop_h.size();//有多少个边界顶点
		unsigned int i = 0;  
		double length = 0.0, rou_angle = 0.0;  
		std::vector<double> vec_angle;
		for (i=0 ; i<n; ++i) {//求出总的长度, 周长 
			length += (OpenMesh::vector_cast<OpenMesh::Vec3d>(mesh_.point(loop[i])) - OpenMesh::vector_cast<OpenMesh::Vec3d>(mesh_.point(loop[(i+1)%n]))).norm();

			OpenMesh::Vec3d d1 = (OpenMesh::vector_cast<OpenMesh::Vec3d>(mesh_.point(loop[i])) - OpenMesh::vector_cast<OpenMesh::Vec3d>(mesh_.point(vh))).normalize();
			OpenMesh::Vec3d d2 = (OpenMesh::vector_cast<OpenMesh::Vec3d>(mesh_.point(loop[(i+1)%n])) - OpenMesh::vector_cast<OpenMesh::Vec3d>(mesh_.point(vh))).normalize();
			rou_angle += acos(dot(d1, d2));
			vec_angle.push_back(acos(dot(d1, d2)));
		}  //std::cout << "angle " << rou_angle << ".\n";//std::cout << "len: " << length << std::endl; //for test	
		if (vec_angle.size() != n) { std::cout << "Error: 容量不一致.\n"; }

		/* // 尝试一, 用单位圆做参数域
		double l = 0.0, angle = 0.0;
		double sum_wij = 0.0;
		OpenMesh::Vec2d sum_wh(0, 0, 0);
		for (i = 0; i < n; ++i) { //顶点vi在one-ring上的有n个邻接点// fix the boundary/one-ring vertices, 				
		angle = l / length * (2.0*M_PI); //std::cout << "ang: " << angle << std::endl;
		mesh_.property(vuv_, loop[i]) = OpenMesh::Vec2d(0.5*cos(angle)+0.5, 0.5*sin(angle)+0.5);
		l += (mesh_.point(loop[i]) - mesh_.point(loop[(i+1)%n])).norm();

		// 求出vi和vj这半边上的系数, 注意这里使用的是mean value coordiante, 而不是harmonic map中的cotangent weight
		TriMesh::VertexHandle v0 = vh; // = mesh_.from_vertex_handle(loop_h[i]);
		TriMesh::VertexHandle v1 = loop[i];// =mesh_.to_vertex_handle(loop_h[i]);
		TriMesh::VertexHandle v2 = loop[(i+1)%n];//这以下两句正确的前提是非边界情况
		TriMesh::VertexHandle v3 = loop[(i+n-1)%n];
		OpenMesh::Vec3d v1v0 = mesh_.point(v1) - mesh_.point(v0);	// 之前这里出现的-1
		OpenMesh::Vec3d v2v0 = mesh_.point(v2) - mesh_.point(v0);
		OpenMesh::Vec3d v3v0 = mesh_.point(v3) - mesh_.point(v0);
		double v1v0v2_angle = acos(dot(v1v0, v2v0) / (v1v0.norm() * v2v0.norm()));//矢量v1v0和v2v0的夹角.
		double v3v0v1_angle = acos(dot(v3v0, v1v0) / (v3v0.norm() * v1v0.norm()));//矢量v3v0和v1v0的夹角.
		double wij = (tan(v1v0v2_angle / 2) + tan(v3v0v1_angle / 2)) / v1v0.norm();
		if (wij < 0) { std::cout << "Error: wij is negative.\n"; return; }

		sum_wij += wij;
		sum_wh += wij * mesh_.property(vuv_, loop[i]);
		}
		mesh_.property(vuv_, vh) = sum_wh * (1.0/sum_wij);//这就是vi顶点在单位圆上的参数化坐标.
		*///std::cout << mesh_.property(vuv_, vh) << ", " << n << ":" ;// for test.
		// 尝试二, 压平, 角度缩放, 但是好像改进的效果不明显啊, 代码看上去简单一点
		double angle_scale_ratio = 2 * M_PI / rou_angle; //缩放比率
		double temp_sum_angle = 0.0, l = 0;
		for (i = 0; i < n; ++i) {
			temp_sum_angle += (vec_angle[i] * angle_scale_ratio);
			l = (OpenMesh::vector_cast<OpenMesh::Vec3d>(mesh_.point(loop[(i+1)%n])) - OpenMesh::vector_cast<OpenMesh::Vec3d>(mesh_.point(vh))).norm(); 
			l *= angle_scale_ratio;
			mesh_.property(vuv_, loop[(i+1)%n]) = OpenMesh::Vec2d(l*cos(temp_sum_angle), l*sin(temp_sum_angle));
		}
		mesh_.property(vuv_, vh) = OpenMesh::Vec2d(0, 0);

		for (i = 0; i < n; ++i)  // 对顶点vh的一邻域里面的每一个面loop_f[i].
		{	//std::cout << "In face " << loop_f[i] << ", ";
			std::vector<TriMesh::VHandle> fvset = mesh_.property(fvset_, loop_f[i]);
			for (std::vector<TriMesh::VHandle>::iterator fvset_it = fvset.begin(), fvset_end = fvset.end(); 
				fvset_it != fvset_end; ++fvset_it) //  对参数化到面loop_f[i]上的每一个顶点.
			{ // TriMesh::VertexHandle *fvset_it是操作的对象.
				if (mesh2_.is_boundary(*fvset_it) == true) continue; //这两种情况的特殊点不参与平滑.
				if (mesh_.property(vp_type_, *fvset_it) == DGP::CREASE_VFT) continue; 

				bool one_ring_vertex_is_valid = true; //先假设顶点*fvset_it的一邻域顶点都在这个参数圆内, 这就是有效.
				double sum_wij = 0.0; OpenMesh::Vec2d sum_wh(0, 0, 0);
				for (TriMesh::VOHIter voh_it(mesh2_, *fvset_it); voh_it; ++voh_it) 
				{ 
					TriMesh::VHandle vj = mesh2_.to_vertex_handle(voh_it.handle());
					std::vector<TriMesh::FHandle>::iterator result = find(loop_f.begin(), loop_f.end(), mesh_.property(vf_, vj)); 
					if (result != loop_f.end()) { //这个顶点vj是否在这个参数圆里?
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
				} // end of for 顶点*fvset_it的一邻域顶点.
				//std::cout << one_ring_vertex_is_valid << "; ";
				if (one_ring_vertex_is_valid == true) {
					mesh_.property(vuv_, *fvset_it) = sum_wh * (1.0 / sum_wij); //这个点被平滑之后的新参数坐标.

					// 给顶点*fvset_it找出它现在落在哪一个面上.
					bool relocated = false;
					for (int ii = 0; ii < n && relocated == false; ++ii) { // loop_f[ii]
						TriMesh::VHandle vx = vh, vy = loop[ii], vz = mesh_.to_vertex_handle(mesh_.next_halfedge_handle(loop_h[ii])); 
						OpenMesh::Vec3d bc 
							= DGP::calc_barycentric_coordinates(mesh_.property(vuv_, vx), mesh_.property(vuv_, vy), mesh_.property(vuv_, vz), mesh_.property(vuv_, *fvset_it));

						if (fabs(bc[0]) < 1.0e-10) bc[0] = 0; if (fabs(bc[1]) < 1.0e-10) bc[1] = 0; if (fabs(bc[2]) < 1.0e-10) bc[2] = 0;
						if (!(bc[0]<0) && !(bc[1]<0) && !(bc[2]<0)) { //要求这个bc的三个分量都是>=0
							//std::cout << "BC: " << bc << std::endl; // for test
							// 下面四行就是每一次参数化所应该获得的或是更改的信息.
							mesh_.property(vf0_, *fvset_it) = vx; mesh_.property(vf1_, *fvset_it) = vy; mesh_.property(vf2_, *fvset_it) = vz;
							mesh_.property(vbc_, *fvset_it) = bc;//std::cout << "a ";
							if (loop_f[i] != loop_f[ii]) { //mesh_.property(vf_, *fvset_it) 就等于loop_f[i].
								// 删去其在原来那个面(loop_f[i])上的记录
								std::vector<TriMesh::VHandle>::iterator to_be_earsed 
									= remove((mesh_.property(fvset_, loop_f[i])).begin(), (mesh_.property(fvset_, loop_f[i])).end(), TriMesh::VHandle((*fvset_it).idx()));
								(mesh_.property(fvset_, loop_f[i])).erase(to_be_earsed, (mesh_.property(fvset_, loop_f[i])).end());
								//std::cout << "b ";
								(mesh_.property(fvset_, loop_f[ii])).push_back(*fvset_it); //std::cout << "bb ";
								mesh_.property(vf_, *fvset_it) = loop_f[ii];//std::cout << "c ";
							}
							relocated = true; 
							//n_vertices_relocated++;//这次smooth para中重新定位的顶点数.
						} //std::cout << "d "; 
					} //std::cout << " ok. ";
					if (relocated == false ) { std::cout << "Error: relocated, it seems don't matter.\n"; }
				} // end of if. 这个顶点*fvset_it的一邻域在有效的圆里面.

			} // end of for.对参数化到面loop_f[i]上的每一个顶点*fvset_it.
		} // end of for.对顶点vh的一邻域里面的每一个面loop_f[i].
	} // end of if 顶点vh不是在边界上 
}
void DecimationModel::local_smooth_parameterization_circledomain() {
	std::cout << "进入函数 local_smooth_parameterization_circledomain.\n";
	for (TriMesh::VIter v_it(mesh_.vertices_begin()), v_end(mesh_.vertices_end()); v_it != v_end; ++v_it) 
	{ //以v_it为圆心建参数域
		local_smooth_parameterization_circledomain(v_it.handle());		
	} // end of for.对基网格上的每一个顶点.
	std::cout << "离开函数 local_smooth_parameterization_circledomain.\n";
} // end of function local_smooth_parameterization_circledomain().
void DecimationModel::local_smooth_parameterization_process() {
	if (do_parameterization_)  {	
		if (do_smooth_) {//////////平滑和中点采样都必须在mesh_.garbage_collection()之前,否则无法使用以前的Handle.
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

					TriMesh::Point vp = vfp0 * bc[0] + vfp1 * bc[1] + vfp2 * bc[2];//被deleted也就是被parameterized的顶点的3d坐标.

					mesh_collapsed_.set_point(v_it, vp);////
				} else {
					mesh_collapsed_.set_point(v_it, mesh_.point(v_it.handle()));//没有被删去的顶点就是基网格的顶点,可以直接取其坐标就可以了.
				}
			} // end of for
		} // end of if (do_smooth_)
	} // end of if ()
}