
// ***************************************************************
//  
// VHierarchy 
// Copyright (C) 2007 - by rencanjiang. All Rights Reserved
// rencanjiang@163.com	
// -------------------------------------------------------------
//  
// version: 1.0
// data: 12/10/2007
// ***************************************************************
//  
// File description:
// ʵ��һ��vertex hierarchy, ���ڼ�¼�򻯵Ĺ��̺�PM�ṹ�Ļָ�.
// hierarchy��һ�����νṹ,�ڵ���VHierarchyNode.
// ***************************************************************

#ifndef dgpstudio_vhierarchy_h
#define dgpstudio_vhierarchy_h

#include <vector>
#include "OpenMeshMeshType.h"

namespace DGP {

class VHierarchyNode {
private: //ע��, ���ﶼ��ʹ��handle�����������Ϊindex����˼��.
	// �����ĸ���basic properties
	int self_node_handle_;//������vertex hierarchy��nodes_�����µ��±�.
	bool active_;

	// ��VertexHandle��Point��ӳ��mappings.
	TriMesh::VertexHandle vh_;

	int point_handle_;//�ڼ�֮��Ὠ��һ��std::vector<TriMesh::Point> vpoints_�������ڴ��ԭģ���е����ж���,�ȴ�base mesh���ٴ�deleted��,��������±�. 
	
	// ������Node��vertex hierarchy�е����ӹ�ϵ. Note:����ָ��node���
	int parent_node_handle_; //ָ���׽ڵ�Ľڵ��Ҳ��������vertex hierarchy�е�������±�.
	int lchild_node_handle_, rchild_node_handle_;
	int fund_cut_node_handle_[2];//(from, to)-> to, (v0, v1)->v1ʱ���vl^, vr^������ԭʼ��(��original mesh�ϵ�)����Ӧ��Node���.

public:
	VHierarchyNode() {
		self_node_handle_ = -1;
		active_ = false; 

		vh_ = TriMesh::VertexHandle(-1);//invalid
		point_handle_ = -1;

		parent_node_handle_ = -1; 
		lchild_node_handle_ = -1; rchild_node_handle_ = -1; 
		fund_cut_node_handle_[0] = -1; fund_cut_node_handle_[1] = -1;
	}
	
	int self_node_handle() { return self_node_handle_; }
	void set_self_node_handle(int _nh) { self_node_handle_ = _nh; }

	bool is_active() const { return active_; }
	void set_active(bool _a) { active_ = _a; }

	TriMesh::VertexHandle vertex_handle() const { return vh_; }
	void set_vertex_handle(TriMesh::VertexHandle _vh) { vh_ = _vh; }

	int point_handle() { return point_handle_; }
	void set_point_handle(int _ph) { point_handle_ = _ph; }

	int parent_node_handle() const { return parent_node_handle_; }
	void set_parent_node_handle(int _ph) { parent_node_handle_ = _ph; }

	int lchild_node_handle() const { return lchild_node_handle_; }
	void set_lchild_node_handle(int _lch) { lchild_node_handle_ = _lch; }
	int rchild_node_handle() const { return rchild_node_handle_; }
	void set_rchild_node_handle(int _rch) { rchild_node_handle_ = _rch; }

	int fund_cut_node_handle0() { return fund_cut_node_handle_[0]; }
	void set_fund_cut_node_handle0(int _vh) { fund_cut_node_handle_[0] = _vh; }
	int fund_cut_node_handle1() { return fund_cut_node_handle_[1]; }
	void set_fund_cut_node_handle1(int _vh) { fund_cut_node_handle_[1] = _vh; }
};

class VHierarchy {
private:
	std::vector<VHierarchyNode> nodes_; //���������±���Ϊnode_handle
	int n_leaves_;//the leaves' num of this vertex hierarchy

public:
	void clear() { nodes_.clear(); } //�������ɭ��, ÿ��Ҫ���½���ɭ��ǰ�������.

	int add_node() { //nodes_������±����0��nodes_.size()-1.
		nodes_.push_back(VHierarchyNode());
		nodes_[nodes_.size() - 1].set_self_node_handle(nodes_.size() - 1); 
		return nodes_.size() - 1;
	}
	VHierarchyNode& node(int _node_handle) {
		return nodes_[_node_handle];
	}
	int size() {		// ɭ�ִ�С
		return nodes_.size();
	}
	int n_leaves() {	// Ҷ�Ӷ���
		return n_leaves_;
	}
	void set_n_leaves(int _n_leaves) { //�������ֻ�ǵ���һ�ε�.
		n_leaves_ = _n_leaves;//���ǳ�ʼΪoriginal mesh�ϵ����ж��㶼����һ��Node��ʱ��.
	}
};
} // end of namespace DGP
#endif //dgpstudio_vhierarchy_h
