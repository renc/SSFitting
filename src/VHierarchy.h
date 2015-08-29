
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
// 实现一个vertex hierarchy, 用于记录简化的过程和PM结构的恢复.
// hierarchy是一个树形结构,节点是VHierarchyNode.
// ***************************************************************

#ifndef dgpstudio_vhierarchy_h
#define dgpstudio_vhierarchy_h

#include <vector>
#include "OpenMeshMeshType.h"

namespace DGP {

class VHierarchyNode {
private: //注意, 这里都是使用handle这个概念来作为index的意思的.
	// 下面四个是basic properties
	int self_node_handle_;//自身在vertex hierarchy的nodes_数组下的下标.
	bool active_;

	// 和VertexHandle和Point的映射mappings.
	TriMesh::VertexHandle vh_;

	int point_handle_;//在简化之后会建立一个std::vector<TriMesh::Point> vpoints_数组用于存放原模型中的所有顶点,先存base mesh的再存deleted的,这就是那下标. 
	
	// 下面是Node在vertex hierarchy中的连接关系. Note:都是指向node序号
	int parent_node_handle_; //指向父亲节点的节点号也就是其在vertex hierarchy中的数组的下标.
	int lchild_node_handle_, rchild_node_handle_;
	int fund_cut_node_handle_[2];//(from, to)-> to, (v0, v1)->v1时候的vl^, vr^这两个原始点(在original mesh上的)所对应的Node序号.

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
	std::vector<VHierarchyNode> nodes_; //这个数组的下标作为node_handle
	int n_leaves_;//the leaves' num of this vertex hierarchy

public:
	void clear() { nodes_.clear(); } //清空整个森林, 每次要重新建立森林前都该清空.

	int add_node() { //nodes_数组的下标的是0到nodes_.size()-1.
		nodes_.push_back(VHierarchyNode());
		nodes_[nodes_.size() - 1].set_self_node_handle(nodes_.size() - 1); 
		return nodes_.size() - 1;
	}
	VHierarchyNode& node(int _node_handle) {
		return nodes_[_node_handle];
	}
	int size() {		// 森林大小
		return nodes_.size();
	}
	int n_leaves() {	// 叶子多少
		return n_leaves_;
	}
	void set_n_leaves(int _n_leaves) { //这个函数只是调用一次的.
		n_leaves_ = _n_leaves;//就是初始为original mesh上的所有顶点都分配一个Node的时候.
	}
};
} // end of namespace DGP
#endif //dgpstudio_vhierarchy_h
