
//-----------------------
// MeshModel

#include "MeshModel.h"
//#include "MainWindow.h" // renc qt version 
//#include "GLWidget.h" // renc qt version 
#include "GUIFrame.h"

MeshModel::MeshModel()
:is_modified(false), draw_mode_(SOLIDSMOOTH), draw_coordiante_system_(true) {

	mesh_toberendered_ = new TriMesh;
	//mesh_.request_face_normals();//改为compile time是用Attributes来申请了
	//mesh_.request_vertex_normals();//改为compile time是用Attributes来申请了	

	line_color_ = OpenMesh::Vec3f(0.f, 0.f, 0.f);//线框颜色初始化为白色
}
//void MeshModel::add_guiframe(MainWindow *_guiframe) { guiframe_ = _guiframe; }
void MeshModel::add_guiframe(GUIFrame *_guiframe) { guiframe_ = _guiframe; }
//////////////////////////////////////////////////////////////////
// 下面都是从CourseExamples中的MeshViewer中拷贝过来的,
//-----------------------------------------------------------------------------
bool MeshModel::open_mesh(const char* _filename)
{
	//if (mesh_.empty() == false) { // renc: 20150829 
	if (mesh_.n_vertices()) {
		mesh_.clear(); 
	} 
	// load mesh
	if (OpenMesh::IO::read_mesh(mesh_, _filename)) {
		// set center and radius
		TriMesh::ConstVertexIter  v_it(mesh_.vertices_begin()), v_end(mesh_.vertices_end());
		//TriMesh::Point            bbMin, bbMax;

		bb.bbMin = bb.bbMax = mesh_.point(v_it);
		for (; v_it!=v_end; ++v_it) {
			bb.bbMin.minimize(mesh_.point(v_it));//取各个元素的最小的值
			bb.bbMax.maximize(mesh_.point(v_it));//取各个元素的最大的值
		}
		bb.center = (bb.bbMin + bb.bbMax) * 0.5;
		bb.diagonal_lenght = (bb.bbMin - bb.bbMax).norm();
		//guiframe_->glWidget()->set_scene( (OpenMesh::Vec3f)(bb.center), 0.5 * bb.diagonal_lenght);//center, radius // renc 20150829
		guiframe_->set_scene( (OpenMesh::Vec3f)(bb.center), 0.5 * bb.diagonal_lenght);//center, radius

		update();
		// info
		std::cerr << mesh_.n_vertices() << " vertices, " << mesh_.n_edges() << " edges, " << mesh_.n_faces()    << " faces\n";
		
		mesh_toberendered_ = &mesh_;
		update_renderface_indices();
		draw_mode_ = SOLIDSMOOTH;//每次载入模型后都将显示模式重新设置回

		filename_ = std::string(_filename);
		return true;
	}

	return false;
}

//-----------------------------------------------------------------------------
// invoked by open_mesh(const char* _filename);in orther to fasten the render
void MeshModel::update_renderface_indices()// update face indices for faster rendering
{
	TriMesh::ConstFaceIter        f_it(mesh_toberendered_->faces_sbegin()), f_end(mesh_toberendered_->faces_end());
	TriMesh::ConstFaceVertexIter  fv_it;

	indices_.clear();
	indices_.reserve(mesh_toberendered_->n_faces()*3);

	for (; f_it!=f_end; ++f_it)
	{	// 遍历所有的面, 记录下每一个面的顶点序号
		indices_.push_back((fv_it=mesh_toberendered_->cfv_iter(f_it)).handle().idx());
		indices_.push_back((++fv_it).handle().idx());
		indices_.push_back((++fv_it).handle().idx());
	}
}

//-----------------------------------------------------------------------------
void MeshModel::draw()
{	
	if (draw_coordiante_system_) {
		glDisable(GL_LIGHTING);
		// 画一个坐标轴
		glPushMatrix();
		glBegin(GL_LINES);
		glColor3f(1, 0, 0);//注释掉这些颜色是因为会影响后面的物体的颜色属性
		glVertex3f(0, 0, 0);
		glVertex3f(1, 0, 0);
		glColor3f(0, 1, 0);
		glVertex3f(0, 0, 0);
		glVertex3f(0, 1, 0);
		glColor3f(0, 0, 1);
		glVertex3f(0, 0, 0);
		glVertex3f(0, 0, 1);
		glEnd();
		glPopMatrix();
		glEnable(GL_LIGHTING);
	}
	
	if (indices_.empty()) {//没有读入面片的时候,调用其父类的显示函数
		int draw_mode_ = SOLIDSMOOTH;
		if (draw_mode_ == WIREFRAME) {
			glDisable(GL_LIGHTING);
			glutWireTeapot(0.5);
		} else if (draw_mode_ == HIDDENLINE) {
			glDisable(GL_LIGHTING);
			glutWireTeapot(0.5);
		} else if (draw_mode_ == SOLIDFLAT) {
			glEnable(GL_LIGHTING);
			glShadeModel(GL_FLAT);
			glutSolidTeapot(0.5);
		} else { //if (draw_mode_ == SOLIDSMOOTH) {
			glEnable(GL_LIGHTING);
			glShadeModel(GL_SMOOTH);
			glutSolidTeapot(0.5);
		} 
	} else if (draw_mode_ == WIREFRAME) { // 读入的model有面片时候, 各种draw mode
		glDisable(GL_LIGHTING);
		GL::glColor(line_color_);
		glPolygonMode(GL_FRONT_AND_BACK, GL_LINE);//启动线框模式

		// 这是原来的渲染方式,使用顶点数组. 但是我在DecimationModel对模型做简化后就发现这种方式不能正确渲染模型了,会出现许多不正确的线边
		// 后来发现是没有更新update_renderface_indices().
		glEnableClientState(GL_VERTEX_ARRAY);  
		GL::glVertexPointer(mesh_toberendered_->points());

		glDrawElements(GL_TRIANGLES, indices_.size(), GL_UNSIGNED_INT, &indices_[0]);

		glDisableClientState(GL_VERTEX_ARRAY);

		/*// 这是我用来取代上面的渲染方式的, 其实就是拿了下面solid flat的渲染方法
		TriMesh::ConstFaceIter f_it(mesh_.faces_begin()), f_end(mesh_.faces_end());
		TriMesh::ConstFaceVertexIter  fv_it;

		glBegin(GL_TRIANGLES);
		for (; f_it!=f_end; ++f_it) { // 这种方式应该学学啊
			GL::glNormal(mesh_.normal(f_it));
			fv_it = mesh_.cfv_iter(f_it.handle()); 
			GL::glVertex(mesh_.point(fv_it));
			++fv_it;
			GL::glVertex(mesh_.point(fv_it));
			++fv_it;
			GL::glVertex(mesh_.point(fv_it));
		}
		glEnd();*/

		glPolygonMode(GL_FRONT_AND_BACK, GL_FILL);//重新启动填充模式

	} else if (draw_mode_ == HIDDENLINE) {
		glDisable(GL_LIGHTING);
		GL::glColor(line_color_);
		
		glDrawBuffer(GL_NONE);
		glDepthRange(0.01, 1.0);
		TriMesh::ConstFaceIter f_it(mesh_toberendered_->faces_begin()), f_end(mesh_toberendered_->faces_end());
		TriMesh::ConstFaceVertexIter  fv_it;

		glBegin(GL_TRIANGLES);
		for (; f_it!=f_end; ++f_it) { // 这种方式应该学学啊
			fv_it = mesh_toberendered_->cfv_iter(f_it.handle()); 
			GL::glVertex(mesh_toberendered_->point(fv_it));
			++fv_it;
			GL::glVertex(mesh_toberendered_->point(fv_it));
			++fv_it;
			GL::glVertex(mesh_toberendered_->point(fv_it));
		}
		glEnd();

		glDrawBuffer(GL_BACK);
		glDepthRange(0.0, 1.0);
		glPolygonMode(GL_FRONT_AND_BACK, GL_LINE);
		glDepthFunc(GL_LEQUAL);
		
		glBegin(GL_TRIANGLES);
		for (f_it = mesh_toberendered_->faces_begin(); f_it!=f_end; ++f_it) { // 这种方式应该学学啊
			fv_it = mesh_toberendered_->cfv_iter(f_it.handle()); 
			GL::glVertex(mesh_toberendered_->point(fv_it));
			++fv_it;
			GL::glVertex(mesh_toberendered_->point(fv_it));
			++fv_it;
			GL::glVertex(mesh_toberendered_->point(fv_it));
		}
		glEnd();

		glPolygonMode(GL_FRONT_AND_BACK, GL_FILL);
		glDepthFunc(GL_LESS);

	} else if (draw_mode_ == SOLIDFLAT) {
		glEnable(GL_LIGHTING);
		glShadeModel(GL_FLAT);//一个面(其中所有的顶点)使用一个法线方向

		TriMesh::ConstFaceIter f_it(mesh_toberendered_->faces_begin()), f_end(mesh_toberendered_->faces_end());
		TriMesh::ConstFaceVertexIter  fv_it;

		glBegin(GL_TRIANGLES);
		for (; f_it!=f_end; ++f_it) { // 这种方式应该学学啊
			GL::glNormal(mesh_toberendered_->normal(f_it));
			fv_it = mesh_toberendered_->cfv_iter(f_it.handle()); 
			GL::glVertex(mesh_toberendered_->point(fv_it));//std::cout << fv_it.handle().idx() << "\t";
			++fv_it;
			GL::glVertex(mesh_toberendered_->point(fv_it));//std::cout << fv_it.handle().idx() << "\t";
			++fv_it;
			GL::glVertex(mesh_toberendered_->point(fv_it));//std::cout << fv_it.handle().idx() << "\n";
		}
		glEnd();
	} else if (draw_mode_ == SOLIDSMOOTH) {
		glEnable(GL_LIGHTING);
		glShadeModel(GL_SMOOTH);//一个面的顶点使用各自的法线方向

		// 这是原来的渲染方式,使用顶点数组. 但是我在DecimationModel对模型做简化后就发现这种方式不能正确渲染模型了,会出现许多不正确的线边
		glEnableClientState(GL_VERTEX_ARRAY); // 这种方式应该学学啊
		glEnableClientState(GL_NORMAL_ARRAY);
		GL::glVertexPointer(mesh_toberendered_->points());
		GL::glNormalPointer(mesh_toberendered_->vertex_normals());

		glDrawElements(GL_TRIANGLES, indices_.size(), GL_UNSIGNED_INT, &indices_[0]);

		glDisableClientState(GL_VERTEX_ARRAY);
		glDisableClientState(GL_NORMAL_ARRAY);

		/*// 这是我用来取代上面的渲染方式的, 其实就是拿了上面solid flat的渲染方法,并稍作修改
		TriMesh::ConstFaceIter f_it(mesh_.faces_begin()), f_end(mesh_.faces_end());
		TriMesh::ConstFaceVertexIter  fv_it;

		glBegin(GL_TRIANGLES);
		for (; f_it!=f_end; ++f_it) {
			fv_it = mesh_.cfv_iter(f_it.handle()); 
			GL::glNormal(mesh_.normal(fv_it));
			GL::glVertex(mesh_.point(fv_it));
			++fv_it;
			GL::glNormal(mesh_.normal(fv_it));
			GL::glVertex(mesh_.point(fv_it));
			++fv_it;
			GL::glNormal(mesh_.normal(fv_it));
			GL::glVertex(mesh_.point(fv_it));
		}
		glEnd();*/
	} else if (draw_mode_ == LINEANDSOLIDFLAT) {
		// 以wire frame形式显示三角形
		glDisable(GL_LIGHTING);
		GL::glColor(line_color_);//白色线框
		glPolygonMode(GL_FRONT_AND_BACK, GL_LINE);//启动线框模式

		// 这是原来的渲染方式,使用顶点数组. 但是我在DecimationModel对模型做简化后就发现这种方式不能正确渲染模型了,会出现许多不正确的线边
		// 后来发现是没有更新update_renderface_indices().
		glEnableClientState(GL_VERTEX_ARRAY);  
		GL::glVertexPointer(mesh_toberendered_->points());

		glDrawElements(GL_TRIANGLES, indices_.size(), GL_UNSIGNED_INT, &indices_[0]);

		glDisableClientState(GL_VERTEX_ARRAY);
		glPolygonMode(GL_FRONT_AND_BACK, GL_FILL);//重新启动填充模式

		// 以"Solid Flat"形式显示三角形
		TriMesh::ConstFaceIter f_it(mesh_toberendered_->faces_begin()), f_end(mesh_toberendered_->faces_end());
		TriMesh::ConstFaceVertexIter  fv_it;

		glEnable(GL_LIGHTING);
		glShadeModel(GL_FLAT);

		glBegin(GL_TRIANGLES);
		for (; f_it!=f_end; ++f_it) { // 这种方式应该学学啊
			GL::glNormal(mesh_toberendered_->normal(f_it));
			fv_it = mesh_toberendered_->cfv_iter(f_it.handle()); 
			GL::glVertex(mesh_toberendered_->point(fv_it));
			++fv_it;
			GL::glVertex(mesh_toberendered_->point(fv_it));
			++fv_it;
			GL::glVertex(mesh_toberendered_->point(fv_it));
		}
		glEnd();
	}
}

//////////////////////////////////////////////////////////////////
void MeshModel::update() {
	// compute face & vertex normals
	mesh_.update_normals();	
	
	this->is_modified = true;
}
/////////////////////////////////////////////////////////////////
