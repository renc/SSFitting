
//-----------------------
// MeshModel

#include "MeshModel.h"
//#include "MainWindow.h" // renc qt version 
//#include "GLWidget.h" // renc qt version 
#include "GUIFrame.h"

MeshModel::MeshModel()
:is_modified(false), draw_mode_(SOLIDSMOOTH), draw_coordiante_system_(true) {

	mesh_toberendered_ = new TriMesh;
	//mesh_.request_face_normals();//��Ϊcompile time����Attributes��������
	//mesh_.request_vertex_normals();//��Ϊcompile time����Attributes��������	

	line_color_ = OpenMesh::Vec3f(0.f, 0.f, 0.f);//�߿���ɫ��ʼ��Ϊ��ɫ
}
//void MeshModel::add_guiframe(MainWindow *_guiframe) { guiframe_ = _guiframe; }
void MeshModel::add_guiframe(GUIFrame *_guiframe) { guiframe_ = _guiframe; }
//////////////////////////////////////////////////////////////////
// ���涼�Ǵ�CourseExamples�е�MeshViewer�п���������,
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
			bb.bbMin.minimize(mesh_.point(v_it));//ȡ����Ԫ�ص���С��ֵ
			bb.bbMax.maximize(mesh_.point(v_it));//ȡ����Ԫ�ص�����ֵ
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
		draw_mode_ = SOLIDSMOOTH;//ÿ������ģ�ͺ󶼽���ʾģʽ�������û�

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
	{	// �������е���, ��¼��ÿһ����Ķ������
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
		// ��һ��������
		glPushMatrix();
		glBegin(GL_LINES);
		glColor3f(1, 0, 0);//ע�͵���Щ��ɫ����Ϊ��Ӱ�������������ɫ����
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
	
	if (indices_.empty()) {//û�ж�����Ƭ��ʱ��,�����丸�����ʾ����
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
	} else if (draw_mode_ == WIREFRAME) { // �����model����Ƭʱ��, ����draw mode
		glDisable(GL_LIGHTING);
		GL::glColor(line_color_);
		glPolygonMode(GL_FRONT_AND_BACK, GL_LINE);//�����߿�ģʽ

		// ����ԭ������Ⱦ��ʽ,ʹ�ö�������. ��������DecimationModel��ģ�����򻯺�ͷ������ַ�ʽ������ȷ��Ⱦģ����,�������಻��ȷ���߱�
		// ����������û�и���update_renderface_indices().
		glEnableClientState(GL_VERTEX_ARRAY);  
		GL::glVertexPointer(mesh_toberendered_->points());

		glDrawElements(GL_TRIANGLES, indices_.size(), GL_UNSIGNED_INT, &indices_[0]);

		glDisableClientState(GL_VERTEX_ARRAY);

		/*// ����������ȡ���������Ⱦ��ʽ��, ��ʵ������������solid flat����Ⱦ����
		TriMesh::ConstFaceIter f_it(mesh_.faces_begin()), f_end(mesh_.faces_end());
		TriMesh::ConstFaceVertexIter  fv_it;

		glBegin(GL_TRIANGLES);
		for (; f_it!=f_end; ++f_it) { // ���ַ�ʽӦ��ѧѧ��
			GL::glNormal(mesh_.normal(f_it));
			fv_it = mesh_.cfv_iter(f_it.handle()); 
			GL::glVertex(mesh_.point(fv_it));
			++fv_it;
			GL::glVertex(mesh_.point(fv_it));
			++fv_it;
			GL::glVertex(mesh_.point(fv_it));
		}
		glEnd();*/

		glPolygonMode(GL_FRONT_AND_BACK, GL_FILL);//�����������ģʽ

	} else if (draw_mode_ == HIDDENLINE) {
		glDisable(GL_LIGHTING);
		GL::glColor(line_color_);
		
		glDrawBuffer(GL_NONE);
		glDepthRange(0.01, 1.0);
		TriMesh::ConstFaceIter f_it(mesh_toberendered_->faces_begin()), f_end(mesh_toberendered_->faces_end());
		TriMesh::ConstFaceVertexIter  fv_it;

		glBegin(GL_TRIANGLES);
		for (; f_it!=f_end; ++f_it) { // ���ַ�ʽӦ��ѧѧ��
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
		for (f_it = mesh_toberendered_->faces_begin(); f_it!=f_end; ++f_it) { // ���ַ�ʽӦ��ѧѧ��
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
		glShadeModel(GL_FLAT);//һ����(�������еĶ���)ʹ��һ�����߷���

		TriMesh::ConstFaceIter f_it(mesh_toberendered_->faces_begin()), f_end(mesh_toberendered_->faces_end());
		TriMesh::ConstFaceVertexIter  fv_it;

		glBegin(GL_TRIANGLES);
		for (; f_it!=f_end; ++f_it) { // ���ַ�ʽӦ��ѧѧ��
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
		glShadeModel(GL_SMOOTH);//һ����Ķ���ʹ�ø��Եķ��߷���

		// ����ԭ������Ⱦ��ʽ,ʹ�ö�������. ��������DecimationModel��ģ�����򻯺�ͷ������ַ�ʽ������ȷ��Ⱦģ����,�������಻��ȷ���߱�
		glEnableClientState(GL_VERTEX_ARRAY); // ���ַ�ʽӦ��ѧѧ��
		glEnableClientState(GL_NORMAL_ARRAY);
		GL::glVertexPointer(mesh_toberendered_->points());
		GL::glNormalPointer(mesh_toberendered_->vertex_normals());

		glDrawElements(GL_TRIANGLES, indices_.size(), GL_UNSIGNED_INT, &indices_[0]);

		glDisableClientState(GL_VERTEX_ARRAY);
		glDisableClientState(GL_NORMAL_ARRAY);

		/*// ����������ȡ���������Ⱦ��ʽ��, ��ʵ������������solid flat����Ⱦ����,�������޸�
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
		// ��wire frame��ʽ��ʾ������
		glDisable(GL_LIGHTING);
		GL::glColor(line_color_);//��ɫ�߿�
		glPolygonMode(GL_FRONT_AND_BACK, GL_LINE);//�����߿�ģʽ

		// ����ԭ������Ⱦ��ʽ,ʹ�ö�������. ��������DecimationModel��ģ�����򻯺�ͷ������ַ�ʽ������ȷ��Ⱦģ����,�������಻��ȷ���߱�
		// ����������û�и���update_renderface_indices().
		glEnableClientState(GL_VERTEX_ARRAY);  
		GL::glVertexPointer(mesh_toberendered_->points());

		glDrawElements(GL_TRIANGLES, indices_.size(), GL_UNSIGNED_INT, &indices_[0]);

		glDisableClientState(GL_VERTEX_ARRAY);
		glPolygonMode(GL_FRONT_AND_BACK, GL_FILL);//�����������ģʽ

		// ��"Solid Flat"��ʽ��ʾ������
		TriMesh::ConstFaceIter f_it(mesh_toberendered_->faces_begin()), f_end(mesh_toberendered_->faces_end());
		TriMesh::ConstFaceVertexIter  fv_it;

		glEnable(GL_LIGHTING);
		glShadeModel(GL_FLAT);

		glBegin(GL_TRIANGLES);
		for (; f_it!=f_end; ++f_it) { // ���ַ�ʽӦ��ѧѧ��
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
