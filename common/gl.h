
#ifndef dgp_common_gl_h
#define dgp_common_gl_h

#include "gl.hh"
#include "OpenMeshMeshType.h"
namespace DGP {
	
	struct Color { // use the OpenMesh::Vec3f not OpenMesh::Vec3d.
		static OpenMesh::Vec3f Black() { return OpenMesh::Vec3f(0, 0, 0); }
		static OpenMesh::Vec3f White() { return OpenMesh::Vec3f(1, 1, 1); }
		static OpenMesh::Vec3f Red() { return OpenMesh::Vec3f(1, 0, 0); }
		static OpenMesh::Vec3f Green() { return OpenMesh::Vec3f(0, 1, 0); }
		static OpenMesh::Vec3f Blue() { return OpenMesh::Vec3f(0, 0, 1); }
		static OpenMesh::Vec3f Yellow() { return OpenMesh::Vec3f(1, 1, 0); }
		static OpenMesh::Vec3f Brown() { return OpenMesh::Vec3f(121.0f/255, 61.0f/255, 0);}
	};

// Note: inline is necessary, or the function is redefined
inline void render_coordiante() {
	// 画一个坐标轴
	glDisable(GL_LIGHTING);
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
//-------------------------------------------------
// Wireframe 的draw mode,
// usage for example, DGP::wireframe_render(mesh_toberendered_->points(), findices_, DGP::Color::Blue());
inline void wireframe_render(const TriMesh::Point * _points, const unsigned int *_fidx, std::vector<unsigned int>::size_type _fidx_size,
							  const OpenMesh::Vec3f &_color) 
{
	glDisable(GL_LIGHTING);
	GL::glColor(_color);
	glPolygonMode(GL_FRONT_AND_BACK, GL_LINE);
	glEnableClientState(GL_VERTEX_ARRAY);
	GL::glVertexPointer(_points);
 
	glDrawElements(GL_TRIANGLES, _fidx_size, GL_UNSIGNED_INT, _fidx);
 
	glDisableClientState(GL_VERTEX_ARRAY);
	glPolygonMode(GL_FRONT_AND_BACK, GL_FILL); 
}
inline void wireframe_render(const TriMesh::Point * _points, const std::vector<unsigned int> &_fidx, const OpenMesh::Vec3f &_color) {
	wireframe_render(_points, &_fidx[0], _fidx.size(), _color);
}
inline void wireframe_render(const TriMesh &_mesh, const OpenMesh::Vec3f & _color, const float _line_width = 1) {
	glDisable(GL_LIGHTING);
	glColor3fv(_color.data());
	glLineWidth(_line_width);
	glPolygonMode(GL_FRONT_AND_BACK, GL_LINE);

	TriMesh::ConstFaceIter f_it(_mesh.faces_begin()), f_end(_mesh.faces_end());
	TriMesh::ConstFaceVertexIter fv_it;

	glBegin(GL_TRIANGLES);
	for (; f_it != f_end; ++f_it) {
		GL::glNormal(_mesh.normal(f_it));
		fv_it = _mesh.cfv_iter(f_it.handle());
		GL::glVertex(_mesh.point(fv_it));
		++fv_it;
		GL::glVertex(_mesh.point(fv_it));
		++fv_it;
		GL::glVertex(_mesh.point(fv_it));
	}
	glEnd();
	glLineWidth(1);
	glPolygonMode(GL_FRONT_AND_BACK, GL_FILL);
}
inline void wireframe_render_poly(const PolyMesh &_mesh, const OpenMesh::Vec3f & _color = DGP::Color::Red(), const float _line_width = 1) {
	glDisable(GL_LIGHTING);
	glColor3fv(_color.data());
	glLineWidth(_line_width);
	glBegin(GL_LINES);
	for (PolyMesh::CEIter e_it(_mesh.edges_begin()), e_end(_mesh.edges_end()); e_it != e_end; ++e_it) {
		GL::glVertex(_mesh.point(_mesh.to_vertex_handle(_mesh.halfedge_handle(e_it, 0))));
		GL::glVertex(_mesh.point(_mesh.to_vertex_handle(_mesh.halfedge_handle(e_it, 1)))); 
	}
	glEnd();
	glLineWidth(1);
}
// -----------------------------------------------------------------------------
// Hidder line draw mode, 尽量使用第二/三个.
// usage for example, DGP::hidderline_render(mesh_toberendered_->points(), findices_, DGP::Color::Blue());
inline void hidderline_render(const TriMesh &_mesh, const OpenMesh::Vec3f &_color) {
	glDisable(GL_LIGHTING);  //这个函数好像有点问题，背面的显示有点问题。
	GL::glColor(_color);// 网格的颜色 

	glDrawBuffer(GL_NONE);
	glDepthRange(0.01, 1.0);
	TriMesh::ConstFaceIter f_it(_mesh.faces_begin()), f_end(_mesh.faces_end());
	TriMesh::ConstFaceVertexIter  fv_it;

	glBegin(GL_TRIANGLES);
	for (; f_it!=f_end; ++f_it) { // 这种方式应该学学啊
		fv_it = _mesh.cfv_iter(f_it.handle()); 
		GL::glVertex(_mesh.point(fv_it));
		++fv_it;
		GL::glVertex(_mesh.point(fv_it));
		++fv_it;
		GL::glVertex(_mesh.point(fv_it));
	}
	glEnd();

	glDrawBuffer(GL_BACK);
	glDepthRange(0.0, 1.0);
	glPolygonMode(GL_FRONT_AND_BACK, GL_LINE);
	glDepthFunc(GL_LEQUAL);

	glBegin(GL_TRIANGLES);
	for (f_it = _mesh.faces_begin(); f_it!=f_end; ++f_it) { // 这种方式应该学学啊
		fv_it = _mesh.cfv_iter(f_it.handle()); 
		GL::glVertex(_mesh.point(fv_it));
		++fv_it;
		GL::glVertex(_mesh.point(fv_it));
		++fv_it;
		GL::glVertex(_mesh.point(fv_it));					
	}
	glEnd();

	glPolygonMode(GL_FRONT_AND_BACK, GL_FILL);
	glDepthFunc(GL_LESS);
}
inline void hidderline_render(const TriMesh::Point * _points, const unsigned int *_fidx, std::vector<unsigned int>::size_type _fidx_size,
							  const OpenMesh::Vec3f &_color) 
{
	glDisable(GL_LIGHTING);
	GL::glColor(_color);
	glEnableClientState(GL_VERTEX_ARRAY);
	GL::glVertexPointer(_points);

	glDrawBuffer(GL_NONE);
	glDepthRange(0.01, 1.0);
	glDrawElements(GL_TRIANGLES, _fidx_size, GL_UNSIGNED_INT, _fidx);

	glDrawBuffer(GL_BACK);
	glDepthRange(0.0, 1.0);
	glPolygonMode(GL_FRONT_AND_BACK, GL_LINE);//glLineWidth(2.0);	
	glDepthFunc(GL_LEQUAL);
	glDrawElements(GL_TRIANGLES, _fidx_size, GL_UNSIGNED_INT, _fidx);//glLineWidth(1.0);	

	glDisableClientState(GL_VERTEX_ARRAY);
	glPolygonMode(GL_FRONT_AND_BACK, GL_FILL);
	glDepthFunc(GL_LESS); 
}
inline void hidderline_render(const std::vector<OpenMesh::Vec3d> &_points, const std::vector<unsigned int> &_fidx, 
							  const OpenMesh::Vec3f &_color) 
{
	hidderline_render(&_points[0], &_fidx[0], _fidx.size(), _color);
}
inline void hidderline_render_poly(const PolyMesh &_mesh, const OpenMesh::Vec3f &_color) {
	glDisable(GL_LIGHTING);   //这个函数好像有点问题，背面的显示有点问题。
	GL::glColor(_color);// 网格的颜色 

	glDrawBuffer(GL_NONE);
	glDepthRange(0.01, 1.0);
	
	PolyMesh::ConstFaceVertexIter  fv_it;

	for (PolyMesh::ConstFaceIter f_it(_mesh.faces_begin()), f_end(_mesh.faces_end()); f_it != f_end; ++f_it) { 
		glBegin(GL_POLYGON);
		for (fv_it = _mesh.cfv_iter(f_it.handle()); fv_it; ++fv_it) {
			GL::glVertex(_mesh.point(fv_it));
		}
		glEnd();
	}	 

	glDrawBuffer(GL_BACK);
	glDepthRange(0.0, 1.0);
	glPolygonMode(GL_FRONT_AND_BACK, GL_LINE);
	glDepthFunc(GL_LEQUAL);

	for (PolyMesh::ConstFaceIter f_it(_mesh.faces_begin()), f_end(_mesh.faces_end()); f_it != f_end; ++f_it) { 
		glBegin(GL_POLYGON);
		for (fv_it = _mesh.cfv_iter(f_it.handle()); fv_it; ++fv_it) {
			GL::glVertex(_mesh.point(fv_it));
		}
		glEnd();
	} 

	glPolygonMode(GL_FRONT_AND_BACK, GL_FILL);
	glDepthFunc(GL_LESS);
}
//-------------------------------------------------
// Solid Flat(flat shade model) draw mode.
// usage for example, DGP::solidflat_render(*mesh_toberendered_);
inline void solidflat_render(const TriMesh &_m) {
	TriMesh::ConstFaceIter f_it(_m.faces_begin()), f_end(_m.faces_end());
	TriMesh::ConstFaceVertexIter  fv_it;

	glEnable(GL_LIGHTING);
	glShadeModel(GL_FLAT);

	glBegin(GL_TRIANGLES);
	for (; f_it!=f_end; ++f_it) { // 这种方式应该学学啊
		GL::glNormal(_m.normal(f_it)); // 之前需要 _m.update_normals();
		fv_it = _m.cfv_iter(f_it.handle()); 
		GL::glVertex(_m.point(fv_it));
		++fv_it;
		GL::glVertex(_m.point(fv_it));
		++fv_it;
		GL::glVertex(_m.point(fv_it));
	}
	glEnd();
}
inline void solidflat_render_poly(const PolyMesh &_m) {
	glEnable(GL_LIGHTING);
	glShadeModel(GL_FLAT);

	PolyMesh::ConstFaceIter f_it(_m.faces_begin()), f_end(_m.faces_end());
	PolyMesh::ConstFaceVertexIter fv_it;

	for (; f_it != f_end; ++f_it) { 
		glBegin(GL_POLYGON); // GL_POLYGON draw only one polygon, this is not like the GL_TRIANGLES draw a serial triangles.
		GL::glNormal(_m.normal(f_it));
		for (fv_it = _m.cfv_iter(f_it.handle()); fv_it; ++fv_it) {
			GL::glVertex(_m.point(fv_it));
		}
		glEnd();
	}	
}
//-------------------------------------------------
// Solid Smooth(smooth shade model) draw mode.
// usage for example, DGP::solidsmooth_render(mesh_toberendered_->points(), mesh_toberendered_->vertex_normals(), findices_); 
inline void solidsmooth_render(const OpenMesh::Vec3d *_points, const OpenMesh::Vec3d* _vnormals, const unsigned int * _fidx, std::vector<unsigned int>::size_type _fidx_size) {
	glEnable(GL_LIGHTING);
	glShadeModel(GL_SMOOTH);

	glEnableClientState(GL_VERTEX_ARRAY); // 这种方式应该学学啊
	glEnableClientState(GL_NORMAL_ARRAY);
	GL::glVertexPointer(_points);
	GL::glNormalPointer(_vnormals);

	glDrawElements(GL_TRIANGLES, _fidx_size, GL_UNSIGNED_INT, _fidx);

	glDisableClientState(GL_VERTEX_ARRAY);
	glDisableClientState(GL_NORMAL_ARRAY); 
}
inline void solidsmooth_render(const TriMesh::Point * _points, const TriMesh::Normal * _vnormals, const std::vector<unsigned int> &_fidx) {
	solidsmooth_render(_points, _vnormals, &_fidx[0], _fidx.size());
}
inline void solidsmooth_render_poly(const PolyMesh & _m) {
	glEnable(GL_LIGHTING);
	glShadeModel(GL_SMOOTH);

	PolyMesh::ConstFaceIter f_it(_m.faces_begin()), f_end(_m.faces_end());
	PolyMesh::ConstFaceVertexIter fv_it;

	for (; f_it != f_end; ++f_it) { 
		glBegin(GL_POLYGON); // GL_POLYGON draw only one polygon, this is not like the GL_TRIANGLES draw a serial triangles.
		for (fv_it = _m.cfv_iter(f_it.handle()); fv_it; ++fv_it) {
			GL::glNormal(_m.normal(fv_it));
			GL::glVertex(_m.point(fv_it));
		}
		glEnd();
	}	
}

inline bool pick(int x, int y, OpenMesh::Vec3f& _p)
{
	GLdouble  modelview[16], projection[16];
	GLint     viewport[4];
	glGetDoublev(GL_MODELVIEW_MATRIX, modelview);
	glGetDoublev(GL_PROJECTION_MATRIX, projection);
	glGetIntegerv(GL_VIEWPORT, viewport);


	// read depth buffer value at (x, y_new)
	float  z;
	int    y_new = viewport[3] - y; // in OpenGL y is zero at the 'bottom'
	glReadPixels(x, y_new, 1, 1, GL_DEPTH_COMPONENT, GL_FLOAT, &z );


	// reverse projection to get 3D point
	double pos[3];
	gluUnProject(x, y_new, z, 
		modelview, 
		projection, 
		viewport, 
		&pos[0], &pos[1], &pos[2]);


	if (z != 1.0f)
	{
		_p = OpenMesh::Vec3f(pos[0], pos[1], pos[2]);
		return true;
	}

	return false;
}

//-------------------------------------------------
} // end of namespace DGP

#endif //dgp_common_gl_h