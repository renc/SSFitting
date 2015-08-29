#include "TrackballManipulator.h"


namespace DGP {
//----------------------------------------------------------------------------

	void TrackballManipulator::init(int _w, int _h) {
		width_ = _w;
		height_ = _h;
		glMatrixMode(GL_MODELVIEW);
		glLoadIdentity();
		glGetDoublev(GL_MODELVIEW_MATRIX, modelview_matrix_);
		set_scene(OpenMesh::Vec3f(0.0, 0.0, 0.0), 1.0);//center, radius
	}
	
void TrackballManipulator::update_reshape(int _w, int _h) {
	width_ = _w;
	height_ = _h;

	update_projection_matrix();
}


//----------------------------------------------------------------------------
// public, 
// (1) 被init函数所调用, 
// (2) and invoked right after load a new model
// 在reshape时候即使width and height也不调用这个set_scene, 所以只是被set_scene调用的view_all也就更不会被调用了.
void TrackballManipulator::set_scene( const OpenMesh::Vec3f& _cog, float _radius )
{
	center_ = _cog;//(0,0,0)
	radius_ = _radius;//1

	near_  = 0.01f * radius_;
	far_   = 100.0f * radius_;
	fovy_ = 45.0;
	update_projection_matrix();

	view_all();
}
void TrackballManipulator::mouse(int button, int state, int x, int y) {
	// mouse press
	if (state == 0)
	{
		last_point_2D_ = OpenMesh::Vec2i(x,y);
		last_point_ok_ = map_to_sphere( last_point_2D_, last_point_3D_ );
		button_down_[button] = true;
	}
	// mouse release
	else
	{
		last_point_ok_ = false;
		button_down_[button] = false;

		//// GLUT: button 3 or 4 -> mouse wheel clicked
		//if (button == 3)       
		//	zoom(0, (int)(last_point_2D_[1] - 0.05*width_));
		//else if (button == 4)
		//	zoom(0, (int)(last_point_2D_[1] + 0.05*width_));
	}

}
void TrackballManipulator::motion(int x, int y) {
	// zoom // 由于右键用于绑定了菜单所以用一二两键模拟第三键
	if (button_down_[0] && button_down_[1])
	{
		zoom(x, y);
	}
	// rotation
	else if (button_down_[0])
	{
		rotation(x, y);
	}
	// translation
	else if (button_down_[1])
	{
		translation(x, y);
	}

	// remeber points
	last_point_2D_ = OpenMesh::Vec2i(x, y);
	last_point_ok_ = map_to_sphere(last_point_2D_, last_point_3D_);
}

//----------------------------------------------------------------------------
// 只是被set_scene函数调用
void TrackballManipulator::view_all()
{  // modelview_matrix_ * center_
	OpenMesh::Vec3f t = OpenMesh::Vec3f(   
		-(modelview_matrix_[0]*center_[0] + 
		modelview_matrix_[4]*center_[1] +
		modelview_matrix_[8]*center_[2] + 
		modelview_matrix_[12] + 0*1),
		-(modelview_matrix_[1]*center_[0] + 
		modelview_matrix_[5]*center_[1] +
		modelview_matrix_[9]*center_[2] + 
		modelview_matrix_[13] + 0*1),
		-(modelview_matrix_[2]*center_[0] + 
		modelview_matrix_[6]*center_[1] +
		modelview_matrix_[10]*center_[2] + 
		modelview_matrix_[14] + 3.0f*radius_*1) 
		);
	std::cout << "For test: t is " << t << " in the view_all." << modelview_matrix_[12] << ", " << modelview_matrix_[13] << ", " << modelview_matrix_[14] << ". radius " << radius_ << ".\n";
	translate( t ); //z方向上移动的负量, 系数3改成更大的话(例如6.0)看上去物体更小,因为离视点更远了
}
// protected, 只被reshape函数和set_scene函数调用
void TrackballManipulator::update_projection_matrix()
{ 
	glMatrixMode( GL_PROJECTION );
	glLoadIdentity();
	gluPerspective(fovy_, (GLfloat)width_/(GLfloat)height_, near_, far_);//这是glut下和reshape函数一致的
	glGetDoublev( GL_PROJECTION_MATRIX, projection_matrix_);
	glMatrixMode( GL_MODELVIEW );
}
//----------------------------------------------------------------------------
// protected, used in virtual traceball,
bool TrackballManipulator::map_to_sphere( const OpenMesh::Vec2i& _v2D, OpenMesh::Vec3f& _v3D )
{ 
	if ((_v2D[0] >= 0) && (_v2D[0] <= width_) &&
		(_v2D[1] >= 0) && (_v2D[1] <= height_)) 
	{
		double x  = (double)(_v2D[0] - 0.5*width_)  / (double)width_;
		double y  = (double)(0.5*height_ - _v2D[1]) / (double)height_;
		double sinx         = sin(M_PI * x * 0.5);
		double siny         = sin(M_PI * y * 0.5);
		double sinx2siny2   = sinx * sinx + siny * siny;

		_v3D[0] = sinx;
		_v3D[1] = siny;
		_v3D[2] = sinx2siny2 < 1.0f ? sqrt(1.0 - sinx2siny2) : 0.0f;

		return true;
	}
	else return false;
}

//-----------------------------------------------------------------------------
void TrackballManipulator::rotation(int x, int y)
{
	if (last_point_ok_) 
	{
		OpenMesh::Vec3f  new_point_3D;

		OpenMesh::Vec2i  new_point_2D = OpenMesh::Vec2i(x, y);
		bool   new_point_ok = map_to_sphere(new_point_2D, new_point_3D);

		if (new_point_ok)
		{
			OpenMesh::Vec3f axis      = (last_point_3D_ % new_point_3D);
			float cos_angle = (last_point_3D_ | new_point_3D);

			if (fabs(cos_angle) < 1.0) 
			{
				rotate(axis, 2.0f*acos(cos_angle) * 180.0f / M_PI);
			}
		}
	}
}

//-----------------------------------------------------------------------------
void TrackballManipulator::translation(int x, int y)
{ // x, y are the current mouse location.
	float dx = x - last_point_2D_[0];
	float dy = y - last_point_2D_[1];

	float z = - ((modelview_matrix_[ 2]*center_[0] + 
		modelview_matrix_[ 6]*center_[1] + 
		modelview_matrix_[10]*center_[2] + 
		modelview_matrix_[14]) /
		(modelview_matrix_[ 3]*center_[0] + 
		modelview_matrix_[ 7]*center_[1] + 
		modelview_matrix_[11]*center_[2] + 
		modelview_matrix_[15]));

	float aspect = (float)width_ / (float)height_;
	float up     = tan(fovy_/2.0f*M_PI/180.f) * near_;
	float right  = aspect*up;

	translate(OpenMesh::Vec3f(
		 2.0f*dx/width_*right/near_*z, 
		-2.0f*dy/height_*up/near_*z, 
		 0.0f));
}


//-----------------------------------------------------------------------------
void TrackballManipulator::zoom(int x, int y)
{
	float dy = y - last_point_2D_[1];
	float h  = height_;
	translate(OpenMesh::Vec3f(0.0f, 0.0f, radius_ * dy * 3.0f / h));
}

//----------------------------------------------------------------------------
// support the above rotation and translation function.
void TrackballManipulator::translate( const OpenMesh::Vec3f& _trans )
{
	glLoadIdentity();
	glTranslated( _trans[0], _trans[1], _trans[2] );
	glMultMatrixd( modelview_matrix_ );
	glGetDoublev( GL_MODELVIEW_MATRIX, modelview_matrix_);
}

//----------------------------------------------------------------------------
void TrackballManipulator::rotate( const OpenMesh::Vec3f& _axis, float _angle )
{
	OpenMesh::Vec3d t( modelview_matrix_[0]*center_[0] + 
		modelview_matrix_[4]*center_[1] +
		modelview_matrix_[8]*center_[2] + 
		modelview_matrix_[12],
		modelview_matrix_[1]*center_[0] + 
		modelview_matrix_[5]*center_[1] +
		modelview_matrix_[9]*center_[2] + 
		modelview_matrix_[13],
		modelview_matrix_[2]*center_[0] + 
		modelview_matrix_[6]*center_[1] +
		modelview_matrix_[10]*center_[2] + 
		modelview_matrix_[14] );

	glLoadIdentity();
	glTranslated(t[0], t[1], t[2]);
	glRotated( _angle, _axis[0], _axis[1], _axis[2]);
	glTranslated(-t[0], -t[1], -t[2]); 
	glMultMatrixd(modelview_matrix_);
	glGetDoublev(GL_MODELVIEW_MATRIX, modelview_matrix_);
}
OpenMesh::Vec3d TrackballManipulator::object_to_eye_coordiante(const OpenMesh::Vec3d &_p) {
	return OpenMesh::Vec3d ( modelview_matrix_[0]*_p[0] + 
		modelview_matrix_[4]*_p[1] +
		modelview_matrix_[8]*_p[2] + 
		modelview_matrix_[12],                   
		modelview_matrix_[1]*_p[0] + 
		modelview_matrix_[5]*_p[1] +
		modelview_matrix_[9]*_p[2] + 
		modelview_matrix_[13],                
		modelview_matrix_[2]*_p[0] + 
		modelview_matrix_[6]*_p[1] +
		modelview_matrix_[10]*_p[2] + 
		modelview_matrix_[14] ); //本来xyz每一项都分别需要除以w的，但是这里没有scale成分, so w=1.
}
} // end of namespace DGP 