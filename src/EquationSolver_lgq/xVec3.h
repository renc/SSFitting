// xVec3.h: interface for the template class xVec3 .
//
//////////////////////////////////////////////////////////////////////

#ifndef __MESHOPTI_VEC3_H__
#define __MESHOPTI_VEC3_H__

#include <math.h>
#include <string.h>


template<class TYPE>
class xVec3
{
public:
	// constructor
	xVec3() { }
	~xVec3() { }
	xVec3(TYPE x, TYPE y, TYPE z) { m_x = x; m_y = y; m_z = z; }
	xVec3(const TYPE e[]) { m_x = e[0]; m_y = e[1]; m_z = e[2]; }
	xVec3(const xVec3& v) { m_x = v.m_x; m_y = v.m_y; m_z = v.m_z; }
	void Set(TYPE x, TYPE y, TYPE z) { m_x = x; m_y = y; m_z = z; }
	void Set(const TYPE e[]) { m_x = e[0]; m_y = e[1]; m_z = e[2]; }
	void Set(const xVec3& v) { m_x = v.m_x; m_y = v.m_y; m_z = v.m_z; }

public:
	template<class TYPE2> xVec3 & operator = (const xVec3<TYPE2>& v) { m_x = (TYPE)(v.m_x); m_y = (TYPE)(v.m_y); m_z = (TYPE)(v.m_z); return (*this); }
	// components
	TYPE&	operator[](int i) { return m_e[i]; }
	TYPE&	X() { return m_x; }
	TYPE&	Y() { return m_y; }
	TYPE&	Z() { return m_z; }

	TYPE operator[](int i) const { return m_e[i]; }
	TYPE X() const { return m_x; }
	TYPE Y() const { return m_y; }
	TYPE Z() const { return m_z; }

	// norm
	double Length() const { return sqrt(m_x*(double)m_x + m_y*(double)m_y + m_z*(double)m_z); }
	double Length2()const { return (m_x*(double)m_x + m_y*(double)m_y + m_z*(double)m_z); }
	// convert to a unit vector
	xVec3& Normalize();
	// type conversion
	operator TYPE*() { return m_e; }
	operator const TYPE*() const { return m_e; }

public:
	union {
		struct { TYPE m_x, m_y, m_z;};
		TYPE m_e[3];
	};
};

template<class TYPE>
inline xVec3<TYPE> & xVec3<TYPE>
::Normalize()
{
	double len = sqrt(m_x*(double)m_x+m_y*(double)m_y+m_z*(double)m_z);
//	if ( len<1.0e-20 )
//		return (*this);

	if ( len<1.0e-20 )
		return (*this);

	TYPE s = (TYPE)(1.0/len);
	m_x *= s; m_y *= s; m_z *= s;
	return (*this);
}

/* assignments */

/* a += b */
template<class TYPE1, class TYPE2>
inline xVec3<TYPE1>& operator += ( xVec3<TYPE1>& a, const xVec3<TYPE2>& b)
{
	a.m_x += (TYPE1)b.m_x; a.m_y += (TYPE1)b.m_y; a.m_z += (TYPE1)b.m_z;
	return a;
}

/* a -= b */
template<class TYPE1, class TYPE2>
inline xVec3<TYPE1>& operator -= ( xVec3<TYPE1>& a, const xVec3<TYPE2>& b)
{
	a.m_x -= (TYPE1)b.m_x; a.m_y -= (TYPE1)b.m_y; a.m_z -= (TYPE1)b.m_z;
	return a;
}

/* a *= s */
template<class TYPE1, class TYPE2>
inline xVec3<TYPE1>& operator *= ( xVec3<TYPE1>& a, TYPE2  s)
{
	a.m_x *= (TYPE1)s; a.m_y *= (TYPE1)s; a.m_z *= (TYPE1)s;
	return a;
}

/* a /= s */
template<class TYPE1, class TYPE2>
inline xVec3<TYPE1>& operator /= ( xVec3<TYPE1>& a, TYPE2  s)
{
	a.m_x = (TYPE1)(a.m_x/s); a.m_y = (TYPE1)(a.m_y/s); a.m_z = (TYPE1)(a.m_z/s);
	return a;
}

/* arithematics */

/* negative : -v */
template<class TYPE>
inline xVec3<TYPE> operator-(const xVec3<TYPE>& v)
{
	return xVec3<TYPE>(-v.m_x, -v.m_y, -v.m_z);
}

/* addition : a+b */
template<class TYPE>
inline xVec3<TYPE> operator + (const xVec3<TYPE>& a, const xVec3<TYPE>& b)
{
	return xVec3<TYPE>(a.m_x+b.m_x, a.m_y+b.m_y, a.m_z+b.m_z);
}

/* subtraction : a-b */
template<class TYPE>
inline xVec3<TYPE> operator - (const xVec3<TYPE>& a, const xVec3<TYPE>& b)
{
	return xVec3<TYPE>(a.m_x-b.m_x, a.m_y-b.m_y, a.m_z-b.m_z);
}

/* scaling : a*s */
template<class TYPE>
inline xVec3<TYPE> operator * (const xVec3<TYPE>& a, TYPE s)
{
	return xVec3<TYPE>(a.m_x*s, a.m_y*s, a.m_z*s);
}

/* scaling : s*a */
template<class TYPE>
inline xVec3<TYPE> operator * (TYPE s, const xVec3<TYPE>& a)
{
	return xVec3<TYPE>(a.m_x*s, a.m_y*s, a.m_z*s);
}

/* division : a/s */
template<class TYPE>
inline xVec3<TYPE> operator / (const xVec3<TYPE>& a, TYPE s)
{
	return xVec3<TYPE>(a.m_x/s, a.m_y/s, a.m_z/s);
}

/* cross product : a^b */
template<class TYPE>
inline xVec3<TYPE> operator ^ (const xVec3<TYPE>& a, const xVec3<TYPE>& b)
{
	return xVec3<TYPE>(a.m_y*b.m_z-a.m_z*b.m_y, a.m_z*b.m_x-a.m_x*b.m_z, a.m_x*b.m_y-a.m_y*b.m_x);
}

/* dot/inner product : a*b */
template<class TYPE>
inline TYPE operator * (const xVec3<TYPE>& a, const xVec3<TYPE>& b)
{
	return (a.m_x*b.m_x+a.m_y*b.m_y+a.m_z*b.m_z);
}

/* equality comparision : a==b ? */
template<class TYPE>
inline bool operator == (const xVec3<TYPE>& a, const xVec3<TYPE>& b)
{
	return (a.m_x==b.m_x&& a.m_y==b.m_y&& a.m_z==b.m_z);
}

/* inequality comparision : a!=b ? */
template<class TYPE>
inline bool operator != (const xVec3<TYPE>& a, const xVec3<TYPE>& b)
{
	return (a.m_x!=b.m_x || a.m_y!=b.m_y || a.m_z!=b.m_z);
}

typedef xVec3<double> xVec3d;
typedef xVec3<float>  xVec3f;
typedef xVec3<int>  xVec3i;
/* arithmetic between different kind of 3D vectors
/* addition : a+b */
/* new
inline xVec3d	operator + (const xVec3d& a, const xVec3f& b)
{
	return xVec3d(a.m_x+b.m_x, a.m_y+b.m_y, a.m_z+b.m_z);
}

inline xVec3d	operator + (const xVec3f& a, const xVec3d& b)
{
	return xVec3d(a.m_x+b.m_x, a.m_y+b.m_y, a.m_z+b.m_z);
}

// subtraction : a-b 
// new
inline xVec3d	operator - (const xVec3d& a, const xVec3f& b)
{
	return xVec3d(a.m_x-b.m_x, a.m_y-b.m_y, a.m_z-b.m_z);
}

inline xVec3d	operator - (const xVec3f& a, const xVec3d& b)
{
	return xVec3d(a.m_x-b.m_x, a.m_y-b.m_y, a.m_z-b.m_z);
}
 
// scaling : a*s 
inline xVec3d	operator * (const xVec3f& a, double s)
{
	return xVec3d(a.m_x*s, a.m_y*s, a.m_z*s);
}

inline xVec3d	operator * (const xVec3d& a, float s)
{
	return xVec3d(a.m_x*s, a.m_y*s, a.m_z*s);
}

// scaling : s*a 
inline xVec3d	operator * (double s, const xVec3f& a)
{
	return xVec3d(a.m_x*s, a.m_y*s, a.m_z*s);
}

inline xVec3d	operator * (float s, const xVec3d& a)
{
	return xVec3d(a.m_x*s, a.m_y*s, a.m_z*s);
}

// division : a/s 
inline xVec3d	operator / (const xVec3f& a, double s)
{
	return xVec3d(a.m_x/s, a.m_y/s, a.m_z/s);
}

inline xVec3d	operator / (const xVec3d& a, float s)
{
	return xVec3d(a.m_x/s, a.m_y/s, a.m_z/s);
}

// cross product : a^b 
inline xVec3d operator ^ (const xVec3d& a, const xVec3f& b)
{
	return xVec3d(a.m_y*b.m_z-a.m_z*b.m_y, a.m_z*b.m_x-a.m_x*b.m_z, a.m_x*b.m_y-a.m_y*b.m_x);
}

inline xVec3d operator ^ (const xVec3f& a, const xVec3d& b)
{
	return xVec3d(a.m_y*b.m_z-a.m_z*b.m_y, a.m_z*b.m_x-a.m_x*b.m_z, a.m_x*b.m_y-a.m_y*b.m_x);
}

// dot/inner product : a*b 
inline double operator * (const xVec3d& a, const xVec3f& b)
{
	return (a.m_x*b.m_x+a.m_y*b.m_y+a.m_z*b.m_z);
}

inline double operator * (const xVec3f& a, const xVec3d& b)
{
	return (a.m_x*b.m_x+a.m_y*b.m_y+a.m_z*b.m_z);
}

*/
#endif // #ifndef __MESHOPTI_VEC3_H__

