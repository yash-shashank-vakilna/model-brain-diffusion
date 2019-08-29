#ifndef QUATERNION_H
#define QUATERNION_H

/*
 *
 * quaternion.h : implements quaternions in 3D
 * BRIAN Software Package Version 3.0
 *
 * $Id: quaternion.h 407 2016-10-02 21:36:46Z frithjof $
 *
 * 0.10 (10/04/16): initial version
 * v406 (28/09/16): bumped to version 3.0
 *
 * adapted from Irrlicht engine
 * Copyright (C) 2002-2012 Nikolaus Gebhardt
 *
 */

/*! \file
    \brief Implements quaternions in 3D.
*/

//! Represents quaternions in 3D.

template <typename T> class quaternion {
	T	x;				//!< vectorial (imaginary) part
	T	y;
	T	z;
	T	w;				//!< real part

	quaternion& set(const T ex, const T ey, const T ez)				//! sets new quaternion based on Euler angles.
		{ const T sr = std::sin(ex*0.5), cr = std::cos(ex*0.5);
		  const T sp = std::sin(ey*0.5), cp = std::cos(ey*0.5);
		  const T sy = std::sin(ez*0.5), cy = std::cos(ez*0.5);
		  const T cpcy = cp*cy, spcy = sp*cy, cpsy = cp*sy, spsy = sp*sy;
		  x = sr*cpcy-cr*spsy; y = cr*spcy+sr*cpsy;
		  z = cr*cpsy-sr*spcy; w = cr*cpcy+sr*spsy; return normalize(); }
public:
	quaternion()									//! constructs empty quaternion.
		 : x(0.0), y(0.0), z(0.0), w(1.0) {}
	quaternion(const T _x, const T _y, const T _z, const T _w)			//! constructs from elements.
		 : x(_x), y(_y), z(_z), w(_w) { }
	quaternion(const T ex, const T ey, const T ez)					//! constructs quaternion from Euler angles (in radians).
		 : x(0), y(0), z(0), w(0) { set(ex,ey,ez); }
	quaternion(const vec3<T>& v)							//! constructs quaternion from Euler angles (in radians).
		{ set(v.x,v.y,v.z); }
	quaternion(const mat3<T>& mat)							//! constructs quaternion from matrix.
		{ *this = mat; }
	quaternion& id()								//! sets quaternion to identity.
		{ x = 0.0; y = 0.0; z = 0.0; w = 1.0; return *this; }
	bool	operator==(const quaternion& b) const					//! returns true if quaternions are equal.
		{ return x == b.x && y == b.y && z == b.z && w == b.w; }
	bool	operator!=(const quaternion& b) const					//! returns true if quaternions are not equal.
		{ return !(*this == b); }
	quaternion& operator=(const quaternion& b) = default;				//! assigns from quaternion b.
	quaternion& operator=(const mat3<T>& m)						//! assigns from matrix m.
		{ const T d = m[0]+m[5]+m[10]+1.0;
		  if (d > 0.0) { const T s = std::sqrt(d)*2.0;
			x = (m[6]-m[9])/s; y = (m[8]-m[2])/s; z = (m[1]-m[4])/s; w = 0.25*s; }
		  else if (m[0] > m[5] && m[0] > m[10])	{
			const T s = std::sqrt(1.0+m[0]-m[5]-m[10])*2.0;
			x = 0.25*s; y = (m[4]+m[1])/s; z = (m[2]+m[8])/s; w = (m[6]-m[9])/s; }
		  else if (m[5]>m[10]) {
			const T s = std::sqrt(1.0+m[5]-m[0]-m[10])*2.0;
			x = (m[4]+m[1])/s; y = 0.25*s; z = (m[9]+m[6])/s; w = (m[8]-m[2])/s; }
		  else { const T s = std::sqrt(1.0+m[10]-m[0]-m[5])*2.0;
			x = (m[8]+m[2])/s; y = (m[9]+m[6])/s; z = 0.25*s; w = (m[1]-m[4])/s; }
		  return normalize(); }
	quaternion operator+(const quaternion& b) const					//! adds two quaternions.
		{ return quaternion(x+b.x, y+b.y, z+b.z, w+b.w); }
	quaternion operator*(const quaternion& b) const					//! multiplies two quaternions.
		{ return quaternion(b.w*x+b.x*w+b.y*z-b.z*y, b.w*y+b.y*w+b.z*x-b.x*z,
			b.w*z+b.z*w+b.x*y-b.y*x, b.w*w-b.x*x-b.y*y-b.z*z); }
	vec3<T>	operator*(const vec3<T>& v) const					//! applies quaternion to vector v.
		{ vec3<T> r(x,y,z), q = cross(r,v)+v*w; return v+cross(r*2.0,q); }
	quaternion operator*(const T s) const						//! multiplies quaternion by scalar s.
		{ return quaternion(s*x,s*y,s*z,s*w); }
	quaternion& operator*=(const T s)						//! multiplies this quaternion by scalar s.
		{ x *= s; y *= s; z *= s; w *= s; return *this; }
	quaternion& operator*=(const quaternion& b)					//! multiplies this quaternion by quaternion b.
		{ return *this = b * *this; }
	T	dot(const quaternion& b) const						//! calculates the dot product.
		{ return x*b.x+y*b.y+z*b.z+w*b.w; }
	quaternion& normalize()								//! normalizes this quaternion.
		{ const T n = SQR(x)+SQR(y)+SQR(z)+SQR(w);
		  if (n != 1.0) *this *= 1.0/std::sqrt(n);
		  return *this; }
	mat3<T>	toMatrix(const vec3<T>& c = vec3<T>(0.0)) const				//! converts quaternion to matrix.
		{ mat3<T> m;
		  m[0] = 1.0-2.0*(y*y+z*z); m[1] = 2.0*(x*y+z*w); m[2] = 2.0*(x*z-y*w);
		  m[4] = 2.0*(x*y-z*w); m[5] = 1.0-2.0*(x*x+z*z); m[6] = 2.0*(z*y+x*w);
		  m[8] = 2.0*(x*z+y*w); m[9] = 2.0*(z*y-x*w); m[10] = 1.0-2.0*(x*x+y*y);
		  m[12] = c.x; m[13] = c.y; m[14] = c.z; m[15] = 1.0; return m; }
	matD<T>	rotMatrix() const							//! converts quaternion to 3x3 rotation matrix.
		{ matD<T> m(3,3);
		  m(0,0) = 1.0-2.0*(y*y+z*z); m(0,1) = 2.0*(x*y+z*w); m(0,2) = 2.0*(x*z-y*w);
		  m(1,0) = 2.0*(x*y-z*w); m(1,1) = 1.0-2.0*(x*x+z*z); m(1,2) = 2.0*(z*y+x*w);
		  m(2,0) = 2.0*(x*z+y*w); m(2,1) = 2.0*(z*y-x*w); m(2,2) = 1.0-2.0*(x*x+y*y);
		  return m; }
	quaternion& inv()								//! inverts this quaternion.
		{ x = -x; y = -y; z = -z; return *this; }
	quaternion& fromAngleAxis(const T ang, const vec3<T>& ax)			//! constructs quaternion from angle (in radians) and axis.
		{ const T s = std::sin(0.5*ang); w = std::cos(0.5*ang); 
		  x = s*ax.x; y = s*ax.y; z = s*ax.z; return *this; }
	void	toAngleAxis(T& ang, vec3<T>& ax) const					//! returns angle (in radians) and axis from this quaternion.
		{ const T s = std::sqrt(SQR(x)+SQR(y)+SQR(z));
		  if (isSmall<T>(s) || w > 1.0 || w < -1.0) { ang = 0.0; ax = 0.0; }
		  else { ang = 2.0*std::acos(w); ax.x = x/s; ax.y = y/s; ax.z = z/s; } }
	vec3<T>	toEuler() const								//! converts quaternion to Euler angle (in radians)
		{ const T t = 2.0*(y*w-x*z);
		  if (std::abs(t-1.0) < 0.000001)
			return vec3<T>(0.0,M_PI/2.0,-2.0*std::atan2(x,w));
		  else if (std::abs(t+1.0) < 0.000001)
			return vec3<T>(0.0,M_PI/-2.0,2.0*std::atan2(x,w));
	 	  else { return vec3<T>(std::atan2(2.0*(y*z+x*w),-x*x-y*y+z*z+w*w),		  
				std::asin(CLAMP(t,-1.0,1.0)),
				std::atan2(2.0*(x*y+z*w),x*x-y*y-z*z+w*w)); } }
	quaternion& rotationFromTo(const vec3<T>& f, const vec3<T>& t)			//! set quaternion to represent a rotation from vector f to t.
		{ vec3<T> v0 = f.normalize(), v1 = t.normalize(), c;
		  T d = ::dot(v0,v1), s = 0.0;
		  if (d >= 1.0) return id();						// v0,v1 are the same
		  else if (d <= -1.0) { c = cross(vec3<T>(1.0,0.0,0.0),v0);		// exactly opposite
			if (norm(c) == 0.0) c = cross(vec3<T>(0.0,1.0,0.0),v0); }
		  else { s = std::sqrt((1.0+d)*2.0); c = cross(v0,v1)/s; }
		  x = c.x; y = c.y; z = c.z; w = s*0.5; return normalize(); }
	void	print()	const								//! prints quaternion b to stdout.
		{ printT(x); printT(y); printT(z); printT(w); printf("\n"); }
};

using fquat = quaternion<float>;
using dquat = quaternion<double>;
using vfquat = std::vector<fquat>;
using vdquat = std::vector<dquat>;

#endif


