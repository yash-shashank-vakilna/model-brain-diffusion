#ifndef GEOMTESTS_H
#define GEOMTESTS_H

/*
 *
 * geomTests.h: (inexact) geometrical tests
 * BRIAN Software Package Version 3.0
 *
 * $Id: geomTests.h 410 2016-10-06 00:48:28Z frithjof $
 *
 * 0.10 (17/11/09): adapted for BRIAN 2.0
 * 0.20 (02/12/12): checkXX usage in intersection tests corrected
 * 0.30 (31/12/12): released version 2.4
 * 0.40 (16/12/13): documented
 * v308 (26/03/16): numerous bug fixes and improvements
 * v406 (28/09/16): bumped to version 3.0
 *
 */

/*! \file
    \brief A collection of geometric tests used in mesh algorithms.
 */

template<typename T>
vec3<T> toCart(const T phi, const T theta)
//! maps from unit spherical to Cartesian coordinates.
{	return vec3<T>(std::cos(theta)*std::sin(phi),
		std::sin(theta)*std::sin(phi), std::cos(phi)); }

template<typename T>
void toPolar(const vec3<T>& v, T& phi, T& theta)
//! maps from Cartesian to unit spherical coordinates.
{	phi = std::acos(v.z); theta = std::atan2(v.y, v.x); }

template<typename T>
T cos(const vec3<T>& u, const vec3<T>& v)
//! returns the cosine of the angle between vectors (u,v).
{	return CLAMP(dot(u,v),T(-1),T(1)); }

template<typename T>
T cos(const vec3<T>& v0, const vec3<T>& v1, const vec3<T>& v2)
//! returns cos of angle at v0.
{	return cos(v1-v0, v2-v0); }

template<typename T>
bool obtuseAngle(const vec3<T>& v0, const vec3<T>& v1, const vec3<T>& v2)
//! returns true if angle at v0 is obtuse.
{	return cos(v0, v1, v2) < T(0); }

template<typename T>
T angle(const vec3<T>& u, const vec3<T>& v)
//! returns angle between vectors (u,v).
{	const T uu = dot(u,u), uv = dot(u,v), vv = dot(v,v), den = uu*vv-uv*uv;
	assertNZ(den >= T(0)); return atan2(std::sqrt(den), uv); }

template<typename T>
T cotan(const vec3<T>& u, const vec3<T>& v)
//! return cotan of angle between vectors (u,v).
{	const T uu = dot(u,u), uv = dot(u,v), vv = dot(v,v), den = uu*vv-uv*uv;
	return checkNZ(den)? uv/std::sqrt(den): std::numeric_limits<T>::max(); }

template<typename T>
T angle(const vec3<T>& a, const vec3<T>& b, const vec3<T>& c)
//! returns angle at a. 
{ return angle(b-a,c-a); }

template<typename T>
vec3<T> rawNormal(const vec3<T>& a, const vec3<T>& b, const vec3<T>& c)
//! returns raw normal.
{ return cross(b-a,c-a); }

template<typename T>
bool normal(const vec3<T>& a, const vec3<T>& b, const vec3<T>& c, vec3<T>& n)
//! returns the plane normal n of triangle (a,b,c).
{	n = rawNormal(a,b,c); const T l = dot(n,n); n = n/T(std::sqrt(l));
	return checkNZ(l); }

template<typename T>
T area(const vec3<T>& a, const vec3<T>& b, const vec3<T>& c)
//! returns area of triangle (a,b,c).
{	const vec3<T> n = rawNormal(a,b,c); return T(0.5)*norm(n); }

template<typename T>
T pointDist(const vec3<T>& a, const vec3<T>& b)
//! returns distance between points (a,b).
{	return a == b? T(0): norm(b-a); }

template<typename T>
bool pointPointCollision(const vec3<T>& a, const vec3<T>& b)
//! returns true if distance between points (a,b) is less than epsilon.
{	return checkEPS<T>(pointDist(a,b)) == false; }

template<typename T>
T geoDist(const vec3<T>& a, const vec3<T>& b)
//! returns geodesic distance between spherical points (a,b).
{	if (a == b) return T(0); const T d = dot(a,b);  
	return std::acos(CLAMP(d,T(-1),T(1))); }					// safeguard against rounding errors

template<typename T>
T planeDist(const vec3<T>& a, const vec3<T>& n, const vec3<T>& p)
//! compute distance of point p to plane (a,n).
{	return a == p? T(0): dot(p-a,n); }

template<typename T>
vec3<T> projectPoint(const vec3<T>& a, const vec3<T>& n, const vec3<T>& p)
//! projects point p onto plane (a,n).
{	return p-planeDist(a,n,p)*n; };

template<typename T>
vec3<T> projectPoint(const vec3<T>& a, const vec3<T>& b, const vec3<T>& c, const vec3<T>& p)
//! returns point p projected onto triangle plane (a,b,c).
{	vec3<T> n; normal(a,b,c,n); return projectPoint(a,n,p); }

template<typename T>
vec3<T> baryWeights(const vec3<T>& a, const vec3<T>& b, const vec3<T>& c, const vec3<T>& p)
//! returns barycentric interpolation weights for point p in triangle (a,b,c).
{	const vec3<T> q(projectPoint(a,b,c,p)), w(area(q,b,c),area(q,a,c),area(q,a,b));
	return w/w.sum(); };

template<typename T>
T interpolateAt(const vec3<T>& a, const vec3<T>& b, const vec3<T>& c, const vec3<T>& p, 
	const T d0, const T d1, const T d2)
//! returns interpolated scalar (d0,d1,d2) in triangle (a,b,c) at p.
{	const vec3<T> w = baryWeights(a,b,c,p); return w.x*d0+w.y*d1+w.z*d2; }    
    
template<typename T>
vec3<T> interpolateAt(const vec3<T>& a, const vec3<T>& b, const vec3<T>& c,
	const vec3<T>& p, const vec3<T>& d0, const vec3<T>& d1, const vec3<T>& d2)
//! returns interpolated vector (d0,d1,d2) in triangle (a,b,c) at p.
{	const vec3<T> w = baryWeights(a,b,c,p); return w.x*d0+w.y*d1+w.z*d2; }    
    
template<typename T>
T bary(const vec3<T>& a, const vec3<T>& b, const vec3<T>& c, const vec3<T>& p)
//! returns barycenter coordinate in triangle (a,b,c) at p.
{	const vec3<T> q(projectPoint(a,b,c,p)); const T t = area(a,b,c);
	return checkEPS(t)? area(q,b,c)/t: T(1); }

template<typename T>
T baryBound(const vec3<T>& a, const vec3<T>& b, const vec3<T>& c, const vec3<T>& p)
//! returns bounded barycenter coordinate of point p in triangle (a,b,c).
{	return std::min(bary(a,b,c,p),T(1)); }

template<typename T>
T VoronoiArea(const vec3<T>& p, const vec3<T>& q, const vec3<T>& r)
//! computes mixed/Voronoi area of triangle (p,q,r).
{	const vec3<T> a = r-q, b = p-r, c = q-p;
	const T al = dot(b,c), bt = dot(a,c), gm = dot(b,a);
	if (al < T(0) || bt < T(0) || gm < T(0)) {
		const T ar = norm(cross(b,c));
		return al < T(0)? T(0.5)*ar: T(0.25)*ar; }
	else return T(0.25)*(dot(b,b)*cotan(c,a)+dot(c,c)*cotan(b,a));
}

template<typename T>
bool pointEdgeCollision(const vec3<T>& a, const vec3<T>& b, const vec3<T>& p)
//! checks if point p lies on the edge (a,b).
{	const vec3<T> pa = p-a, bp = b-p, ba = b-a, cp = cross(pa,ba);
	const T dpa = norm(pa), dbp = norm(bp), dba = norm(ba), dpl = norm(cp);
	return checkEPS<T>(dba-dpa-dbp)? false: checkEPS<T>(dpl/dba) == false;
}

template<typename T>
T edgeParameter(const vec3<T>& a, const vec3<T>& b, const vec3<T>& p)
//! returns fractional distance of p which is on the edge (a,b)
{	const T dba = pointDist(a,b), dpa = pointDist(a,p);
	const T t = dba > T(0)? dpa/dba: T(1); return std::min(t,T(1));
}

template<typename T>
bool pointInTriangle(const vec3<T>& a, const vec3<T>& b, const vec3<T>& c,
	const vec3<T>& p)
//! checks if point p is inside triangle (a,b,c).
{	const vec3<T> u = b-a, v = c-a, w = p-a;
	const T uv = dot(u,v), wv = dot(w,v), vv = dot(v,v), uu = dot(u,u);
	const T wu = dot(w,u), den = uv*uv-uu*vv; assertNZ(den);
	const T s = (uv*wv-vv*wu)/den, t = (uv*wu-uu*wv)/den;
	return (s >= T(0) && s <= T(1)) && (t >= T(0) && s+t <= T(1));
}

template<typename T>
bool edgeLineIntersection(const vec3<T>& a, const vec3<T>& n0, const vec3<T>& p,
	const vec3<T>& q, vec3<T>& i)
//! checks if line (a,n0) and edge (p,q) intersect and returns intersection point in i.
{	const vec3<T> n1 = q-p, n10 = cross(n0,n1), p10 = p-a;
	const T den = dot(n10,n10); if (checkNZ(den) == false) return false;		// lines are parallel
	const T s = det3(p10,n1,n10)/den, t = det3(p10,n0,n10)/den;
	if (t < T(0) || t > T(1)) return false;						// no intersection
	const vec3<T> p11 = a+s*n0; i = p+t*n1;
	return checkEPS<T>(pointDist(p11, i)) == false;
}

template<typename T>
bool edgeEdgeIntersection(const vec3<T>& a, const vec3<T>& b, const vec3<T>& p,
	const vec3<T>& q, vec3<T>& i)
//! checks if edges (a,b) and (p,q) intersect and returns intersection point in i.
{	const vec3<T> n1 = q-p, n0 = b-a, n10 = cross(n0,n1), p10 = p-a;
	const T den = dot(n10,n10); if (checkNZ(den) == false) return false;		// lines are parallel
	const T s = det3(p10,n1,n10)/den, t = det3(p10,n0,n10)/den;
	if (s < T(0) || s > T(1) || t < T(0) || t > T(1)) return false;			// no intersection
	const vec3<T> p11 = a+s*n0; i = p+t*n1;
	return checkEPS<T>(pointDist(p11, i)) == false;
}

template<typename T>
bool triangleEdgeIntersection(const vec3<T>& a, const vec3<T>& b, const vec3<T>& c,
	const vec3<T>& p, const vec3<T>& q, vec3<T>& i)
//! checks if triangle (a,b,c) intersects with edge (p,q) and returns intersection point i.
{	vec3<T> n; i = p;
	if (pointPointCollision(a,p) || pointPointCollision(b,p) || 
		pointPointCollision(c,p)) return true;
	if (pointPointCollision(a,q) || pointPointCollision(b,q) ||
		pointPointCollision(c,q)) return true;
	if (edgeEdgeIntersection(p,q,a,b,i) || edgeEdgeIntersection(p,q,b,c,i) ||
		edgeEdgeIntersection(p,q,c,a,i)) return true;
	if (normal(a,b,c,n)) {
	const T dp = planeDist(a,n,p), dq = planeDist(a,n,q), eps = safeMin<T>();
	if ((dp > eps && dq > eps) || (dp < -eps && dq < -eps)) return false;
	else if (std::abs(dp) > eps && std::abs(dq) > eps) {
		const T t = dp/(dp-dq); assert(t >= T(0) && t <= T(1)); i = p+t*(q-p); }
	else if (std::abs(dp) < eps && std::abs(dq) > eps) i = p-dp*n;
	else if (std::abs(dp) > eps && std::abs(dq) < eps) i = q-dq*n;
	return pointInTriangle(a,b,c,i); };
	return false;
}

template<typename T>
bool edgeTriangleIntersection(const vec3<T>& a, const vec3<T>& b, const vec3<T>& c,
	const vec3<T>& p, const vec3<T>& q, vec3<T>& i)
//! checks if edge (p,q) intersects with triangle (a,b,c) and returns intersection point i.
{	vec3<T> n; i = p;
	if (pointPointCollision(a,p) || pointPointCollision(b,p) || 
		pointPointCollision(c,p)) return true;
	if (pointPointCollision(a,q) || pointPointCollision(b,q) || 
		pointPointCollision(c,q)) return true;
	if (edgeEdgeIntersection(a,b,p,q,i) || edgeEdgeIntersection(b,c,p,q,i) || 
		edgeEdgeIntersection(c,a,p,q,i)) return true;
	if (normal(a,b,c,n))  {
		const T dp = planeDist(a,n,p), dq = planeDist(a,n,q), eps = safeMin<T>();
		if ((dp > eps && dq > eps) || (dp < -eps && dq < -eps)) return false;
		else if (std::abs(dp) > eps && std::abs(dq) > eps) {
			const T t = dp/(dp-dq); assert(t >= T(0) && t <= T(1)); i = p+t*(q-p); }
		else if (std::abs(dp) < eps && std::abs(dq) > eps) i = p-dp*n;
		else if (std::abs(dp) > eps && std::abs(dq) < eps) i = q-dq*n;
		return pointInTriangle(a,b,c,i); };
	return false;
}

template<typename T>
bool pointInTetra(const vec3<T>& a, const vec3<T>& b, const vec3<T>& c,
	const vec3<T>& d, const vec3<T>& p)
//! checks if point p is inside tetrahedron (a,b,c,d).
{	vec3<T> n, q; T dd, dp, eps = safeMin<T>(); int i = 0;
	bool s1 = false, s2 = false, s3 = false, s4 = false;
	if (pointPointCollision(a,p) || pointPointCollision(b,p) || 
		pointPointCollision(c,p) || pointPointCollision(d,p)) return true;
	if (pointEdgeCollision(a,b,p) || pointEdgeCollision(a,c,p) || 
		pointEdgeCollision(a,d,p)) return true;
	if (pointEdgeCollision(b,c,p) || pointEdgeCollision(b,d,p) || 
		pointEdgeCollision(c,d,p)) return true;
	if (normal(a,b,c,n)) { dd = planeDist(a,n,d); dp = planeDist(a,n,p);
		s1 = (dd > T(0) && dp > T(0)) || (dd < T(0) && dp < T(0));
		if (s1 == false && std::abs(dp) < eps) { q = p-dp*n;
			if (pointInTriangle(a,b,c,q)) return true; } } else i++;
	if (normal(b,c,d,n)) { dd = planeDist(b,n,a); dp = planeDist(b,n,p);
		s2 = (dd > T(0) && dp > T(0)) || (dd < T(0) && dp < T(0));
		if (s2 == false && std::abs(dp) < eps) { q = p-dp*n;
			if (pointInTriangle(b,c,d,q)) return true; } } else i++;
	if (normal(c,d,a,n)) { dp = planeDist(c,n,p); dd = planeDist(c,n,b);
		s3 = (dd > T(0) && dp > T(0)) || (dd < T(0) && dp < T(0));
		if (s3 == false && std::abs(dp) < eps) { q = p-dp*n;
			if (pointInTriangle(c,d,a,q)) return true; } } else i++;
	if (normal(d,a,b,n)) { dp = planeDist(d,n,p); dd = planeDist(d,n,c);
		s4 = (dd > T(0) && dp > T(0)) || (dd < T(0) && dp < T(0));
		if (s4 == false && std::abs(dp) < eps) { q = p-dp*n; 
			if (pointInTriangle(d,a,b,q)) return true; } } else i++;
	if (i == 4) { q = b; const T db = pointDist(a, b), dc = pointDist(a, c);
		if (dc > db) { q = c; }; dd = pointDist(a, d);
		if (dc > dc && dd > db) { q = d; }; return pointEdgeCollision(a, q, p); };
	return s1 && s2 && s3 && s4;
}
#endif

