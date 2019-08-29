#ifndef SIMPMESH_H
#define SIMPMESH_H

/*
 *
 * simplifyMesh.h: simplify triangular meshes
 * BRIAN Software Package Version 3.0
 *
 * $Id: simplifyMesh.h 436 2016-11-06 00:42:43Z frithjof $
 *
 * 0.10 (19/11/09): initial version
 * 0.30 (31/12/12): released version 2.4
 * 0.40 (19/12/13): documented
 *
 */

/*! \file
    \brief Templated classes for simplifying triangular meshes.
*/

// forward definitions
template<typename T> struct Vertex;
template<typename T> struct Edge;
template<typename T> struct Triangle;

//! Helper class that provides a link to a vertex

template<typename T>
struct LinkVertex {
	Vertex<T>*	vertex;				//!< the vertex
	unsigned int	number;				//!< vertex id
	bool		isElement;			//!< if vertex exists
	
	LinkVertex(Vertex<T>* const v)							//! allocates a record for vertex v.
		 : vertex(v), number(1), isElement(true) { }
};

//! Helper class that provides a link to an edge

template<typename T>
struct LinkEdge {
	Vertex<T>*	vertices[2];			//!< the vertices of this edge
	bool		isElement;			//!< if edge exists
	
	LinkEdge(Vertex<T>* const v0, Vertex<T>* const v1)				//! allocates a record for vertices (v0,v1).
	 : isElement(true) { assert(v0 && v1); vertices[0] = v0; vertices[1] = v1; }
};

//! Implements a set of links to a vertices and edges

template<typename T>
struct Link {
	std::deque<LinkVertex<T>> vertices;		//!< set of vertex links
	std::deque<LinkEdge<T>> edges;			//!< set of edge links
	Link(Vertex<T>* const v, Vertex<T>* const dummy = 0);
	bool	isHalfDisk(const unsigned int sedge, const unsigned int svertex, const Vertex<T>* const stopvertex);
	void	intersect(Vertex<T>* const v, Vertex<T>* const dummy = 0);
	void	except(Edge<T>* const e, Vertex<T>* const dummy = 0);
	bool	emptyEdges() const;
	bool	empty() const;
};

template<typename T>
Link<T>::Link(Vertex<T>* const v, Vertex<T>* const dummy)
//! allocates a link between vertex v and vertex dummy.
 : vertices(), edges()
{	for (const auto vi : v->triangles) {
		Vertex<T>* tmp[2]; unsigned int j = 0;
		if (vi->vertices[0] != v) tmp[j++] = vi->vertices[0];
		if (vi->vertices[1] != v) tmp[j++] = vi->vertices[1];
		if (vi->vertices[2] != v) tmp[j++] = vi->vertices[2];
		assert(j == 2); j = 0;
		while (j < vertices.size() && tmp[0] != vertices[j].vertex) ++j;
		if (j == vertices.size()) vertices.push_back(LinkVertex<T>(tmp[0])); else ++vertices[j].number;
		j = 0;
		while (j < vertices.size() && tmp[1] != vertices[j].vertex) ++j;
		if (j == vertices.size()) vertices.push_back(LinkVertex<T>(tmp[1])); else ++vertices[j].number;
		edges.push_back(LinkEdge<T>(tmp[0], tmp[1])); }
	if (dummy) {
		vertices.push_back(LinkVertex<T>(dummy));
		for (const auto e : v->edges) {
			if (e->getDegree() != 2) {
				if (e->vertices[0] != v) edges.push_back(LinkEdge<T>(dummy, e->vertices[0]));
				else edges.push_back(LinkEdge<T>(dummy, e->vertices[1])); } }	}
}

template<typename T>
bool Link<T>::isHalfDisk(const unsigned int startedge, const unsigned int startvertex, const Vertex<T>* const stopvertex)
//! checks the configuration of this link.
{	unsigned int curedge = startedge; unsigned int curvertex = startvertex;
	curvertex = (curvertex+1)%2;
	if (!edges[curedge].isElement) return false;
	edges[curedge].isElement = false;
	while (edges[curedge].vertices[curvertex] != stopvertex) {
		unsigned int i = 0; if (i == curedge) ++i;
		while (edges[curedge].vertices[curvertex] != edges[i].vertices[0] &&
		       edges[curedge].vertices[curvertex] != edges[i].vertices[1]) {
			++i; if(i == curedge) ++i; }
		curvertex = (edges[curedge].vertices[curvertex] == edges[i].vertices[0])? 0: 1;
		curedge = i; curvertex = (curvertex+1)%2;
		if (!edges[curedge].isElement) return false;
		edges[curedge].isElement = false; }
	return true;
}

template<typename T>
void Link<T>::intersect(Vertex<T>* const v0, Vertex<T>* const dummy)
//! used for assessing topology preservation.
{	for (auto& v : vertices) {
		if (v.isElement) { v.isElement = false;
			if (dummy && dummy == v.vertex) v.isElement = true;
			else if (v.vertex != v0) {
				for (const auto t : v.vertex->triangles) {
					if (t->vertices[0] == v0) v.isElement = true;
					else if (t->vertices[1] == v0) v.isElement = true;
					else if (t->vertices[2] == v0) v.isElement = true; } } } };
	for (auto& e : edges) { if (e.isElement == false) continue;
		e.isElement = false;
		if (e.vertices[0] != v0 && e.vertices[1] != v0) {
			if (dummy && dummy == e.vertices[0]) {
				for (const auto e0 : v0->edges) {
					if (e0->vertices[0] == e.vertices[1] && e0->getDegree() != 2) e.isElement = true;
					else if (e0->vertices[1] == e.vertices[1] && e0->getDegree() != 2) e.isElement = true; } }
			else if (dummy && dummy == e.vertices[1]) {
				for (const auto e0 : v0->edges) {
					if (e0->vertices[0] == e.vertices[0] && e0->getDegree() != 2) e.isElement = true;
					else if (e0->vertices[1] == e.vertices[0] && e0->getDegree() != 2) e.isElement = true; } }
			else {	for (const auto t0 : e.vertices[0]->triangles) {
					bool tmp = false;
					if (t0->vertices[0] == v0) tmp = true;
					else if (t0->vertices[1] == v0) tmp = true;
					else if (t0->vertices[2] == v0) tmp = true;
					if (tmp == true) {
						if (t0->vertices[0] == e.vertices[1]) e.isElement = true;
						else if (t0->vertices[1] == e.vertices[1]) e.isElement = true;
						else if (t0->vertices[2] == e.vertices[1]) e.isElement = true; } } } } };
}

template<typename T>
void Link<T>::except(Edge<T>* const e, Vertex<T>* const dummy)
//! used for assessing topology preservation.
{	for (const auto t0 : e->vertices[0]->triangles) {
		for (const auto t1 : e->vertices[1]->triangles) { if (t0 != t1) continue;
			Vertex<T>* tmp = 0;
			if (t0->vertices[0] != e->vertices[0] && t0->vertices[0] != e->vertices[1]) tmp = t0->vertices[0];
			else if (t0->vertices[1] != e->vertices[0] && t0->vertices[1] != e->vertices[1]) tmp = t0->vertices[1];
			else if (t0->vertices[2] != e->vertices[0] && t0->vertices[2] != e->vertices[1]) tmp = t0->vertices[2];
			for (auto& v : vertices) if (v.vertex == tmp) v.isElement = !v.isElement; } }
	if (dummy) { for (auto& v : vertices) if (v.vertex == dummy) v.isElement = !v.isElement; }
}

template<typename T>
bool Link<T>::emptyEdges() const
//! checks if edge list has all empty elements.
{	for (const auto& e : edges) if (e.isElement) return false;
	return true;
}

template<typename T>
bool Link<T>::empty() const
//! checks if vertex and edge lists have all empty elements.
{	for (const auto& v : vertices) if (v.isElement) return false;
	for (const auto& e : edges) if (e.isElement) return false;
	return true;
}

//! Implements the quadric error metric for mesh simplification.

template<typename T>
class Quadric {
	T _aa, _ab, _ac, _bb, _bc, _cc, _ad, _bd, _cd, _dd;	//!< parameters of the quadric
	T _area;					//!< triangle area
public:
	Quadric()									//! allocates an empty quadric.
		 :_aa(T(0)), _ab(T(0)), _ac(T(0)), _bb(T(0)), _bc(T(0)), _cc(T(0)),
			_ad(T(0)), _bd(T(0)), _cd(T(0)), _dd(T(0)), _area(T(0)) { }
	Quadric(const vec3<T> &p0, const vec3<T> &p1, const vec3<T> &p2);
	T	operator()(const vec3<T> &p) const;
	Quadric<T>& operator+=(const Quadric<T>& rhs);
	Quadric<T>& operator*=(const T& s);
	Quadric<T>  operator*(const T& s) const;
	T	area() const								//! returns the area of this triangle.
		{ return _area; }
	bool	placeNewVertex1(Edge<T>* const e) const;
	bool	placeNewVertex2(Edge<T>* const e) const;
	void	placeNewVertex3(Edge<T>* const e) const;
};

template<typename T>
Quadric<T>::Quadric(const vec3<T> &p0, const vec3<T> &p1, const vec3<T> &p2) :
	_aa(T(0)), _ab(T(0)), _ac(T(0)), _bb(T(0)), _bc(T(0)), _cc(T(0)),
	_ad(T(0)), _bd(T(0)), _cd(T(0)), _dd(T(0)), _area(T(0))
//! allocates a quadric for triangle (p0,p1,p2).
{	vec3<T> n = cross(p1-p0,p2-p0); T tmp = dot(n,n);
	if (checkEPS(tmp)) { tmp = std::sqrt(tmp); n /= tmp; }
	else { n = 0; tmp = 0; }
	const T d = -dot(n,p0);
	_aa = n.x*n.x; _ab = n.x*n.y; _ac = n.x*n.z; _bb = n.y*n.y; _bc = n.y*n.z; _cc = n.z*n.z;
	_ad = n.x*d; _bd = n.y*d; _cd = n.z*d; _dd = d*d; _area = 0.5f*tmp;
}

template<typename T>
T Quadric<T>::operator()(const vec3<T>& p) const 
//! applies point p to a quadric.
{	return p.x*p.x*_aa+2*p.x*p.y*_ab+2*p.x*p.z*_ac+p.y*p.y*_bb+
		2*p.y*p.z*_bc+p.z*p.z*_cc+2*p.x*_ad+2*p.y*_bd+2*p.z*_cd+_dd;
}

template<typename T>
Quadric<T>& Quadric<T>::operator+=(const Quadric<T>& rhs) 
//! adds a quadric to this one.
{	_aa += rhs._aa; _ab += rhs._ab; _ac += rhs._ac; _bb += rhs._bb;
	_bc += rhs._bc; _cc += rhs._cc;	_ad += rhs._ad; _bd += rhs._bd;
	_cd += rhs._cd; _dd += rhs._dd;	_area += rhs._area; return *this;
}

template<typename T>
Quadric<T>& Quadric<T>::operator*=(const T& s) 
//! multiplies this quadric by a scalar.
{	assert(s >= 0); _aa *= s; _ab *= s; _ac *= s; _bb *= s;
	_bc *= s; _cc *= s; _ad *= s; _bd *= s; _cd *= s; _dd *= s;
	_area *= s; return *this;
}

template<typename T>
Quadric<T> Quadric<T>::operator*(const T& s) const
//! multiplies a quadric and a scalar.
{	assert(s >= 0); Quadric<T> a(*this); a *= s; return a;
}

template<typename T>
bool Quadric<T>::placeNewVertex1(Edge<T>* const e) const 
//! moves an edge and recompute quadric error metric.
{	const T tmp = (_bb*_cc-_bc*_bc)*_aa+(_ac*_bc-_ab*_cc)*_ab+(_ab*_bc-_ac*_bb)*_ac;
	if (std::sqrt(std::numeric_limits<T>::epsilon()) > std::abs(tmp)) return false;
	T invA[9];
	invA[0] = (_bb*_cc-_bc*_bc)/tmp; invA[1] = (_ac*_bc-_ab*_cc)/tmp;
	invA[2] = (_ab*_bc-_ac*_bb)/tmp; invA[3] = invA[1];
	invA[4] = (_aa*_cc-_ac*_ac)/tmp; invA[5] = (_ab*_ac-_aa*_bc)/tmp;
	invA[6] = invA[2]; invA[7] = invA[5]; invA[8] = (_aa*_bb-_ab*_ab)/tmp;
	e->np.x = -(invA[0]*_ad+invA[1]*_bd+invA[2]*_cd);
	e->np.y = -(invA[3]*_ad+invA[4]*_bd+invA[5]*_cd);
	e->np.z = -(invA[6]*_ad+invA[7]*_bd+invA[8]*_cd);
	return true;
}

template<typename T>
bool Quadric<T>::placeNewVertex2(Edge<T>* const e) const 
//! moves an edge and recompute quadric error metric.
{	vec3<T> v0(e->vertices[0]->point); vec3<T> v1(e->vertices[1]->point); vec3<T> Av0;
	Av0.x = _aa*v0.x+_ab*v0.y+_ac*v0.z; Av0.y = _ab*v0.x+_bb*v0.y+_bc*v0.z; Av0.z = _ac*v0.x+_bc*v0.y+_cc*v0.z;
	vec3<T> d = v1-v0; vec3<T> Ad;
	Ad.x = _aa*d.x+_ab*d.y+_ac*d.z; Ad.y = _ab*d.x+_bb*d.y+_bc*d.z; Ad.z = _ac*d.x+_bc*d.y+_cc*d.z;
	const T tmp = dot(d,Ad);
	if (std::sqrt(std::numeric_limits<T>::epsilon()) <= std::abs(tmp)) {
		T t = -(dot(d,Av0)+_ad*d.x+_bd*d.y+_cd*d.z)/tmp;
		if (t < 0) t = 0; else if (t > 1.0) t = T(1.0);
		e->np = v0+t*d; return true;
	} else return false;
}

template<typename T>
void Quadric<T>::placeNewVertex3(Edge<T>* const e) const
//! moves an edge and recompute quadric error metric.
{	T costBegin = operator()(e->vertices[0]->point); T costEnd = operator()(e->vertices[1]->point);
	if (costBegin < costEnd) e->np = e->vertices[0]->point;
	else if (costBegin > costEnd) e->np = e->vertices[1]->point;
	else e->np = 0.5f*(e->vertices[0]->point+e->vertices[1]->point);
}

//! Represents an edge contraction record in mesh simplification.

template<typename T>
struct Edge {
	Vertex<T>* vertices[2];				//!< vertices connected by this edge.
	vec3<T> np;					//!< new position in edge contraction.
	T	cost;					//!< cost of this edge contraction.
	unsigned int heap;				//!< location on heap.
	bool	tp;					//!< the contraction of the edge preserves topology
	bool	cp;					//!< the contraction of the edge is consistent
	
	Edge(Vertex<T>* const v0, Vertex<T>* const v1)					//! allocates a contraction record for edge (v0,v1).
		 : np(), cost(T(0)), heap(0), tp(false), cp(false)
		{ assert(v0 && v1); vertices[0] = v0; vertices[1] = v1; }
	void resize(Vertex<T>* const vn, Vertex<T>* const vo)				//! replaces vertex vo by vn.
		{ const unsigned int i = vertices[0] == vo? 0 : 1; assert(vertices[i] == vo);
		  vertices[i]->remove(this); vertices[i] = vn; vertices[i]->generate(this); }
	unsigned int getDegree() const;
	bool	isTopologyPreserving();
	bool	isConsistent();
	bool	isCoherentOriented() const;
	void	setCost(bool check);
	T	length() const								//! returns edge length.
		{ return norm(vertices[1]->point-vertices[0]->point); }
private:
	bool	isTopologyPreservingAccelerated() const;
	bool	isBoundaryLinkEmpty() const;
};

template<typename T>
unsigned int Edge<T>::getDegree() const
//! returns the degree of this edge.
{	unsigned int degree = 0;
	for (const auto& t : vertices[0]->triangles) {
		if (t->vertices[0] == vertices[1] || t->vertices[1] == vertices[1] || t->vertices[2] == vertices[1]) ++degree; }
	return degree;
}

template<typename T>
bool Edge<T>::isTopologyPreservingAccelerated() const
//! accelerated path for special case 2: deg ab = 2, ord ab = ord a = ord b = 0.
{	std::vector<Vertex<T>* > linkA;
	for (const auto& e : vertices[0]->edges) {
		const unsigned int j = e->vertices[0] != vertices[0]? 0: 1;
		linkA.push_back(e->vertices[j]); }
	unsigned int k = 0; std::vector<Vertex<T>* > linkAB(2);
	for (const auto& e : vertices[1]->edges) {
		unsigned int j = e->vertices[0] != vertices[1]? 0: 1; 
		for (const auto& l : linkA)
			if (e->vertices[j] == l) { if (1 < k) return false; linkAB[k] = l; ++k; break; } }
	assert(k == 2); assert(linkAB[0] != linkAB[1]);
	unsigned int exists = 0;
	for (const auto& t : linkAB[0]->triangles) { unsigned int j;
		if (t->vertices[0] == linkAB[0]) j = 0;
		else if (t->vertices[1] == linkAB[0]) j = 1;
		else j = 2;
		assert(linkAB[0] == t->vertices[j]);
		if (t->vertices[(j+1)%3] == linkAB[1]) {
			if (t->vertices[(j+2)%3] == vertices[0]) {
				if (exists == 0) exists = 1;
				else if (exists == 2) return false;
				else assert(exists != 1); }
			else if (t->vertices[(j+2)%3] == vertices[1]) {
				if (exists == 0) exists = 2;
				else if (exists == 1) return false;
				else assert(exists != 2); } }
		else if (t->vertices[(j+2)%3] == linkAB[1]) {
			if (t->vertices[(j+1)%3] == vertices[0]) {
				if (exists == 0) exists = 1;
				else if (exists == 2) return false;
				else assert(exists != 1); }
			else if (t->vertices[(j+1)%3] == vertices[1]) {
				if (exists == 0) exists = 2;
				else if (exists == 1) return false;
				else assert(exists != 2); } } }
	return true;
}

template<typename T>
bool Edge<T>::isBoundaryLinkEmpty() const 
//! check for special case 3.
{	std::deque<Vertex<T>* > linkA;
	for (const auto& e : vertices[0]->edges) {
		if (e->getDegree() == 2) continue;
		const unsigned int j = e->vertices[0] != vertices[0]? 0: 1;
		linkA.push_back(e->vertices[j]); };
	for (const auto& e : vertices[1]->edges) {
		if (e->getDegree() == 2) continue;
		const unsigned int j = e->vertices[0] != vertices[1]? 0: 1;
		const auto it = std::find(linkA.begin(), linkA.end(), e->vertices[j]);
		if (it != linkA.end()) return false; };
	return true;
}

template<typename T>
bool Edge<T>::isTopologyPreserving() 
//! checks if edge contraction is preserves topology.
{	if (vertices[0]->order > 2 || vertices[1]->order > 2) return false;
	const unsigned int degree = getDegree(); assert(degree > 0);
	if(degree == 2) {						// accelerated path
		if (vertices[0]->order == 0 && vertices[1]->order == 0)	// special case 2: deg ab = 2, ord ab = ord a = ord b = 0
			return isTopologyPreservingAccelerated();
		if (vertices[0]->order == 0) {				// case 2: deg ab = 2, ord ab = ord a = 0
			Link<T> link(vertices[0]); link.intersect(vertices[1]);
			link.except(this); return link.empty(); }
		else if(vertices[1]->order == 0) {			// case 2: deg ab = 2, ord ab = ord b = 0
			Link<T> link(vertices[1]); link.intersect(vertices[0]);
			link.except(this); return link.empty(); }
		else return false; }					// no case
	else {	if (vertices[0]->order == 1) {				// case 3: ord ab = ord a = 1
			Vertex<T> dummy; Link<T> link(vertices[0], &dummy);
			link.intersect(vertices[1], &dummy); link.except(this, &dummy);
			return (link.empty() && isBoundaryLinkEmpty()); }
		else if (vertices[1]->order == 1) {			// case 3: ord ab = ord b = 1
			Vertex<T> dummy; Link<T> link(vertices[1], &dummy);
			link.intersect(vertices[0], &dummy); link.except(this, &dummy);
			return (link.empty() && isBoundaryLinkEmpty()); }
		else return false; }	 				// no case
}

template<typename T>
bool Edge<T>::isConsistent()
//! checks if edge contraction is consistent.
{	if (vertices[0]->point == vertices[1]->point) return true;
	for (const auto& t : vertices[0]->triangles)
		if (!t->isConsistent(vertices[0], vertices[1], np)) return false;
	for (const auto& t : vertices[1]->triangles)
		if (!t->isConsistent(vertices[1], vertices[0], np)) return false;
	return true;
}

template<typename T>
bool Edge<T>::isCoherentOriented() const
//! checks if edge contraction retains the orientation.
{	unsigned int degree = 0; Triangle<T>* tr[2];
	for (const auto& t : vertices[0]->triangles)
		if (t->vertices[0] == vertices[1] ||
		    t->vertices[1] == vertices[1] ||
		    t->vertices[2] == vertices[1]) {
			if (degree < 2) tr[degree] = t;
			degree++; }
	if (degree < 2) return true; else if (2 < degree) return false;
	unsigned int i = 0, j = 0;
	while (tr[0]->vertices[i] != vertices[0]) ++i;
	while (tr[1]->vertices[j] != vertices[0]) ++j;
	if (tr[0]->vertices[(i+1)%3]  == vertices[1])
		return (tr[1]->vertices[(j+1)%3] != vertices[1]);
	else return (tr[1]->vertices[(j+1)%3] == vertices[1]);
}

template<typename T>
void Edge<T>::setCost(bool check)
//! computes the cost of an edge contraction.
{	Quadric<T> q = vertices[0]->quadric; q += vertices[1]->quadric;
	if (vertices[0]->point == vertices[1]->point) np = vertices[0]->point;
	else if (vertices[0]->order == vertices[1]->order) {
		if (check) q.placeNewVertex3(this);
		else if (vertices[0]->order == 0) {
			if (!q.placeNewVertex1(this) && !q.placeNewVertex2(this))
				q.placeNewVertex3(this); }
		else if (!q.placeNewVertex2(this)) q.placeNewVertex3(this); }
	else np = (vertices[0]->order > vertices[1]->order)? vertices[0]->point: vertices[1]->point;
	cost = q(np); if (cost < 0 || checkEPS(q.area()) == false) cost = 0;
	const T alpha = 0.025f; const T beta = T(M_PI/6.0); T angle = T(M_PI);
	vec3<T> uv = vertices[1]->point-vertices[0]->point; cost = std::sqrt(cost+alpha*norm2(uv));
	for (const auto& t : vertices[0]->triangles) {
		if (t->vertices[0] == vertices[1] || t->vertices[1] == vertices[1] || t->vertices[2] == vertices[1]) {
			unsigned int j = 2;
			if (t->vertices[0] != vertices[0] && t->vertices[0] != vertices[1]) j = 0;
			else if (t->vertices[1] != vertices[0] && t->vertices[1] != vertices[1]) j = 1;
			T a = t->getAngle(j, (j+1)%3, (j+2)%3);
			if (a > T(0.5f*M_PI)) a = std::abs(T(M_PI-a));
			if (a < angle) angle = a; } };
	if (angle < beta) cost *= angle/beta;
}

//! Represents a vertex in mesh simplification.

template<typename T>
struct Vertex {
	vec3<T>		point;				//!< position of this vertex
	unsigned int	number;				//!< vertex number
	int		order;				//!< vertex order
	unsigned int	id;				//!< vertex id in base mesh
	Quadric<T>	quadric;			//!< quadric related to this vertex
	std::deque<Edge<T>*> edges;			//!< attached edges
	std::deque<Triangle<T>*> triangles;		//!< attached triangles
	
	Vertex(vec3<T> p = 0, unsigned int i = 0)					//! allocates an empty vertex record.
		 : point(p), number(0), order(-1), id(i), quadric(), edges(), triangles() { }
	void generate(Edge<T>* const e)							//! adds an edge to this vertex.
		{ edges.push_back(e); }
	void generate(Triangle<T>* const t)						//! adds a triangle to this vertex.
		{ triangles.push_back(t); }
	void setRef(Edge<T>* const en, Edge<T>* const eo)				//! replaces edge eo by en.
		{ for (auto& e : edges) { if (e == eo) { e = en; break; } } }
	void setRef(Triangle<T>* const tn, Triangle<T>* const to)			//! replaces triangle to by tn.
		{ for (auto& t : triangles) { if (t == to) { t = tn; break; } } }
	void setRef(Vertex<T>* const vn, Vertex<T>* const vo);
	void remove(Edge<T>* const e);
	void remove(Triangle<T>* const t);
	void setOrder();
};

template<typename T>
void Vertex<T>::setRef(Vertex<T>* const vn, Vertex<T>* const vo) 
//! replaces vertex vo by vn in neighborhood.
{	for (auto& e : edges) {
		if (e->vertices[0] == vo) e->vertices[0] = vn;
		else if (e->vertices[1] == vo) e->vertices[1] = vn; }
	for (auto& t : triangles) {
		if (t->vertices[0] == vo) t->vertices[0] = vn;
		else if (t->vertices[1] == vo) t->vertices[1] = vn;
		else if (t->vertices[2] == vo) t->vertices[2] = vn; }
}

template<typename T>
void Vertex<T>::remove(Edge<T>* const e) 
//! removes edge e from neighborhood.
{	assert(edges.size()); unsigned int i = 0;
	while (e != edges[i]) { ++i; assert(i < edges.size()); }
	if (i != edges.size()-1) std::swap(edges[i], edges.back());
	edges.pop_back();
}

template<typename T>
void Vertex<T>::remove(Triangle<T>* const t) 
//! removes triangle t from neighborhood.
{	assert(triangles.size()); unsigned int i = 0;
	while (t != triangles[i]) { ++i; assert(i < triangles.size()); }
	if (i != triangles.size()-1) std::swap(triangles[i], triangles.back());
	triangles.pop_back();
}

template<typename T>
void Vertex<T>::setOrder()
//! computes the order of this vertex. 
{	assert(triangles.size()); Link<T> link(this);
	unsigned int k = 0; LinkVertex<T>* boundary[2];
	for (auto& v : link.vertices)
		if (v.number != 2) { if (k < 2) boundary[k] = &v; ++k; }
	if (k == 0) {
		const bool tmp = link.isHalfDisk(0, 0, link.edges[0].vertices[0]);
		if (tmp && link.emptyEdges()) order = 0; else order = 2; }
	else if (k == 2 && boundary[0]->number == boundary[1]->number) {
		unsigned int i = 0, j = 0; bool tmp = true;
		while (i < link.edges.size() && j < boundary[0]->number && tmp) {
			if (link.edges[i].isElement == true) {
				if (link.edges[i].vertices[0] == boundary[0]->vertex) {
					tmp = link.isHalfDisk(i, 0, boundary[1]->vertex); ++j; }
				else if(link.edges[i].vertices[1] == boundary[0]->vertex) {
					tmp = link.isHalfDisk(i, 1, boundary[1]->vertex); ++j; } }
			++i; }
		order = (tmp && link.emptyEdges())? 1: 2; }
	else order = 2;
}

//! Represents a triangle in mesh simplification.

template<typename T> struct Triangle {
	Vertex<T>* vertices[3];				//!< vertices of this triangle
	
	Triangle(Vertex<T> *v0, Vertex<T> *v1, Vertex<T> *v2)				//! allocates a record for vertices (v0,v1,v2).
		{ assert(v0 && v1 && v2); vertices[0] = v0; vertices[1] = v1; vertices[2] = v2; }
	void	resize(Vertex<T> *vn, Vertex<T> *vo);
	T	getAngle(unsigned int i, unsigned int j, unsigned int k) const;
	bool	isConsistent(Vertex<T> *v0, Vertex<T> *v1, vec3<T> &pn);
	T	getArea() const; 
};

template<typename T>
void Triangle<T>::resize(Vertex<T> *vn, Vertex<T> *vo)
//! replaces vertex vo by vn.
{	unsigned int i = 2;
	if (vertices[0] == vo) i = 0;
	else if (vertices[1] == vo) i = 1; assert(vertices[i] == vo);
	vertices[i]->remove(this); vertices[i] = vn; vertices[i]->generate(this);
}

template<typename T>
T Triangle<T>::getAngle(unsigned int i, unsigned int j, unsigned int k) const 
//! computes angle at i.
{	assert(i < 3 && j < 3 && k < 3 && i != j && i != k && j != k);
	vec3<T> ba = vertices[j]->point-vertices[i]->point, ca = vertices[k]->point-vertices[i]->point;
	double tmp = dot(ba,ba)*dot(ca,ca);
	if (checkEPS(tmp)) tmp = dot(ba,ca)/std::sqrt(tmp); else tmp = 1.0;
	if (tmp < -1.0) tmp = -1.0; else if (1.0 < tmp) tmp = 1.0;
	return T(acos(tmp));
}

template<typename T>
bool Triangle<T>::isConsistent(Vertex<T> *v0, Vertex<T> *v1, vec3<T> &pn) 
//! checks if edge contraction at (v0,v1) is consistent.
{	unsigned int i;
	if (vertices[0] == v0) i = 0; 
	else if (vertices[1] == v0) i = 1; 
	else i = 2; assert(vertices[i] == v0);
	if (vertices[(i+1)%3] == v1 || vertices[(i+2)%3] == v1) return true;
	if (vertices[(i+1)%3]->point == vertices[(i+2)%3]->point) return true;
	vec3<T> p0 = pn; vec3<T> a(vertices[i]->point);
	vec3<T> b(vertices[(i+1)%3]->point); vec3<T> c(vertices[(i+2)%3]->point);
	vec3<T> np(p0), v21 = c-b, v01 = a-b, npv1 = np-b, n1 = cross(v21,v01); n1 = cross(n1,v21);
	double ret = dot(n1,npv1); return ret > 0;
}

template<typename T>
T Triangle<T>::getArea() const 
//! returns triangle area.
{	vec3<T> u = vertices[1]->point-vertices[0]->point, v = vertices[2]->point-vertices[0]->point;
	return 0.5f*norm(cross(u,v));
}

//! Implements a spatial data structure for checking collisions in mesh simplification.

template<typename T, typename TX, typename GP> class collisionCache : public spatialCache<GP> {
	void	triangleBounds(const GP t);
	std::list<Vertex<T>*> collectVertices();
	std::list<Edge<T>*> collectEdges();
	bool	hit1EdgeTri(const vec3<T> &p, const vec3<T> &q, TX &tmin);
	bool	hit1TriEdge(const vec3<T> &a, const vec3<T> &b, const vec3<T> &c, TX &tmin);
	bool	hit1TetraVtx(const vec3<T> &a, const vec3<T> &b, const vec3<T> &c, const vec3<T> &d, TX &tmin);
	bool	testHit1(const Vertex<T> *v0, const vec3<T> &np, const Vertex<T> *v1, TX &tmin);
	bool	hit2TriVtx(const Edge<T> *e, TX &tmin);
	bool	cont1VtxEdge(const Edge<T> *e, TX &tmin);
	bool	cont1VtxTri(const Edge<T> *e, TX &tmin);
	bool	cont1EdgeEdge(const Edge<T> *e, TX &tmin);
public:
	collisionCache(vec3<T>& _min, vec3<T>& _max, T len)				//! allocates an empty cache with bounds (min,max) and average edge length len.
		 : spatialCache<GP>(_min, _max, len) { }
	~collisionCache() { }
	bool	test(const Edge<T>* const e);
	bool	test(const vec3<T> &p, const vec3<T> &q);
};

template<typename T, typename TX, typename GP>
bool collisionCache<T,TX,GP>::test(const vec3<T> &p, const vec3<T> &q)
//! checks edge p-q for collision.
{	TX tmin; if (hit1EdgeTri(p, q, tmin) == false) return false;
	printf("collision at p %f %f %f q %f %f %f\n", p.x, p.y, p.z, q.x, q.y, q.z);
	hit1EdgeTri(p, q, tmin); return true;						// test once again, for debugging
}

template<typename T, typename TX, typename GP>
void collisionCache<T,TX,GP>::triangleBounds(const GP t)
//! determines bounding box for a triangle.
{	this->boundTo3(t->vertices[0]->point, t->vertices[1]->point, t->vertices[2]->point);
	this->setBounds();
}

template<typename T, typename TX, typename GP>
std::list<Vertex<T>*> collisionCache<T,TX,GP>::collectVertices()
// makes a list of all vertices in the bounding box.
{	this->setBounds(); std::list<Vertex<T>*> vl;
	for (unsigned int iz = this->izmin; iz <= this->izmax; iz++)  {
		for (unsigned int iy = this->iymin; iy <= this->iymax; iy++)  {
			for (unsigned int ix = this->ixmin; ix <= this->ixmax; ix++)  {
				const unsigned int i = this->index(ix, iy, iz);
				if (this->data == nullptr || this->data[i] == nullptr) continue; // no simplices in this box
				for (const auto& li : *(this->data[i])) {
					for (unsigned int j = 0; j < 3; j++)  {
						Vertex<T>* vi = li->vertices[j]; bool f = false;
						for (const auto v : vl) if (v == vi) { f = true; break; };
						if (f == false) vl.push_front(vi); } } } } };
	return vl;
}

template<typename T, typename TX, typename GP>
std::list<Edge<T>*> collisionCache<T,TX,GP>::collectEdges()
// makes a list of all edges in the bounding box.
{	this->setBounds(); std::list<Edge<T>*> el;
	for (unsigned int iz = this->izmin; iz <= this->izmax; iz++)  {
		for (unsigned int iy = this->iymin; iy <= this->iymax; iy++)  {
			for (unsigned int ix = this->ixmin; ix <= this->ixmax; ix++)  {
				unsigned int i = this->index(ix, iy, iz);
				if (this->data == nullptr || this->data[i] == nullptr) continue; // no simplices in this box
				for (const auto& li : *(this->data[i])) {
					for (unsigned int j = 0; j < 3; j++)  {
						auto ei = new Edge<T>(li->vertices[j], li->vertices[(j+1)%3]); bool f = false;
						for (const auto e : el) if (e == ei) { f = true; break; };
						if (f == false) el.push_front(ei); else delete ei; } } } } };
	return el;
}

template<typename T, typename TX, typename GP>
bool collisionCache<T,TX,GP>::hit1EdgeTri(const vec3<T>& e1, const vec3<T>& e2, TX &tmin)
//! tests hit collision 1st kind: edge-triangle intersection.
{	tmin = 1; bool collision = false; this->boundTo2(e1, e2); 
	const vec3<TX> p(e1); vec3<TX> q(e2); vec3<TX> i;
	const auto tl = this->collectTriangles();					// make a list of triangles in the bounding box of p-q
	for (const auto t : tl)  {							// check edge-triangle intersection for all cells in the bounding box
		const vec3<TX> a(t->vertices[0]->point), b(t->vertices[1]->point), c(t->vertices[2]->point);
		if (p == a || p == b || p == c || q == a || q == b || q == c) continue;
		if (edgeTriangleIntersection(a, b, c, p, q, i) == false) continue;
		if (edgeTriangleIntersection(a, b, c, p, q, i) == false) continue;
		const TX ti = edgeParameter(p, q, i);
		if (ti < tmin) tmin = ti;
		collision = true; if (tmin == 0) break; };
	return collision;
}

template<typename T, typename TX, typename GP>
bool collisionCache<T,TX,GP>::hit1TriEdge(const vec3<T>& t1, const vec3<T>& t2, const vec3<T>& t3, TX &tmin)
//! tests hit collision 1st kind: triangle-edge intersection.
{	tmin = 1; bool collision = false; this->boundTo3(t1, t2, t3); 
	const vec3<TX> a(t1), b(t2), c(t3); vec3<TX> i;
	const auto el = collectEdges();							// make a list of edges in the bounding box of a-b-c
	for (const auto e : el) {							// check edge-triangle intersection for all cells in the bounding box
		const vec3<TX> p(e->vertices[0]->point), q(e->vertices[1]->point);
		if (p == a || p == b || p == c || q == a || q == b || q == c) continue;
		if (triangleEdgeIntersection(a, b, c, p, q, i) == false) continue;
		const TX ti = baryBound(b, a, c, i); if (ti < tmin) tmin = ti;
		collision = true; if (tmin == 0) break; };
	for (auto& e : el) delete e;
	return collision;
}

template<typename T, typename TX, typename GP>
bool collisionCache<T,TX,GP>::hit1TetraVtx(const vec3<T> &t1, const vec3<T> &t2, const vec3<T> &t3, const vec3<T> &t4, TX &tmin)
//! test hit collision 1st kind: tetra-vertex intersection.
{	tmin = 1; bool collision = false; this->boundTo4(t1, t2, t3, t4); 
	const vec3<TX> a(t1), b(t2), c(t3), d(t4);
	const auto vl = collectVertices();						// make a list of vertices in the bounding box of a-b-c-d
	for (const auto v : vl) { const vec3<TX> p(v->point);				// check point-tetra inclusion for all cells in the bounding box
		if (p == a || p == b || p == c || p == d) continue;
		if (pointInTetra(a, b, c, d, p) == true) { collision = true; break; } };
	return collision;
}

template<typename T, typename TX, typename GP>
bool collisionCache<T,TX,GP>::testHit1(const Vertex<T>* vi, const vec3<T> &np, const Vertex<T>* vj, TX &tmin)
//! tests hit collision 1st kind: all one-sided intersections.
{	TX t = std::numeric_limits<float>::max(); tmin = 1; bool collision = false;
	if (hit1EdgeTri(vi->point, np, t))  {						// test edge - triangle intersections
		if (t < tmin) tmin = t;
		collision = true; if (tmin == 0) return collision; };
	for (const auto es : vi->edges)  {						// test triangle - edge intersections
		const unsigned int s = es->vertices[0] == vi? 1 : 0;
		Vertex<T> *vs = es->vertices[s]; if (vs->point == vj->point) continue;
		if (hit1TriEdge(vi->point, np, vs->point, t))  {
			if (t < tmin) tmin = t;
			collision = true; if (tmin == 0) return collision; } };
	for (const auto vs : vi->triangles)  {						// test tetra - vertex intersections
		Vertex<T> **tv = vs->vertices;
		if (tv[0]->point == vj->point || tv[1]->point == vj->point || tv[2]->point == vj->point) continue;
		if (hit1TetraVtx(tv[0]->point, tv[1]->point, tv[2]->point, np, t))  {
			if (t < tmin) tmin = t;
			collision = true; if (tmin == 0) return collision; } };
	return collision;
}

template<typename T, typename TX, typename GP>
bool collisionCache<T,TX,GP>::hit2TriVtx(const Edge<T> *e, TX &tmin)
//! test hit collisions 2nd kind: triangle-vertex intersection. 
{	tmin = 1; bool collision = false; Vertex<T> *v0 = e->vertices[0]; Vertex<T> *v1 = e->vertices[1];
	for (const auto t : v0->triangles) {
		if (t->vertices[0] != v1 && t->vertices[1] != v1 && t->vertices[2] != v1) continue;		
		this->boundTo4(t->vertices[0]->point, t->vertices[1]->point, t->vertices[2]->point, e->np);
		const vec3<TX> a(t->vertices[0]->point), b(t->vertices[1]->point), c(t->vertices[2]->point); vec3<TX> d(e->np);
		const auto vl = collectVertices();
		for (const auto v : vl) { const vec3<TX> p(v->point);			// for a list of vertices in the bounding box of a-b-c-vc
			if (p == a || p == b || p == c || p == d) continue;
			if (pointInTetra(a, b, c, d, p)) { collision = true; return collision; } } };
	return collision;
}

template<typename T, typename TX, typename GP>
bool collisionCache<T,TX,GP>::cont1VtxEdge(const Edge<T> *e, TX &tmin)
//! test contraction collision 1st kind: vertex-edge intersection.
{	tmin = 1; bool collision = false; Vertex<T> *v0 = e->vertices[0]; Vertex<T> *v1 = e->vertices[1];
	vec3<TX> a(v0->point); vec3<TX> b(v1->point); vec3<TX> c(e->np); vec3<TX> n, i;
	if (normal(a, b, c, n) == false) return false;					// triangle plane normal
	for (const auto es : v0->edges) {						// for all edges starting at v0...
		const unsigned int s = es->vertices[0] == v0 ? 1 : 0;
		Vertex<T>* const vs = es->vertices[s]; if (vs == v1) continue;		// find the other vertex to form edge v0-vs
		vec3<TX> p(vs->point); if (checkEPS(planeDist(a,n,p))) continue;	// p is not in plane
		if (a == b || a == c || p == b || p == c) continue;
		if (edgeEdgeIntersection(a, p, b, c, i) == false) continue;		// check intersection v0-vs and v1-vc
		TX ti = baryBound(c, b, a, i); 
		if (std::abs(ti) < tmin) tmin = std::abs(ti);
		collision = true; if (tmin == 0) return collision; };
	for (const auto es : v1->edges) {						// for all edges starting at v1...
		const unsigned int s = es->vertices[0] == v1 ? 1 : 0;
		Vertex<T>* const vs = es->vertices[s]; if (vs == v0) continue;		// find the other vertex to form edge v1-vs
		vec3<TX> p(vs->point); if (checkEPS(planeDist(a,n,p))) continue;	// p is not in plane
		if (b == a || b == c || p == a || p == c) continue;
		if (edgeEdgeIntersection(b, p, a, c, i) == false) continue;		// check intersection v1-vs and v0-vc
		TX ti = baryBound(c, a, b, i);
		if (std::abs(ti) < tmin) tmin = std::abs(ti);
		collision = true; if (tmin == 0) return collision; };
	return collision;
}

template<typename T, typename TX, typename GP>
bool collisionCache<T,TX,GP>::cont1VtxTri(const Edge<T> *e, TX &tmin)
//! test contraction collision 1st kind: vertex-triangle intersection.
{	tmin = 1; bool collision = false; unsigned int i1 = 0, i2 = 0; vec3<TX> i;
	Vertex<T> *v0 = e->vertices[0]; Vertex<T> *v1 = e->vertices[1];
	vec3<TX> p(v0->point); vec3<TX> q(v1->point); vec3<TX> r(e->np);
	for (const auto t : v0->triangles) {						// for all triangles at v0... 
		if (t->vertices[0] == v1 || t->vertices[1] == v1 || t->vertices[2] == v1) continue;
		if (t->vertices[0] == v0) { i1 = 1; i2 = 2; }				// find the other vertices
		else if (t->vertices[1] == v0) { i1 = 0; i2 = 2; }
		else if (t->vertices[2] == v0) { i1 = 0; i2 = 1; };
		vec3<TX> b(t->vertices[i1]->point), c(t->vertices[i2]->point);
		if (q == p || q == b || q == c || r == p || r == b || r == c) continue;
		if (edgeTriangleIntersection(p, b, c, q, r, i) == false) continue;	// test intersection of triangle v0-vs-vs' with edge v1-vc
		TX ti = baryBound(r, q, p, i); 
		if (std::abs(ti) < tmin) tmin = std::abs(ti);
		collision = true; if (tmin == 0) return collision; };
	for (const auto t : v1->triangles) { 						// for all triangles at v1... 
		if (t->vertices[0] == v0 || t->vertices[1] == v0 || t->vertices[2] == v0) continue;
		if (t->vertices[0] == v1) { i1 = 1; i2 = 2; }				// find the other vertices
		else if (t->vertices[1] == v1) { i1 = 0; i2 = 2; }
		else if (t->vertices[2] == v1) { i1 = 0; i2 = 1; };
		vec3<TX> b(t->vertices[i1]->point), c(t->vertices[i2]->point);
		if (p == q || p == b || p == c || r == q || r == b || r == c) continue;
		if (edgeTriangleIntersection(q, b, c, p, r, i) == false) continue; 	// test intersection of triangle v1-vs-vs' with edge v0-vc
		TX ti = baryBound(r, p, q, i);
		if (std::abs(ti) < tmin) tmin = std::abs(ti);
		collision = true; if (tmin == 0) return collision; };
	return collision;
}

template<typename T, typename TX, typename GP>
bool collisionCache<T,TX,GP>::cont1EdgeEdge(const Edge<T> *e, TX &tmin)
//! tests contraction collision 1st kind: edge-edge intersection.
{	tmin = 1; bool collision = false; Vertex<T> *v0 = e->vertices[0]; Vertex<T> *v1 = e->vertices[1];
	vec3<TX> a(e->np); vec3<TX> b(v0->point); vec3<TX> p(v1->point); vec3<TX> i0, i1, n0, n1;
	if (checkEPS(area(a,b,p)) == false) return false; 				// degenerated case
	for (const auto es0 : v0->edges) {						// for all edges starting at v0...
		const unsigned int s0 = es0->vertices[0] == v0 ? 1 : 0;
		const Vertex<T> *vs0 = es0->vertices[s0]; if (vs0 == v1) continue;	// find the other vertex to form triangle vc-v0-vs0
		const vec3<TX> q(vs0->point); if (normal(a, b, q, n0) == false) continue; // triangle plane normal
		for (const auto es1 : v1->edges) {					// for all edges starting at v1... 
			const unsigned int s1 = es1->vertices[0] == v1 ? 1 : 0;
			const Vertex<T> *vs1 = es1->vertices[s1]; if (vs1 == v0) continue;
			const vec3<TX> r(vs1->point);					// find the other vertex to form triangle vc-v1-vs1
			if (normal(a, p, r, n1) == false) continue;			// triangle plane normal
			const vec3<TX> n = cross(n1,n0);
			if (checkEPS(dot(n,n)) == false) continue;			// parallel planes
			if (b == r || b == p || q == p || q == r) continue;
			if (edgeLineIntersection(a, n, b, q, i0) == false) continue;	// compute intersection of ABQ & APR
			if (edgeLineIntersection(a, n, p, r, i1) == false) continue;
			TX b0 = bary(a, b, q, i0); if (b0 < 0 || b0 > 1) continue;
			TX b1 = bary(a, p, r, i1); if (b1 < 0 || b1 > 1) continue;
			TX l0 = pointDist(a, i0); TX l1 = pointDist(a, i1);		// compute collision time
			TX den = b1*l0-b0*l1; if (den == 0) continue;
			TX ti = (l0-l1)/den; if (ti <= 0 || ti >= 1) continue;
			if (ti < tmin) tmin = ti;
			collision = true; if (tmin == 0) return collision; } };
	return collision;
}

template<typename T, typename TX, typename GP>
bool collisionCache<T,TX,GP>::test(const Edge<T> *e)
//! test collision of Edge e with mesh.
{	Vertex<T> *v0 = e->vertices[0]; Vertex<T> *v1 = e->vertices[1];
	TX t = std::numeric_limits<float>::max();
	if (testHit1(v0, e->np, v1, t)) return true;
	if (testHit1(v1, e->np, v0, t)) return true;
	if (hit1TriEdge(v0->point, e->np, v1->point, t)) return true;
	if (hit2TriVtx(e, t)) return true;
	if (cont1VtxEdge(e, t)) return true;
	if (cont1VtxTri(e, t)) return true;
	if (cont1EdgeEdge(e, t)) return true;
	return false;
}

//! Implements a context for simplification of triangular meshes by edge contraction.

template<typename T, typename TX> class SimplicialComplex
{
	bool	isCoherentOriented() const						//! checks that mesh is properly oriented.
		{ for (const auto& e : edges) {
			if (!e.isCoherentOriented()) return false; };
		  return true; }
	T	edgeLength()								//! computes average edge length.
		{ const auto n = edges.size(); T el = 0;
		  for (const auto& e : edges) el += e.length();
		  return n? el/n: el; }
	void	setCache();
	bool	contract(Edge<T>* const e);
	unsigned int countIntersections();
	void	remove(Vertex<T>* const v);
	void	remove(Edge<T>* const e);
	void	remove(Triangle<T>* const t);
	void	removeEdges(Vertex<T>* const v0, Vertex<T>* const v1);
	void	resizeTriangles(Vertex<T>* const v0, Vertex<T>* const v1, vec3<T> const &np);
	void	setEdge(Vertex<T>* const v0, Vertex<T>* const v1);
	void	update(Vertex<T>* const v0);
	void	updateEdge(Edge<T>* const e);
	SimplicialComplex(const SimplicialComplex& s) = delete;
	SimplicialComplex& operator=(const SimplicialComplex& s) = delete;
	SimplicialComplex& operator=(SimplicialComplex&& s) = delete;
public:
	std::deque<Vertex<T>> vertices;			//!< set of vertices
	std::deque<Edge<T>> edges;			//!< set of edges
	std::deque<Triangle<T>> triangles;		//!< set of triangles
	std::deque<cmd>* history;			//!< simplification history
	Heap<Edge<T>> heap;				//!< contains the sequence of edges considered for contraction
	collisionCache<T,TX,Triangle<T>* > *cc;		//!< used to check for edge collisions
	bool	check;					//!< flags that collision check is wanted
	unsigned int ntri;				//!< target number of triangles
	
	SimplicialComplex(triangleMesh& m, bool _check, cmdlist* h = nullptr);
	~SimplicialComplex() { delete cc; }
	void	simplify(const unsigned int nf, const bool verbose = false);
	void	saveCompact(triangleMesh& m);
	void	save(triangleMesh& m);
};

template<typename T, typename TX>
SimplicialComplex<T,TX>::SimplicialComplex(triangleMesh& m, bool _check, cmdlist* h)
//! allocates a context for simplifying triangle mesh m, including checks, and recording the simplification operations.
	 : vertices(), edges(), triangles(), history(h), heap(), cc(0), check(_check), ntri(0)
{	const unsigned int nv = m.nVertices(); uvecD ids(nv+1);
	ids = std::numeric_limits<unsigned int>::max(); unsigned int id = 0;
	for (unsigned int i = 0; i < nv; i++) { const vertex *v = m.vertices(i);
		if (v) { vertices.push_back(Vertex<T>(v->pos,i)); ids[i] = id++; }; };
	for (unsigned int i = 0; i < m.nPrimitives(); i++) {
		const triangle* t = m.tr(i); if (t == nullptr) continue;
		triangles.push_back(Triangle<T>(&vertices[ids[t->vtx(0)]],
				&vertices[ids[t->vtx(1)]], &vertices[ids[t->vtx(2)]])); }
	for (unsigned int i = 0; i < triangles.size(); i++) {
		triangles[i].vertices[0]->triangles.push_back(&triangles[i]);
		triangles[i].vertices[1]->triangles.push_back(&triangles[i]);
		triangles[i].vertices[2]->triangles.push_back(&triangles[i]); };
	for (unsigned int i = 0; i < vertices.size(); i++) {
		if (vertices[i].triangles.size() == 0) remove(&vertices[i]); };
	for (unsigned int i = 0; i < triangles.size(); i++) {
		setEdge(triangles[i].vertices[0], triangles[i].vertices[1]);
		setEdge(triangles[i].vertices[1], triangles[i].vertices[2]);
		setEdge(triangles[i].vertices[2], triangles[i].vertices[0]); };
	for (unsigned int i = 0; i < vertices.size(); i++) {
		if (vertices[i].order < 3) vertices[i].setOrder(); };
	for (unsigned int i = 0; i < ITOU(triangles.size()); i++) {
		const Quadric<T> q(triangles[i].vertices[0]->point, triangles[i].vertices[1]->point, triangles[i].vertices[2]->point);
		const T w0 = T(triangles[i].getAngle(0, 1, 2)*q.area()/M_PI);
		const T w1 = T(triangles[i].getAngle(1, 2, 0)*q.area()/M_PI);
		const T w2 = T(triangles[i].getAngle(2, 0, 1)*q.area()/M_PI);
		triangles[i].vertices[0]->quadric += q*w0;
		triangles[i].vertices[1]->quadric += q*w1;
		triangles[i].vertices[2]->quadric += q*w2; };
	for (unsigned int i = 0; i < edges.size(); i++) edges[i].setCost(check);
	if (check) setCache();
}

template<typename T, typename TX>
void SimplicialComplex<T,TX>::saveCompact(triangleMesh& m)
//! saves remaining vertices and triangles in a mesh - compact version.
{	m.vertices.clearNodes(); m.primitives.clearNodes();
	for (unsigned int i = 0; i < vertices.size(); i++) { 
		vertex* v = new vertex(vertices[i].point);
		m.addVertex(v,i); vertices[i].number = i; };
	for (unsigned int i = 0; i < triangles.size(); i++) {
 		triangle* t = new triangle(triangles[i].vertices[0]->number,
			triangles[i].vertices[1]->number, triangles[i].vertices[2]->number);
   		m.addTriangle(t,i); }
}

template<typename T, typename TX>
void SimplicialComplex<T,TX>::save(triangleMesh& m)
// saves remaining vertices and triangles in a mesh - keeps original ids.
{	m.vertices.clearNodes(); m.primitives.clearNodes();
	for (unsigned int i = 0; i < vertices.size(); i++) {
		vertex* v = new vertex(vertices[i].point); m.addVertex(v, vertices[i].id); };
	for (unsigned int i = 0; i < triangles.size(); i++) {
 		triangle* t = new triangle(triangles[i].vertices[0]->id,
				triangles[i].vertices[1]->id, triangles[i].vertices[2]->id);
   		m.addTriangle(t,i); }
}

template<typename T, typename TX>
void SimplicialComplex<T,TX>::setCache()
//! setup collision cache for this mesh.
{	vec3<T> bmin, bmax; bmin = HUGE; bmax = -HUGE;
	for (unsigned int i = 0; i < vertices.size(); i++) {				// compute bounding box
		vec3<T> p = vertices[i].point;
		if (p.x < bmin.x) bmin.x = p.x;
		if (p.y < bmin.y) bmin.y = p.y;
		if (p.z < bmin.z) bmin.z = p.z;
		if (p.x > bmax.x) bmax.x = p.x;
		if (p.y > bmax.y) bmax.y = p.y;
		if (p.z > bmax.z) bmax.z = p.z;	};
	cc = new collisionCache<T,TX,Triangle<T>* >(bmin, bmax, edgeLength());		// set up spatial cache
	for (unsigned int i = 0; i < triangles.size(); i++)				// add triangles to spatial cache
		cc->addTriangle(&triangles[i]);
	ntri = ITOU(triangles.size()); 
	if (countIntersections()) throw rtException("Mesh has intersections");		// cache must be free of intersections
}

template<typename T, typename TX>
unsigned int SimplicialComplex<T,TX>::countIntersections()
//! returns the number of intersections in this mesh.
{	unsigned int n = 0;
	for (unsigned int i = 0; i < edges.size(); i++)
		n += cc->test(edges[i].vertices[0]->point, edges[i].vertices[1]->point);
	return n;
}

template<typename T, typename TX>
bool SimplicialComplex<T,TX>::contract(Edge<T>* const e) 
//! contracts edge e.
{	Vertex<T>* v0 = e->vertices[0]; Vertex<T>* v1 = e->vertices[1];
	if (check && cc->test(e)) { heap.remove(e); return false; };					// check if new edge intersects mesh
	if (history) history->push_back(cmd('c', e->vertices[1]->id));
	resizeTriangles(v0, v1, e->np);
	if (v0->order < v1->order) { v0->number = v1->number; v0->order = v1->order; }
	v0->quadric += v1->quadric; remove(e); removeEdges(v0, v1);
	if (v0 == &vertices.back()) { remove(v1); std::swap(v0, v1); } else remove(v1); update(v0);
	if (check) { if (triangles.size() < ntri/2) { cc->clear(); setCache(); }; };			// update spatial cache
	return true;
}

template<typename T, typename TX>
void SimplicialComplex<T,TX>::remove(Vertex<T>* const v) 
//! removes vertex v.
{	assert(v->edges.size() == 0 && v->triangles.size() == 0);
	if (history) history->push_back(cmd('a', v->id));
	if (v != &vertices.back()) { std::swap(*v, vertices.back()); v->setRef(v, &vertices.back()); }
	vertices.pop_back();
}

template<typename T, typename TX>
void SimplicialComplex<T,TX>::remove(Edge<T>* const e) 
//! removes edge e.
{	e->vertices[0]->remove(e); e->vertices[1]->remove(e); heap.remove(e);
	if (history) history->push_back(cmd('l', e->vertices[0]->id, e->vertices[1]->id));
	if (e != &edges.back()) { std::swap(*e, edges.back()); heap.setRef(e->heap, e);
		e->vertices[0]->setRef(e, &edges.back()); e->vertices[1]->setRef(e, &edges.back()); }
	edges.pop_back();
}

template<typename T, typename TX>
void SimplicialComplex<T,TX>::remove(Triangle<T>* const t) 
//! removes triangle t.
{	if (check) cc->removeTriangle(t);
	t->vertices[0]->remove(t); t->vertices[1]->remove(t); t->vertices[2]->remove(t);
	if (t != &triangles.back()) {
		if (check) cc->removeTriangle(&triangles.back());
		std::swap(*t, triangles.back());
		t->vertices[0]->setRef(t, &triangles.back()); t->vertices[1]->setRef(t, &triangles.back()); t->vertices[2]->setRef(t, &triangles.back());
		if (check) cc->addTriangle(t);
	}
	triangles.pop_back();
}

template<typename T, typename TX>
void SimplicialComplex<T,TX>::removeEdges(Vertex<T>* const v0, Vertex<T>* const v1) 
//! removes edge (v0,v1).
{	for (unsigned int i = 0; i < v1->edges.size(); i++) {
		unsigned int k1 = v1->edges[i]->vertices[0] != v1? 0 : 1; unsigned int j = 0;
		while (j < v0->edges.size()) { const unsigned int k2 = v0->edges[j]->vertices[0] != v0? 0 : 1;
			if (v1->edges[i]->vertices[k1] == v0->edges[j]->vertices[k2]) { remove(v1->edges[i]); --i; break; }
			else ++j; } }
	while (v1->edges.size() != 0) {
		if (history) history->push_back(cmd('l', v1->edges[0]->vertices[0]->id, v1->edges[0]->vertices[1]->id));
		int k = v1->edges[0]->vertices[0]->id == v1->id? 1: 0;
		if (history) history->push_back(cmd('u', v1->edges[0]->vertices[k]->id, v0->id));
		v1->edges[0]->resize(v0, v1); }						// resize the other edges
}

template<typename T, typename TX>
void SimplicialComplex<T,TX>::resizeTriangles(Vertex<T>* const v0, Vertex<T>* const v1, vec3<T> const &np)
//! replaces position np in vertex v0.
{	unsigned int i = 0;
	while (i < v1->triangles.size()) {
		if (v1->triangles[i]->vertices[0] == v0 ||
		    v1->triangles[i]->vertices[1] == v0 ||
		    v1->triangles[i]->vertices[2] == v0) remove(v1->triangles[i]);
		else ++i; }
	if (check) { for (unsigned int i = 0; i < v0->triangles.size(); i++) cc->removeTriangle(v0->triangles[i]); };
	v0->point = np;
	if (check) { for (unsigned int i = 0; i < v0->triangles.size(); i++) cc->addTriangle(v0->triangles[i]); };
	while (v1->triangles.size() != 0) {					// resizing the other tr
		Triangle<T>* t = v1->triangles[0]; if (check) cc->removeTriangle(t);
		v1->triangles[0]->resize(v0, v1); if (check) cc->addTriangle(t); };
}

template<typename T, typename TX>
void SimplicialComplex<T,TX>::setEdge(Vertex<T>* const v0, Vertex<T>* const v1)
//! allocates an edge record for vertices (v0,v1).
{	for (unsigned int i = 0; i < v0->edges.size(); i++) { 
		unsigned int j = v0->edges[i]->vertices[0] != v0? 0: 1;
		if (v0->edges[i]->vertices[j] == v1) return;
	}
	edges.push_back(Edge<T>(v0, v1)); v0->edges.push_back(&edges.back()); v1->edges.push_back(&edges.back());
}

template<typename T, typename TX>
void SimplicialComplex<T,TX>::updateEdge(Edge<T>* const e)
//! updates configuration at edge e.
{	e->setCost(check);
	if (e->heap == 0) { e->tp = e->cp = true; heap.insert(e); }
	else heap.update(e);
}

template<typename T, typename TX>
void SimplicialComplex<T,TX>::update(Vertex<T>* const v)
//! updates configuration at vertex v.
{	for (unsigned int i = 0; i < v->edges.size(); i++) updateEdge(v->edges[i]);
	for (unsigned int i = 0; i < v->triangles.size(); i++) {
		Vertex<T>* tmp[2]; unsigned int j = 0;
		if (v->triangles[i]->vertices[0] != v) tmp[j++] = v->triangles[i]->vertices[0];
		if (v->triangles[i]->vertices[1] != v) tmp[j++] = v->triangles[i]->vertices[1];
		if (v->triangles[i]->vertices[2] != v) tmp[j++] = v->triangles[i]->vertices[2];
		assert(j == 2);
		for (unsigned int j = 0; j < tmp[0]->edges.size(); j++) {
			const unsigned int k = tmp[0]->edges[j]->vertices[0] != tmp[0]? 0: 1;
			if (tmp[0]->edges[j]->vertices[k] == tmp[1]) { updateEdge(tmp[0]->edges[j]); break; } } }
	for (unsigned int i = 0; i < v->edges.size(); i++) {
		const unsigned int j = v->edges[i]->vertices[0] != v? 0 : 1;
		for (unsigned int k = 0; k < v->edges[i]->vertices[j]->edges.size(); k++)
			if (v->edges[i]->vertices[j]->edges[k]->heap == 0)
					updateEdge(v->edges[i]->vertices[j]->edges[k]); }
}

template<typename T, typename TX>
void SimplicialComplex<T,TX>::simplify(const unsigned int nf, const bool verbose)
//! simplifies triangular mesh to target number of faces nf.
{	const bool before = isCoherentOriented();
	for (unsigned int i = 0; i < edges.size(); ++i) {
		Edge<T> *e = &edges[i];	e->tp = e->isTopologyPreserving(); e->cp = e->isConsistent();
		if (e->tp && e->cp) heap.insert(e); else e->heap = 0; };
	while (triangles.size() > nf && heap.size() > 1) {
		Edge<T> *e = heap.getRef(1); e->tp = e->isTopologyPreserving(); e->cp = e->isConsistent();
		if (verbose) { printf("%zd triangles %f\r", triangles.size(), e->cost); fflush(stdout); }
		if (e->tp && e->cp) contract(e);
		else heap.remove(e); };
	const bool after = isCoherentOriented();
	if (before && after == false) throw rtException("Mesh simplification resulted in an orientation exception");
	if (check && countIntersections()) throw rtException("Mesh simplification resulted in mesh intersection");
}

#endif
