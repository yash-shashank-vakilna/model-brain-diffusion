#ifndef ELEMENT_H
#define ELEMENT_H

/*
 *
 * element.h: definitions & classes for finite elements modeling
 * BRIAN Software Package Version 3.0
 *
 * $Id: element.h 504 2017-03-16 23:28:00Z kruggel $
 *
 * 0.10 (18/02/11): first version UH & FK
 * 0.20 (14/05/11): first working FK
 * 0.30 (23/10/11): first release for BRIAN 2.2
 * 0.31 (31/12/12): released version 2.4
 * 0.40 (16/12/13): documented
 * v406 (28/09/16): bumped to version 3.0
 *
 */

/*! \file
    \brief Definitions & classes for finite elements modeling.
*/

#define allDOF(d)		(unsigned int d = 0; d < dof; d++)
#define allNodes(i)		(unsigned int i = 0; i < element<T>::nn; i++)
#define allGPoints(i)		(unsigned int i = 0; i < element<T>::ng; i++)

/*! \enum elementType
    \brief Symbolic constants for element types in volumetric meshes.
*/
enum class elementType {
	hexIso8,					///< 8-node hexahedral element
	tetIso4,					///< 4-node tetrahedral element
	hexIso20,					///< 20-node hexahedral element
	tetIso10					///< 10-node tetrahedral element
};

using uvec = std::vector<unsigned int>;

//! Base class of all finite elements

template <typename T> class element
{	// base class for all finite elements
protected:
	unsigned int nn;			//!< number of nodes
	unsigned int ng;			//!< number of Gauss points
	const vec3<T>* gp;			//!< Gauss points
	const vec3<T>* nd;			//!< Gauss nodes
	vecD<T>	weight;				//!< weights at Gauss points
	vecD<T>	sf;				//!< shape functions
	vecD<T>	du, dv, dw;			//!< derivatives of shape functions
	vecD<T>	dx, dy, dz; 			//!< global derivatives
	const cell& c;				//!< element description from volume mesh
	const volumeMesh& m;			//!< vertex graph
	unsigned int dof;			//!< number of DOFs
	T	det;				//!< determinant of Jacobian
	std::vector<vec3<T>> dp;		//!< vector of displacements
	const std::vector<uvec>* fc;

	vec3<T>	getPos(const unsigned int i) const					//! returns world position of node i.
		{ return m.pos(c.vtx(i)); }
	vecD<T>	getNormal(const uvec& f) const						//! return face normal
		{ const vec3<T> e0 = getPos(f[2])-getPos(f[0]), e1 = getPos(f[1])-getPos(f[0]);
		  const vec3<T> n = cross(e0,e1).normalize();
		  vecD<T> v(3); v[0] = n.x; v[1] = n.y; v[2] = n.z; return v; }
	T	getArea(const uvec& f) const						//! returns face area.
		{ const vec3<T> e0 = getPos(f[2])-getPos(f[0]), e1 = getPos(f[1])-getPos(f[0]);
		  T a = norm(cross(e0,e1).normalize()); if (f.size() == 3) return a*T(0.5);
		  const vec3<T> e2 = getPos(f[3])-getPos(f[1]), e3 = getPos(f[2])-getPos(f[1]);
		  a += norm(cross(e2,e3).normalize()); return a*T(0.5); }
	matD<T>	materialMatrix() const;
	matD<T>	stressMatrix() const;
	matD<T>	shapeMatrix() const;
	virtual void localShapeFunctionsAt(const vec3<T>& p) = 0;
	element(const unsigned int _nn, const unsigned int _ng, const vec3<T>* _gp,
		const vec3<T>* _nd, const cell& _c, const volumeMesh& _m, const unsigned int _dof,
		const std::vector<uvec>* _fc)						//! allocates element with specified properties.
		/*! Parameter nn corresponds to the number of nodes, ng to the number of Gauss points, c contains 
		    the element description from volume mesh m, with dof degrees of freedom per node.
		 */
		 : nn(_nn), ng(_ng), gp(_gp), nd(_nd), weight(ng), sf(nn),
		du(nn), dv(nn), dw(nn), dx(nn), dy(nn), dz(nn), 
		c(_c), m(_m), dof(_dof), det(T(0)), dp(nn), fc(_fc) { }
public:
	element(const element& e)							//! copies element properties from model e.
		 : nn(e.nn), ng(e.ng), gp(e.gp), nd(e.nd), weight(e.weight), sf(e.sf),
		du(e.du), dv(e.dv), dw(e.dw), dx(e.dx), dy(e.dy), dz(e.dz),
		c(e.c), m(e.m), dof(e.dof), det(e.det), dp(e.dp), fc(e.fc) { }
	element& operator=(const element& e)						//! assigns element properties from model e.
		{ nn = e.nn; ng = e.ng; gp = e.gp; nd = e.nd; weight = e.weight; sf = e.sf;
		  du = e.du; dv = e.dv; dw = e.dw; dx = e.dx; dy = e.dy; dz = e.dz;
		  c = e.c; m = e.m; dof = e.dof; det = e.det; dp = e.dp; fc = e.fc;
		  return *this; }
	virtual ~element() { }
	matD<T>	stiffnessMatrix();
	matD<T>	massMatrix();
	void	collectNodalStresses(const vecD<T>& u, matD<T>& s);
	void	setDisplacement(std::vector<vec3<T>>& _dp)
		{ dp = _dp; };
	void	globalShapeFunctions(const vec3<T>& p);
	matD<T>	viscousMatrix();
	void	collectNodalForces(const vecD<T>& u, vecD<T>& f);
	matD<T>	inertiaMatrix(const std::vector<vec3<T>>& v);
};

template <typename T>
matD<T> element<T>::stressMatrix() const
//! setup local stress/strain matrix.
{	matD<T> B(dof*2,dof*nn); B = T(0);
	for allNodes(i) {
		B(0,dof*i+0) = B(3,dof*i+1) = B(5,dof*i+2) = dx(i);			// Eq. 18.16
		B(1,dof*i+1) = B(3,dof*i+0) = B(4,dof*i+2) = dy(i);
		B(2,dof*i+2) = B(4,dof*i+1) = B(5,dof*i+0) = dz(i); };
	return B;
}

template <typename T>
matD<T> element<T>::shapeMatrix() const
//! setup local shape function matrix.
{	matD<T> H(dof,dof*nn); H = T(0);
	for allNodes(i) { for allDOF(d) H(d,dof*i+d) = sf(i); };
	return H;
}

template <typename T>
matD<T> element<T>::materialMatrix() const
//! setup isotropic material matrix.
{	const T nu = c.getPossionRatio(), E = c.getYoungModulus();
	if (nu < T(0) || nu > T(0.5)) throw optException("Material has invalid elasticity coefficient");
	const T fac = E/((T(1)+nu)*(T(1)-T(2)*nu));
	matD<T> C(dof*2,dof*2); C = 0; 							// see Eq. 16.35
	C(0,0) = C(1,1) = C(2,2) = T(1)-nu; C(3,3) = C(4,4) = C(5,5) = T(0.5)-nu;
	C(0,1) = C(1,0) = C(2,0) = C(0,2) = C(2,1) = C(1,2) = nu;
	C *= fac; return C;
}

template <typename T>
void element<T>::globalShapeFunctions(const vec3<T>& p)
//! compute global derivatives dx,dy,dz at point p.
{	localShapeFunctionsAt(p);							// set shape functions and derivatives
	vec3<T> u = T(0), v = T(0), w = T(0); 						// u,v,w are rows of Jacobian
	for allNodes(i) { vec3<T> q = getPos(i)+dp[i];					// position (from mesh) and update
		u += du[i]*q; v += dv[i]*q; w += dw[i]*q; };
    	const vec3<T> j0 = cross(v,w), j1 = cross(w,u), j2 = cross(u,v);		// columns of the inverse Jacobian
	det = u.x*j0.x+v.x*j1.x+w.x*j2.x;						// determinant
	if (det <= T(0)) throw rtException("Element has negative determinant");
	for allNodes(i) { u = (du[i]*j0+dv[i]*j1+dw[i]*j2)/det; 			// set global shape function derivatives
		dx[i] = u.x; dy[i] = u.y; dz[i] = u.z; };
}

template <typename T>
matD<T> element<T>::stiffnessMatrix()
//! assemble local element stiffness matrix K. Optimized for isotropic elastic behaviour.
{	const T nu = c.getPossionRatio(), E = c.getYoungModulus();
	if (nu < T(0) || nu > T(1)) throw optException("Material has invalid elasticity coefficient");
	const T fac = E/((T(1)+nu)*(T(1)-T(2)*nu));
	matD<T> C(2*dof,dof), K(dof*nn, dof*nn); K = T(0);
	assert(dof == 3);
	for allGPoints(k) { globalShapeFunctions(gp[k]); const T wk = det*weight(k);	// volume at k (in m^3)
		const T e0 = fac*(T(1)-nu)*wk, e1 = fac*(T(0.5)-nu)*wk, e2 = fac*nu*wk;
		for allNodes(j) { 							// AFEM Eq. 15.12
			C(0,0) = dx(j)*e0; C(0,1) = dy(j)*e2; C(0,2) = dz(j)*e2;
			C(1,0) = dx(j)*e2; C(1,1) = dy(j)*e0; C(1,2) = dz(j)*e2;
			C(2,0) = dx(j)*e2; C(2,1) = dy(j)*e2; C(2,2) = dz(j)*e0;
			C(3,0) = dy(j)*e1; C(3,1) = dx(j)*e1; C(3,2) = T(0);
			C(4,0) = T(0);     C(4,1) = dz(j)*e1; C(4,2) = dy(j)*e1;
			C(5,0) = dz(j)*e1; C(5,1) = T(0);     C(5,2) = dx(j)*e1;
			for allNodes(i) {
				K(i*dof+0,j*dof+0) += dx(i)*C(0,0)+dy(i)*C(3,0)+dz(i)*C(5,0);
				K(i*dof+0,j*dof+1) += dx(i)*C(0,1)+dy(i)*C(3,1)+dz(i)*C(5,1);
				K(i*dof+0,j*dof+2) += dx(i)*C(0,2)+dy(i)*C(3,2)+dz(i)*C(5,2);	
				K(i*dof+1,j*dof+0) += dy(i)*C(1,0)+dx(i)*C(3,0)+dz(i)*C(4,0);
				K(i*dof+1,j*dof+1) += dy(i)*C(1,1)+dx(i)*C(3,1)+dz(i)*C(4,1);
				K(i*dof+1,j*dof+2) += dy(i)*C(1,2)+dx(i)*C(3,2)+dz(i)*C(4,2);	
				K(i*dof+2,j*dof+0) += dz(i)*C(2,0)+dy(i)*C(4,0)+dx(i)*C(5,0);
				K(i*dof+2,j*dof+1) += dz(i)*C(2,1)+dy(i)*C(4,1)+dx(i)*C(5,1);
				K(i*dof+2,j*dof+2) += dz(i)*C(2,2)+dy(i)*C(4,2)+dx(i)*C(5,2); } } };	
	return K;
}

template <typename T>
matD<T> element<T>::massMatrix()
//! computes the local mass matrix.
{	const T rho = c.getDensity(); matD<T> M(dof*nn,dof*nn); M = T(0);
	for allGPoints(k) { globalShapeFunctions(gp[k]);
		const T wk = rho*det*weight(k);						// weight at k (in kg)	
		for allNodes(i) { const T wik = sf(i)*wk;
			for allNodes(j) { const T m = wik*sf(j);			// M has unit (kg), to be multiplied by acceleration (m/s^2)
				M(i*dof+0,j*dof+0) += m;
				M(i*dof+1,j*dof+1) += m;
				M(i*dof+2,j*dof+2) += m; } } };				// NB: eq 3 is 0 in a fluid model (pressure)
	return M;
}

template <typename T>
matD<T> element<T>::viscousMatrix()
//! computes the local dissipative viscous and pressure gradient terms for a fluid problem (see Yang, Eq. 51)
// this routine is either elegant or efficient but not both.
{	matD<T> K11(nn,nn), K12(nn,nn), K13(nn,nn), K22(nn,nn), K23(nn,nn), K33(nn,nn);
	matD<T> L1(nn,nn), L2(nn,nn), L3(nn,nn), L4(nn,nn);
	K11 = T(0); K12 = T(0); K13 = T(0); K22 = T(0); K23 = T(0); K33 = T(0);
	L1 = T(0); L2 = T(0); L3 = T(0), L4 = T(0);
	const T nu = c.getViscosity(), tau{1e3};
	for allGPoints(k) { globalShapeFunctions(gp[k]); const T wk = det*weight(k);	// volume at k (in m^3)
		for allNodes(i) { for allNodes(j) {
			K11(i,j) += nu*dx(i)*dx(j)*wk; K22(i,j) += nu*dy(i)*dy(j)*wk;	// see Eq. 52-57
			K33(i,j) += nu*dz(i)*dz(j)*wk; K12(i,j) += nu*dx(i)*dy(j)*wk;
			K13(i,j) += nu*dx(i)*dz(j)*wk; K23(i,j) += nu*dy(i)*dz(j)*wk;	// K has unit (kg/s), to be multiplied by velocity (m/s)
			L1(i,j) -= dx(i)*sf(j)*wk; L2(i,j) -= dy(i)*sf(j)*wk;		// L has unit (m^2), to be multiplied by pressure (kg/m s)
			L3(i,j) -= dz(i)*sf(j)*wk; };
			L4(i,i) -= sf(i)*sf(i)*wk/tau; } };
	const matD<T> Kd = K11+K22+K33;
	assert(dof == 4); matD<T> K(dof*nn,dof*nn); K = T(0);
	for allNodes(i) for allNodes(j) {
		K(i*dof+0,j*dof+0) = Kd(i,j)+K11(i,j); K(i*dof+0,j*dof+1) = K12(i,j);	// 1st row: K11,K12,K13,L1
		K(i*dof+0,j*dof+2) = K13(i,j); K(i*dof+0,j*dof+3) = L1(i,j);
		K(i*dof+1,j*dof+0) = K12(j,i); K(i*dof+1,j*dof+1) = Kd(i,j)+K22(i,j);	// 2nd row: K12^T,K22,K13,L2
		K(i*dof+1,j*dof+2) = K23(i,j); K(i*dof+1,j*dof+3) = L2(i,j);
		K(i*dof+2,j*dof+0) = K13(j,i); K(i*dof+2,j*dof+1) = K23(j,i);		// 3rd row: K13^T,K23^T,K33,L3
		K(i*dof+2,j*dof+2) = Kd(i,j)+K33(i,j); K(i*dof+2,j*dof+3) = L3(i,j);
		K(i*dof+3,j*dof+0) = L1(j,i); K(i*dof+3,j*dof+1) = L2(j,i);		// 4th row: L1^T,L3^T,L3^T,L4
		K(i*dof+3,j*dof+2) = L3(j,i); K(i*dof+3,j*dof+3) = L4(i,j); };
	return K;
}

template <typename T>
matD<T> element<T>::inertiaMatrix(const std::vector<vec3<T>>& v)
//! computes the local inertia matrix for a fluid problem (see Yang, Eq 69,71-73)
{	const T rho = c.getDensity();
	assert(dof == 4); matD<T> N(dof*nn,dof*nn); N = T(0);
	for allGPoints(k) { globalShapeFunctions(gp[k]);
		const T wk = rho*det*weight(k);						// weight at k (in kg)
		for allNodes(i) { const T wik = sf(i)*wk;
			for allNodes(j) {
				const T aj = dot(v[j],vec3<T>(dx(j),dy(j),dz(j)));	// velocity in spatial directions (1/s)
				N(i*dof+0,j*dof+0) += wik*aj;				// N has unit (kg/s), to be multiplied by velocity (m/s)
				N(i*dof+1,j*dof+1) += wik*aj;
				N(i*dof+2,j*dof+2) += wik*aj; } } };			// NB: eq 3 is 0 (pressure)
	return N;
}

template <typename T>
void element<T>::collectNodalStresses(const vecD<T>& u, matD<T>& S)
//! compute stress S at nodes, given displacements &dp.
{	vecD<T> ul(dof*nn); const matD<T> C = materialMatrix();
	for allNodes(i) { for allDOF(d) ul[i*dof+d] = u(c.vtx(i)*dof+d); };		// get displacements for this element
	for allNodes(i) { globalShapeFunctions(nd[i]);					// set dx, dy, dz & det
		const vecD<T> tau = (C*stressMatrix())*ul;				// stress at Gauss point, dof*2
		for (unsigned int j = 0; j < tau.N; j++) S(c.vtx(i),j) += tau(j); };	// return results in node-wise rows of S
}

template <typename T>
void element<T>::collectNodalForces(const vecD<T>& u, vecD<T>& f)
//! given a element-local vector of velocities and pressures u, computes contribution to nodal forces in f.
{	const T nu = c.getViscosity(); matD<T> G(3,3);
	for (const auto& k : *fc) { if (k.size() < 3) continue;				// for each face
		const vecD<T> n = getNormal(k);	const T a = getArea(k);			// get face normal and area
		for (const auto i : k) { const unsigned int id = c.vtx(i);		// for all nodes of this face
			globalShapeFunctions(nd[i]);
			const vec3<T> v(u(i*dof+0),u(i*dof+1),u(i*dof+2));		// velocity from fluid solution u
			G(0,0) = dx(i)*v.x*T(2); G(0,1) = dx(i)*v.y+dy(i)*v.x; G(1,0) = G(0,1);
			G(1,1) = dy(i)*v.y*T(2); G(0,2) = dx(i)*v.z+dz(i)*v.x; G(2,0) = G(0,2);
			G(2,2) = dz(i)*v.z*T(2); G(1,2) = dy(i)*v.z+dz(i)*v.y; G(2,1) = G(1,2);
			const vecD<T> fi = G*((nu*a)*n)-(u(i*dof+3)*a)*n; 		// compute force at node i in face k
			for (unsigned int d = 0; d < 3; d++) f[id*3+d] += fi(d); } };	// update force in global vector
}

//! Implements a tetrahedral element with four nodes.

template <typename T> class tetElement4 : public element<T> {
	static const vec3<T> nds[];			//!< local element nodes
	static const vec3<T> gps[];			//!< Gauss points
	static const std::vector<uvec> faces;

	void	localShapeFunctionsAt(const vec3<T> &p);
public:
	tetElement4<T>(const cell& c, const volumeMesh& m, const unsigned int dof)	//! allocates a tetrahedral element with four nodes.
		 : element<T>(4, 1, gps, nds, c, m, dof, &faces)
		{ element<T>::weight = T(1.0/6.0); }
};

template <typename T>
const std::vector<uvec> tetElement4<T>::faces = {
	{ 0, 1, 2 }, { 0, 3, 1 }, { 0, 2, 3 }, { 1, 3, 2} };

template <typename T>
const vec3<T> tetElement4<T>::nds[]  = { vec3<T>(0,0,0), vec3<T>(1,0,0), vec3<T>(0,1,0), vec3<T>(0,0,1) }; 

template <typename T>
const vec3<T> tetElement4<T>::gps[]  = { vec3<T>(0.25,0.25,0.25) }; 

template <typename T>
void tetElement4<T>::localShapeFunctionsAt(const vec3<T>& p) 				//! interpolate shape functions at p.
{	// derivatives are independent of p
	element<T>::sf[0] = T(1.0-p.x-p.y-p.z);
	element<T>::sf[1] = p.x; element<T>::sf[2] = p.y; element<T>::sf[3] = p.z;
	element<T>::du = T(0); element<T>::du[0] = T(-1); element<T>::du[1] = T(1);
	element<T>::dv = T(0); element<T>::dv[0] = T(-1); element<T>::dv[2] = T(1);
	element<T>::dw = T(0); element<T>::dw[0] = T(-1); element<T>::dw[3] = T(1);
}

//! Implements a tetrahedral element with ten nodes.

template <typename T> class tetElement10 : public element<T> {
	static constexpr T t10a = T((5.0+3.0*std::sqrt(5.0))/20.0);
	static constexpr T t10b = T((5.0-std::sqrt(5.0))/20.0);
	static const vec3<T> nds[];			//!< local element nodes
	static const vec3<T> gps[];			//!< Gauss points
	static const std::vector<uvec> faces;

	void	localShapeFunctionsAt(const vec3<T> &p);
public:
	tetElement10<T>(const cell& c, const volumeMesh& m, const unsigned int dof)	//! allocates a tetrahedral element with ten nodes.
		 : element<T>(10, 4, gps, nds, c, m, dof, &faces)
		{ element<T>::weight = T(1.0/6.0); }
};

template <typename T>
const std::vector<uvec> tetElement10<T>::faces = {
	{ 0, 1, 2 }, { 0, 3, 1 }, { 0, 2, 3 }, { 1, 3, 2} };

template <typename T>
const vec3<T> tetElement10<T>::nds[] = {
	vec3<T>(0,0,0), vec3<T>(1,0,0), vec3<T>(0,1,0), vec3<T>(0,0,1), vec3<T>(0.5,0,0),
	vec3<T>(0.5,0.5,0), vec3<T>(0,0.5,0), vec3<T>(0,0,0.5), vec3<T>(0.5,0,0.5), vec3<T>(0,0.5,0.5) }; 

template <typename T>
const vec3<T> tetElement10<T>::gps[] = {
	vec3<T>(t10a,t10b,t10b), vec3<T>(t10b,t10a,t10b), vec3<T>(t10b,t10b,t10a), vec3<T>(t10b,t10b,t10b) };

template <typename T>
void tetElement10<T>::localShapeFunctionsAt(const vec3<T>& p) 				//! interpolate shape functions at p.
{	T t = 1-p.x-p.y-p.z, u = p.x, v = p.y, w = p.z;
	element<T>::sf[0] = t*(2*t-1); element<T>::sf[1] = u*(2*u-1);			// shape functions at corners
	element<T>::sf[2] = v*(2*v-1); element<T>::sf[3] = w*(2*w-1);
	element<T>::sf[4] = 4*t*u; element<T>::sf[5] = 4*u*v; element<T>::sf[6] = 4*v*t; // shape functions at mid-edges
	element<T>::sf[7] = 4*t*w; element<T>::sf[8] = 4*u*w; element<T>::sf[9] = 4*v*w;
	element<T>::du = 0; element<T>::du[0] = 1-4*t;					// derivatives
	element<T>::du[1] = 4*u-1; element<T>::du[4] = 4*(t-u);
	element<T>::du[5] = 4*v; element<T>::du[6] = -4*v;
	element<T>::du[7] = -4*w; element<T>::du[8] = 4*w;
	element<T>::dv = 0; element<T>::dv[0] = element<T>::du[0];
	element<T>::dv[2] = 4*v-1; element<T>::dv[4] = -4*u; element<T>::dv[5] = 4*u;
	element<T>::dv[6] = 4*(t-v); element<T>::dv[7] = -4*w; element<T>::dv[9] = 4*w;
	element<T>::dw = 0; element<T>::dw[0] = element<T>::du[0];
	element<T>::dw[3] = 4*w-1; element<T>::dw[4] = -4*u; element<T>::dw[6] = -4*v; 
	element<T>::dw[7] = 4*(t-w); element<T>::dw[8] = 4*u; element<T>::dw[9] = 4*v;
};
//! Implements a hexahedral element with eight nodes.

template <typename T> class hexElement8 : public element<T> {
	static constexpr T h8 = T(std::sqrt(1.0/3.0));
	static const vec3<T> nds[];			//!< local element nodes
	static const vec3<T> gps[];			//!< Gauss points
	static const std::vector<uvec> faces;

	void	localShapeFunctionsAt(const vec3<T> &p);
public:
	hexElement8<T>(const cell& c, const volumeMesh& m, const unsigned int dof)	//! allocates a hexahedral element with eight nodes.
		 : element<T>(8, 8, gps, nds, c, m, dof, &faces)
		{ element<T>::weight = T(1.0); }
}; 

template <typename T>
const std::vector<uvec> hexElement8<T>::faces = {
	{ 0, 1, 2, 3 }, { 0, 4, 5, 1 }, { 0, 3, 7, 4 },
	{ 4, 7, 6, 5 }, { 2, 6, 7, 3 },	{ 1, 5, 6, 2 } };

template <typename T>
const vec3<T> hexElement8<T>::nds[]  = {
	vec3<T>(-1,-1,-1), vec3<T>( 1,-1,-1), vec3<T>( 1, 1,-1), vec3<T>(-1, 1,-1), 
        vec3<T>(-1,-1, 1), vec3<T>( 1,-1, 1), vec3<T>( 1, 1, 1), vec3<T>(-1, 1, 1) };

template <typename T>
const vec3<T> hexElement8<T>::gps[]  = {
	vec3<T>(-h8,-h8,-h8), vec3<T>( h8,-h8,-h8), vec3<T>( h8, h8,-h8), vec3<T>(-h8, h8,-h8), 
        vec3<T>(-h8,-h8, h8), vec3<T>( h8,-h8, h8), vec3<T>( h8, h8, h8), vec3<T>(-h8, h8, h8) };

template <typename T>
void hexElement8<T>::localShapeFunctionsAt(const vec3<T>& p) 				//! interpolate shape functions at p.
{	for allNodes(i) { T a = element<T>::nd[i].x*p.x, b = element<T>::nd[i].y*p.y, c = element<T>::nd[i].z*p.z;
		element<T>::sf[i] = T(0.125)*(T(1)+a)*(T(1)+b)*(T(1)+c);
		element<T>::du[i] = T(0.125)*element<T>::nd[i].x*(T(1)+b+c+b*c);
		element<T>::dv[i] = T(0.125)*element<T>::nd[i].y*(T(1)+a+c+a*c);
		element<T>::dw[i] = T(0.125)*element<T>::nd[i].z*(T(1)+a+b+a*b); }
} 

//! Implements a hexahedral element with twenty nodes.

template <typename T> class hexElement20 : public element<T> {
	static constexpr T w20a = T(8.0/9.0);
	static constexpr T w20b = T(5.0/9.0);
	static constexpr T f20 = T(5.0);
	static constexpr T h20 = T(std::sqrt(0.6));
	static const vec3<T> nds[];		//!< local element nodes
	static const vec3<T> gps[];		//!< Gauss points
	static const std::vector<uvec> faces;

	void	localShapeFunctionsAt(const vec3<T> &p);
public:
	hexElement20<T>(const cell& c, const volumeMesh& m, const unsigned int dof)	//! allocates a hexahedral element with twenty nodes.
		 : element<T>(20, 27, gps, nds, c, m, dof, &faces)
		{ element<T>::weight[0] = w20a*w20a*w20a*f20;
		  for (unsigned int i = 1; i < 7; i++) element<T>::weight[i] = w20a*w20a*w20b*f20; 
		  for (unsigned int i = 7; i < 19; i++) element<T>::weight[i] = w20a*w20b*w20b*f20; 
		  for (unsigned int i = 19; i < 27; i++) element<T>::weight[i] = w20b*w20b*w20b*f20; }
};

template <typename T>
const std::vector<uvec> hexElement20<T>::faces = {
	{ 0, 1, 2, 3 }, { 0, 4, 5, 1 }, { 0, 3, 7, 4 },
	{ 4, 7, 6, 5 }, { 2, 6, 7, 3 },	{ 1, 5, 6, 2 } };

template <typename T>
const vec3<T> hexElement20<T>::nds[] = {
	vec3<T>(-1,-1,-1), vec3<T>( 1,-1,-1), vec3<T>( 1, 1,-1), vec3<T>(-1, 1,-1), 
	vec3<T>(-1,-1, 1), vec3<T>( 1,-1, 1), vec3<T>( 1, 1, 1), vec3<T>(-1, 1, 1),
	vec3<T>( 0,-1,-1), vec3<T>( 1, 0,-1), vec3<T>( 0, 1,-1), vec3<T>(-1, 0,-1),
	vec3<T>( 0,-1, 1), vec3<T>( 1, 0, 1), vec3<T>( 0, 1, 1), vec3<T>(-1, 0, 1),
	vec3<T>(-1,-1, 0), vec3<T>( 1,-1, 0), vec3<T>( 1, 1, 0), vec3<T>(-1, 1, 0) }; 

template <typename T>
const vec3<T> hexElement20<T>::gps[] = {
	vec3<T>(0,0,0), vec3<T>(-h20,0,0), vec3<T>(h20,0,0), vec3<T>(0,-h20,0),
	vec3<T>(0,h20,0), vec3<T>(0,0,-h20), vec3<T>(0,0,h20), vec3<T>(-h20,-h20,0),
	vec3<T>(-h20,h20,0), vec3<T>(h20,-h20,0), vec3<T>(h20,h20,0), vec3<T>(-h20,0,-h20),
	vec3<T>(-h20,0,h20), vec3<T>(h20,0,-h20), vec3<T>(h20,0,h20), vec3<T>(0,-h20,-h20),
	vec3<T>(0,-h20,h20), vec3<T>(0,h20,-h20), vec3<T>(0,h20,h20), vec3<T>(-h20,-h20,-h20),
	vec3<T>(-h20,h20,-h20), vec3<T>(h20,-h20,-h20), vec3<T>(h20,h20,-h20), vec3<T>(-h20,-h20,h20),
	vec3<T>(-h20,h20,h20), vec3<T>(h20,-h20,h20), vec3<T>(h20,h20,h20) };

template <typename T>
void hexElement20<T>::localShapeFunctionsAt(const vec3<T>& p) 				//! interpolate shape functions at p.
{	T u = p.x, um = T(1)-u, up = T(1)+u, u2 = T(1)-u*u;
	T v = p.y, vm = T(1)-v, vp = T(1)+v, v2 = T(1)-v*v;
	T w = p.z, wm = T(1)-w, wp = T(1)+w, w2 = T(1)-w*w;
	element<T>::sf[0] = T(0.125)*um*vm*wm*(-u-v-w-T(2));				// shape functions at corners
	element<T>::sf[1] = T(0.125)*up*vm*wm*( u-v-w-T(2));
	element<T>::sf[2] = T(0.125)*up*vp*wm*( u+v-w-T(2));
	element<T>::sf[3] = T(0.125)*um*vp*wm*(-u+v-w-T(2));
	element<T>::sf[4] = T(0.125)*um*vm*wp*(-u-v+w-T(2));
	element<T>::sf[5] = T(0.125)*up*vm*wp*( u-v+w-T(2));
	element<T>::sf[6] = T(0.125)*up*vp*wp*( u+v+w-T(2));
	element<T>::sf[7] = T(0.125)*um*vp*wp*(-u+v+w-T(2));
	element<T>::sf[8] =  T(0.25)*u2*vm*wm; element<T>::sf[9] =  T(0.25)*v2*up*wm;  	// shape functions at mid-edges
	element<T>::sf[10] = T(0.25)*u2*vp*wm; element<T>::sf[11] = T(0.25)*v2*um*wm;
	element<T>::sf[12] = T(0.25)*u2*vm*wp; element<T>::sf[13] = T(0.25)*v2*up*wp;
	element<T>::sf[14] = T(0.25)*u2*vp*wp; element<T>::sf[15] = T(0.25)*v2*um*wp;
	element<T>::sf[16] = T(0.25)*w2*um*vm; element<T>::sf[17] = T(0.25)*w2*up*vm;
	element<T>::sf[18] = T(0.25)*w2*up*vp; element<T>::sf[19] = T(0.25)*w2*um*vp;
	element<T>::du[0] = T(-0.125)*(vm*wm-T(2)*u*vm*wm-v*vm*wm-w*vm*wm-T(2)*vm*wm);	// derivatives
	element<T>::du[1] =  T(0.125)*(vm*wm+T(2)*u*vm*wm-v*vm*wm-w*vm*wm-T(2)*vm*wm);
	element<T>::du[2] =  T(0.125)*(vp*wm+T(2)*u*vp*wm+v*vp*wm-w*vp*wm-T(2)*vp*wm);
	element<T>::du[3] = T(-0.125)*(vp*wm-T(2)*u*vp*wm+v*vp*wm-w*vp*wm-T(2)*vp*wm);
	element<T>::du[4] = T(-0.125)*(vm*wp-T(2)*u*vm*wp-v*vm*wp+w*vm*wp-T(2)*vm*wp);
	element<T>::du[5] =  T(0.125)*(vm*wp+T(2)*u*vm*wp-v*vm*wp+w*vm*wp-T(2)*vm*wp);
	element<T>::du[6] =  T(0.125)*(vp*wp+T(2)*u*vp*wp+v*vp*wp+w*vp*wp-T(2)*vp*wp);
	element<T>::du[7] = T(-0.125)*(vp*wp-T(2)*u*vp*wp+v*vp*wp+w*vp*wp-T(2)*vp*wp);
	element<T>::du[8] = T(-0.5)*u*vm*wm; element<T>::du[9] =  T(0.25)*(wm-v*v*wm);
	element<T>::du[10] = T(-0.5)*u*vp*wm; element<T>::du[11] = -T(0.25)*(wm-v*v*wm);
	element<T>::du[12] = T(-0.5)*u*vm*wp; element<T>::du[13] =  T(0.25)*(wp-v*v*wp);
	element<T>::du[14] = T(-0.5)*u*vp*wp; element<T>::du[15] = -T(0.25)*(wp-v*v*wp);
	element<T>::du[16] = -T(0.25)*(vm-w*w*vm); element<T>::du[17] =  T(0.25)*(vm-w*w*vm);
	element<T>::du[18] =  T(0.25)*(vp-w*w*vp); element<T>::du[19] = -T(0.25)*(vp-w*w*vp);
	element<T>::dv[0] = T(-0.125)*(um*wm-T(2)*v*um*wm-u*um*wm-w*um*wm-T(2)*um*wm);
	element<T>::dv[1] = T(-0.125)*(up*wm-T(2)*v*up*wm+u*up*wm-w*up*wm-T(2)*up*wm); 
	element<T>::dv[2] =  T(0.125)*(up*wm+T(2)*v*up*wm+u*up*wm-w*up*wm-T(2)*up*wm); 
	element<T>::dv[3] =  T(0.125)*(um*wm+T(2)*v*um*wm-u*um*wm-w*um*wm-T(2)*um*wm); 
	element<T>::dv[4] = T(-0.125)*(um*wp-T(2)*v*um*wp-u*um*wp+w*um*wp-T(2)*um*wp); 
	element<T>::dv[5] = T(-0.125)*(up*wp-T(2)*v*up*wp+u*up*wp+w*up*wp-T(2)*up*wp); 
	element<T>::dv[6] =  T(0.125)*(up*wp+T(2)*v*up*wp+u*up*wp+w*up*wp-T(2)*up*wp); 
	element<T>::dv[7] =  T(0.125)*(um*wp+T(2)*v*um*wp-u*um*wp+w*um*wp-T(2)*um*wp); 
	element<T>::dv[8] = -T(0.25)*(wm-u*u*wm); element<T>::dv[9] = T(-0.5)*v*up*wm;
	element<T>::dv[10] =  T(0.25)*(wm-u*u*wm); element<T>::dv[11] = T(-0.5)*v*um*wm;
	element<T>::dv[12] = -T(0.25)*(wp-u*u*wp); element<T>::dv[13] = T(-0.5)*v*up*wp;
	element<T>::dv[14] =  T(0.25)*(wp-u*u*wp); element<T>::dv[15] = T(-0.5)*v*um*wp;
	element<T>::dv[16] = -T(0.25)*(um-w*w*um); element<T>::dv[17] = -T(0.25)*(up-w*w*up); 
	element<T>::dv[18] =  T(0.25)*(up-w*w*up); element<T>::dv[19] =  T(0.25)*(um-w*w*um);
	element<T>::dw[0] = T(-0.125)*(um*vm-T(2)*w*um*vm-u*um*vm-v*um*vm-T(2)*um*vm); 
	element<T>::dw[1] = T(-0.125)*(up*vm-T(2)*w*up*vm+u*up*vm-v*up*vm-T(2)*up*vm);  
	element<T>::dw[2] = T(-0.125)*(up*vp-T(2)*w*up*vp+u*up*vp+v*up*vp-T(2)*up*vp);  
	element<T>::dw[3] = T(-0.125)*(um*vp-T(2)*w*um*vp-u*um*vp+v*um*vp-T(2)*um*vp);  
	element<T>::dw[4] =  T(0.125)*(um*vm+T(2)*w*um*vm-u*um*vm-v*um*vm-T(2)*um*vm);  
	element<T>::dw[5] =  T(0.125)*(up*vm+T(2)*w*up*vm+u*up*vm-v*up*vm-T(2)*up*vm);  
	element<T>::dw[6] =  T(0.125)*(up*vp+T(2)*w*up*vp+u*up*vp+v*up*vp-T(2)*up*vp);  
	element<T>::dw[7] =  T(0.125)*(um*vp+T(2)*w*um*vp-u*um*vp+v*um*vp-T(2)*um*vp);  
	element<T>::dw[8] = -T(0.25)*(vm-u*u*vm); element<T>::dw[9] = -T(0.25)*(up-v*v*up);
	element<T>::dw[10] = -T(0.25)*(vp-u*u*vp); element<T>::dw[11] = -T(0.25)*(um-v*v*um);
	element<T>::dw[12] =  T(0.25)*(vm-u*u*vm); element<T>::dw[13] =  T(0.25)*(up-v*v*up);
	element<T>::dw[14] =  T(0.25)*(vp-u*u*vp); element<T>::dw[15] =  T(0.25)*(um-v*v*um);
	element<T>::dw[16] = T(-0.5)*w*um*vm; element<T>::dw[17] = T(-0.5)*w*up*vm;
	element<T>::dw[18] = T(-0.5)*w*up*vp; element<T>::dw[19] = T(-0.5)*w*um*vp;
}

#endif

