// This file is part of NURBS, a simple NURBS library.
// github repo: https://github.com/aijm/NURBS
// Copyright (C) 2018 Jiaming Ai <aichangeworld@gmail.com>


#ifndef NURBSCURVE_H
#define NURBSCURVE_H

#include <igl/opengl/glfw/Viewer.h>
using namespace Eigen;
using namespace std;
struct NURBSCurve
{
	NURBSCurve(){}
	/*input format:
	_n       : P_0,P_1,...,P_n; _n is the final index
	_k       : order of BSpline
	_controlP: P_0,P_1,...,P_n; (n+1) by 2 or 3
	_knots   : t_0,t_1,...,t_(n+k); */
	NURBSCurve(int _n, int _k, MatrixXd _controlP, VectorXd _knots, bool _isRational = false);

	NURBSCurve(const NURBSCurve &curve){
		isRational = curve.isRational;
		n = curve.n;
		k = curve.k;
		knots = curve.knots;
		controlPw = curve.controlPw;
	}

	NURBSCurve& operator=(const NURBSCurve &curve){
		isRational = curve.isRational;
		n = curve.n;
		k = curve.k;
		knots = curve.knots;
		controlPw = curve.controlPw;
		return *this;
	}
	// load
	bool loadNURBS(string);
	// save
	bool saveNURBS(string);
	
	// find the knot interval of t by binary searching
	int find_ind(double t)const;

	// evaluate the coordinate of curvePoint with parameter t  
	MatrixXd eval(double t)const;

	// chord length parameterization
	static VectorXd parameterize(const MatrixXd &points);
	// basis function N_(i,p)(t)
	double basis(int i, int p, double t, const VectorXd &knotvector);

	// interpolate by bspline of degree 3
	void interpolate(const MatrixXd &points);

	// interpolate with appointed knot vector
	void interpolate(const MatrixXd &points, const VectorXd &knotvector);

	// pia fit by B-spline of degree 3
	void piafit(const MatrixXd &points,int max_iter_num=100, double eps=1e-5);

	// pia fit with appointed knot vector
	void piafit(const MatrixXd &points, const VectorXd &knotvector, int max_iter_num = 100, double eps = 1e-5);

	// given Q_0,...,Q_m, fit by B-spline with control points P_0,...,P_n
	void lspiafit(const MatrixXd & points, const int &n_cpts);

	// lspia fit with appointed knot vector
	void lspiafit(const MatrixXd & points, const int &n_cpts, const VectorXd & knotvector);

	// kont insertion
	bool insert(double t);

	// display by libigl
	void draw(igl::opengl::glfw::Viewer& viewer, bool showpolygon=true,bool showsurface=true,double resolution = 0.01);

	// draw controlpolygon
	void drawControlPolygon(igl::opengl::glfw::Viewer &viewer);

	// draw NURBS surface
	void drawSurface(igl::opengl::glfw::Viewer &viewer, double resolution = 0.01);


	bool isRational = false;
	int n; // P_0,P_1,...,P_n; _n is the final index
	int k; // order of BSpline
	VectorXd knots;
	MatrixXd controlP;
	MatrixXd controlPw;

};
#endif // !NURBSCURVE_H




