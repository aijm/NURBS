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
	/*input format:
	_n       : P_0,P_1,...,P_n; _n is the final index
	_k       : order of BSpline
	_controlP: P_0,P_1,...,P_n; (n+1) by 2 or 3
	_knots   : t_0,t_1,...,t_(n+k); */
	NURBSCurve(int _n, int _k, MatrixXd _controlP, VectorXd _knots, bool _isRational = false);

	// find the knot interval of t by binary searching
	int find_ind(double t);

	// evaluate the coordinate of curvePoint with parameter t  
	MatrixXd eval(double t);

	// display by libigl
	void show(igl::opengl::glfw::Viewer& viewer, double resolution = 0.01);

	bool isRational;
	int n; // P_0,P_1,...,P_n; _n is the final index
	int k; // order of BSpline
	VectorXd knots;
	MatrixXd controlP;
	MatrixXd controlPw;

};
#endif // !NURBSCURVE_H




