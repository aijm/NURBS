// This file is part of NURBS, a simple NURBS library.
// github repo: https://github.com/aijm/NURBS
// Copyright (C) 2018 Jiaming Ai <aichangeworld@gmail.com>

#ifndef BSPLINECURVE_H
#define BSPLINECURVE_H

#include <igl/opengl/glfw/Viewer.h>
using namespace Eigen;
using namespace std;
struct BSplineCurve
{
	/*input format:
	_n       : P_0,P_1,...,P_n; _n is the final index
	_k       : order of BSpline
	_controlP: P_0,P_1,...,P_n; (n+1) by 2 or 3
	_knots   : t_0,t_1,...,t_(n+k); */
	BSplineCurve(int _n, int _k, MatrixXd _controlP, vector<double> _knots);

	// find the knot interval of t by binary searching
	int find_ind(double t);

	// evaluate the coordinate of curvePoint with parameter t  
	MatrixXd eval(double t);
	
	// display by libigl
	void show(igl::opengl::glfw::Viewer& viewer, double resolution = 0.01);
	

	int n; // P_0,P_1,...,P_n; _n is the final index 
	int k; // order of BSpline
	vector<double> knots;
	MatrixXd controlP;
};


#endif // !BSPLINECURVE_H


