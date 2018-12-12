// This file is part of NURBS, a simple NURBS library.
// github repo: https://github.com/aijm/NURBS
// Copyright (C) 2018 Jiaming Ai <aichangeworld@gmail.com>

#ifndef NURBSSURFACE_H
#define NURBSSURFACE_H

#include <igl/opengl/glfw/Viewer.h>
#include "NURBSCurve.h"
using namespace Eigen;
using namespace std;
struct NURBSSurface
{
	NURBSSurface(){}
	/*input format:
	_order       : // _order(0):u direction; _order(1): v direction
	_controlP: _controlP[i] represents u direction control point ,matrix (m+1) by 3 or 4
				  v
				  |
	_controlP[n]: | P_0n P_1n ... P_mn
	_controlP[i]: | ...
	_controlP[1]: | P_01 P_11 ... P_m1
	_controlP[0]: | P_00 P_10 ... P_m0
				   ------------------------> u

	_uknots   : u_0,u_1,...,u_(m+u_order)
	_vknots   : v_0,v_1,...,v_(n+v_order)
	_isRational:          */

	NURBSSurface(VectorXi _order, vector<MatrixXd> _controlP, VectorXd _uknots, VectorXd _vknots, bool _isRational = false);

	// load
	bool loadNURBS(string);

	// save 
	bool saveNURBS(string);

	// find the knot interval of t by binary searching
	int find_ind(double t, int k, int n, const VectorXd& knots);

	// calculate coordinate of curve point with parameter u & v
	MatrixXd eval(double u, double v);

	MatrixXd eval(double t, const MatrixXd &_controlP, const VectorXd &knots);

	// knot insertion
	bool insert(double s, char dir='u');
	// kont insertion
	bool insert(double s, double t);

	// display by libigl
	void draw(igl::opengl::glfw::Viewer& viewer, bool showpolygon=true,bool showsurface=true,double resolution = 0.01);

	// draw controlpolygon
	void drawControlPolygon(igl::opengl::glfw::Viewer &viewer);

	// draw NURBS surface
	void drawSurface(igl::opengl::glfw::Viewer &viewer, double resolution = 0.01);

	// surface skinning
	void skinning(const vector<NURBSCurve> &curves, igl::opengl::glfw::Viewer &viewer);

	bool isRational = false;
	int u_order; // order of u direction
	int v_order; // order of v direction
	int u_num; // the final index of u direction control point
	int v_num; // the final index of v direction control point
	VectorXd uknots; // knots vector u direction : u_0, u_1, ..., u_(m + u_order)
	VectorXd vknots; // knots vector v direction : v_0,v_1,...,v_(n+v_order)
	vector<MatrixXd> controlP;
	vector<MatrixXd> controlPw;
	int dimension; // the dimension of control point 2 or 3 or 4

	// use libigl mesh(V,F) structure to show the surface
	MatrixXd mesh_V;
	MatrixXi mesh_F;
};


#endif // !NURBSSURFACE_H




