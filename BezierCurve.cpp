// This file is part of NURBS, a simple NURBS library.
// github repo: https://github.com/aijm/NURBS
// Copyright (C) 2018 Jiaming Ai <aichangeworld@gmail.com>

#include "BezierCurve.h"

BezierCurve::BezierCurve(Eigen::MatrixXd P)
{
	assert(P.rows() > 1 && (P.cols() == 2 || P.cols() == 3));
	controlP = P;
}

// evaluate the coordinate of curvePoint with parameter t  
Eigen::MatrixXd BezierCurve::eval(double t)
{
	Eigen::MatrixXd temp = controlP;
	int n = temp.rows() - 1; // P_0,...,P_n
	for (int k = 1; k <= n; k++)
		for (int i = 0; i <= n - k; i++)
			temp.row(i) = (1 - t)*temp.row(i) + t*temp.row(i + 1);
	return temp.row(0);
}

// display by libigl
void BezierCurve::show(igl::opengl::glfw::Viewer& viewer, double resolution)
{
	const int n = 1 / resolution;
	// plot BezierCurve
	for (int i = 0; i < n; i++)
	{
		Eigen::MatrixXd P1 = eval(1.0*i / n);

		Eigen::MatrixXd P2 = eval(1.0*(i + 1) / n);
		viewer.data().add_edges(P1, P2, Eigen::RowVector3d(1, 1, 1));
	}
	// plot control points
	viewer.data().add_points(controlP, Eigen::RowVector3d(1, 0, 0));
	// plot control polygon
	for (int i = 0; i < controlP.rows() - 1; i++)
	{
		viewer.data().add_edges(
			controlP.row(i),
			controlP.row(i + 1),
			Eigen::RowVector3d(1, 0, 0));
	}
}







