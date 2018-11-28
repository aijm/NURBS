// This file is part of NURBS, a simple NURBS library.
// github repo: https://github.com/aijm/NURBS
// Copyright (C) 2018 Jiaming Ai <aichangeworld@gmail.com>

#include "BSplineCurve.h"

/*input format:
_n       : P_0,P_1,...,P_n; _n is the final index
_k       : order of BSpline
_controlP: P_0,P_1,...,P_n; (n+1) by 2 or 3
_knots   : t_0,t_1,...,t_(n+k); */
BSplineCurve::BSplineCurve(
	int _n, int _k, 
	MatrixXd _controlP, 
	vector<double> _knots)
{
	assert(_k >= 1 && _n >= _k - 1 && _controlP.rows() == _n + 1 && _knots.size() == _n + _k + 1);
	assert(_controlP.cols() == 2 || _controlP.cols() == 3);
	n = _n;
	k = _k;
	controlP = _controlP;
	knots = _knots;
}

int BSplineCurve::find_ind(double t)
{
	// find the knot interval of t by binary searching
	if (t == knots[n + 1]) return n;
	int low = k - 1;
	int high = n + 1;
	assert(t >= knots[low] && t < knots[high]);

	int mid = (low + high) / 2;
	while (t < knots[mid] || t >= knots[mid + 1])
	{
		if (t < knots[mid]) high = mid;
		else low = mid;
		mid = (low + high) / 2;
	}
	return mid;
}

// evaluate the coordinate of curvePoint with parameter t  
MatrixXd BSplineCurve::eval(double t)
{
	// find the knot interval of t by binary searching
	int L = find_ind(t); //[t_L,t_(L+1)]  

	// P_(L-k+1),..,P_L control the interval [t_L,t_(L+1)]
	MatrixXd temp = controlP.block(L - k + 1, 0, k, controlP.cols());

	// de-Boor algorithm
	for (int r = 1; r <= k - 1; r++)
		for (int i = L - k + 1 + r; i <= k - 1; i++)
		{
			double factor = (t - knots[i]) / (knots[i + k - r] - knots[i]);
			int start = i - (L - k + 1 + r);
			temp.row(start) = (1.0 - factor)*temp.row(start) + factor*temp.row(start + 1);
		}
	MatrixXd curvePoint = temp.row(0);
	return curvePoint;
}

// display by libigl
void BSplineCurve::show(igl::opengl::glfw::Viewer& viewer, double resolution)
{
	double left = knots[k - 1];
	double right = knots[n + 1];
	const int num = (right - left) / resolution;
	double new_resolution = (right - left) / num;
	// plot BSplineCurve
	for (int i = 0; i < num; i++)
	{
		Eigen::MatrixXd P1 = eval(left + 1.0*i* new_resolution);

		Eigen::MatrixXd P2 = eval(left + 1.0*(i + 1)*new_resolution);
		viewer.data().add_edges(P1, P2, Eigen::RowVector3d(1, 1, 1));
	}
	// plot control points
	viewer.data().add_points(controlP, Eigen::RowVector3d(1, 0, 0));
	// plot control polygon
	for (int i = 0; i < n; i++)
	{
		viewer.data().add_edges(
			controlP.row(i),
			controlP.row(i + 1),
			Eigen::RowVector3d(1, 0, 0));
	}
}


