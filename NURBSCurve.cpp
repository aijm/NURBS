// This file is part of NURBS, a simple NURBS library.
// github repo: https://github.com/aijm/NURBS
// Copyright (C) 2018 Jiaming Ai <aichangeworld@gmail.com>

#include "NURBSCurve.h"

/*input format:
_n       : P_0,P_1,...,P_n; _n is the final index 
_k       : order of BSpline
_controlP: P_0,P_1,...,P_n; (n+1) by 2 or 3
_knots   : t_0,t_1,...,t_(n+k); */
NURBSCurve::NURBSCurve(int _n, int _k, MatrixXd _controlP, VectorXd _knots, bool _isRational)
{
	assert(_k >= 1 && _n >= _k - 1 && _controlP.rows() == _n + 1 && _knots.size() == _n + _k + 1);

	n = _n;
	k = _k;
	knots = _knots;
	controlPw = _controlP;
	isRational = _isRational;

	// if (_isRational) {
	// 	assert(_controlP.cols() == 3 || _controlP.cols() == 4);
	// 	controlP = controlPw.rowwise().hnormalized();
	// }
	// else {
	// 	assert(_controlP.cols() == 2 || _controlP.cols() == 3);
	// 	controlP = controlPw;
	// }

}

// load
bool NURBSCurve::loadNURBS(string name){
	ifstream in(name.c_str());
	if(!in){
		return false;
	}
	in >> isRational;
	in >> n >> k;
	int dimension=3;
	in >> dimension;

	controlPw = MatrixXd(n+1,dimension);
	knots = VectorXd(n+k+1);
	// to do.
	for(int i=0;i < controlPw.rows();i++){
		for(int j=0;j< controlPw.cols();j++){
			in >> controlPw(i,j);
		}
	}
	for(int i=0;i<knots.size();i++){
		in >> knots(i);
	}
	return true;
}
// save
bool NURBSCurve::saveNURBS(string name){
	if(isRational){
		name += ".cptw";
	}else{
		name += ".cpt";
	}
	ofstream out(name.c_str());
	if(!out){
		return false;
	}
	out<< isRational<<endl;
	out<< n <<" "<< k<<endl;
	out<< controlPw.cols()<<endl;
	out<< controlPw <<endl;
	out<< knots.transpose();
	return true;
}

// find the knot interval of t by binary searching
int NURBSCurve::find_ind(double t)
{
	
	if (t == knots(n + 1)) return n;
	int low = 0;
	int high = n+k;
	assert(t >= knots(low) && t < knots(high));

	int mid = (low + high) / 2;
	while (t < knots(mid) || t >= knots(mid + 1))
	{
		if (t < knots(mid)) high = mid;
		else low = mid;
		mid = (low + high) / 2;
	}
	return mid;
}

// evaluate the coordinate of curvePoint with parameter t  
MatrixXd NURBSCurve::eval(double t)
{
	// find the knot interval of t by binary searching
	int L = find_ind(t); //[t_L,t_(L+1)] 
	assert(L >= k - 1 && L <= n + 1);

	// P_(L-k+1),..,P_L control the interval [t_L,t_(L+1)]
	MatrixXd temp = controlPw.block(L - k + 1, 0, k, controlPw.cols());

	// de-Boor algorithm
	for (int r = 1; r <= k - 1; r++)
		for (int i = L - k + 1 + r; i <= L; i++) 
		{
			double factor = (t - knots(i)) / (knots(i + k - r) - knots(i));
			int start = i - (L - k + 1 + r);
			temp.row(start) = (1.0 - factor)*temp.row(start) + factor*temp.row(start + 1);
		}

	return temp.row(0);
}

VectorXd NURBSCurve::parameterize(const MatrixXd & points)
{
	
	int K = points.rows() - 1;
	// chord length parametrization, u_0,...,u_K
	VectorXd params(points.rows());
	params(0) = 0.0;
	for (int i = 1; i <= K; i++) {
		params(i) = params(i - 1) + (points.row(i) - points.row(i - 1)).norm();
	}
	params = params / params(K);
	params(K) = 1.0;
	return params;
}

double NURBSCurve::basis(int i,int p, double t,const VectorXd &knotvector)
{
	//cout << "knotvector:\n" << knotvector.transpose() << endl;
	//int p = knotvector.size() - 1;
	assert(p >= 1);
	if (p == 1) {
		if (t >= knotvector(i)&& t < knotvector(i+1)) {
			return 1.0;
		}
		else {
			return 0.0;
		}
	}
	
	double a = knotvector(i+p-1) - knotvector(i);
	double b = knotvector(i+p) - knotvector(i+1);
	a = (a == 0.0) ? 0.0 : (t - knotvector(i)) / a;
	b = (b == 0.0) ? 0.0 : (knotvector(i+p) - t) / b;
	return a*basis(i,p-1,t,knotvector) + b*basis(i+1,p-1,t,knotvector);
	
}

void NURBSCurve::interpolate(const MatrixXd &points)
{
	// points: P_0, ..., P_K
	//assert(points.rows() >= 4 && (points.cols() == 2 || points.cols() == 3));
	int K = points.rows() - 1;
	// chord length parametrization, u_0,...,u_K
	VectorXd params(points.rows());
	params(0) = 0.0;
	for (int i = 1; i <= K; i++) {
		params(i) = params(i-1) + (points.row(i) - points.row(i - 1)).norm();
	}
	params = params / params(K);
	params(K) = 1.0;

	knots = VectorXd::Zero(K + 7); // u_0,u_0,u_0, u_0,...,u_K, u_K,u_K,u_K
	
	knots.block(3, 0, K + 1, 1) = params;
	knots(K + 6) = 1.0; knots(K + 5) = 1.0; knots(K + 4) = 1.0;
	//cout << "knots:\n" << knots << endl;
	isRational = false;
	n = K + 2; // X_0,X_1,...,X_(K+2)
	k = 4; // order 4
	controlP = MatrixXd::Zero(K + 3, points.cols());
	controlP.row(0) = points.row(0);
	controlP.row(K + 2) = points.row(K);

	// solve linear equation: AX = b
	VectorXd d0 = (1.0 / (params(0) - params(1)) + 1.0 / (params(0) - params(2)))*points.row(0)
				+ (1.0 / (params(1) - params(0)) + 1.0 / (params(2) - params(1)))*points.row(1)
				+ (1.0 / (params(2) - params(0)) + 1.0 / (params(2) - params(1)))*points.row(2);

	VectorXd dK = (1.0 / (params(K) - params(K-1)) + 1.0 / (params(K) - params(K-2)))*points.row(K)
		+ (1.0 / (params(K-1) - params(K)) + 1.0 / (params(K-2) - params(K-1)))*points.row(K-1)
		+ (1.0 / (params(K-1) - params(K-2)) - 1.0 / (params(K) - params(K-2)))*points.row(K-2);
	MatrixXd X(K + 1, points.cols()); //X_1,...,X_(K+1)
	MatrixXd A = MatrixXd::Zero(K + 1, K + 1);
	
	MatrixXd b = points;
	b.row(0) += (params(1) - params(0)) / 3 * d0;
	b.row(K) -= (params(K) - params(K - 1)) / 3 * dK;

	A(0, 0) = A(K, K) = 1.0;
	for (int i = 1; i <= K - 1; i++) {
		A(i, i - 1) = basis(i, 4, params(i), knots);
		A(i, i) = basis(i + 1, 4, params(i), knots);
		//cout << "iiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiii:" << i << endl;
		A(i, i + 1) = basis(i + 2, 4, params(i), knots);

	}
	//cout << "A:\n" << A << endl;
	X = A.inverse()*b;
	//cout << "X:\n" << X << endl;
	controlP.block(1, 0, K + 1, points.cols()) = X;
	controlPw = controlP;

}

// interpolate with appointed knot vector
void NURBSCurve::interpolate(const MatrixXd &points, const VectorXd &knotvector)
{
	// points: P_0, ..., P_K
	//assert(points.rows() >= 4 && (points.cols() == 2 || points.cols() == 3));
	int K = points.rows() - 1;
	assert(knotvector.size() == points.rows() + 6);
	
	knots = knotvector;
	VectorXd params = knots.block(3, 0, K + 1, 1);
	/*knots.block(3, 0, K + 1, 1) = params;
	knots(K + 6) = 1.0; knots(K + 5) = 1.0; knots(K + 4) = 1.0;*/
	//cout << "knots:\n" << knots << endl;
	isRational = false;
	n = K + 2; // X_0,X_1,...,X_(K+2)
	k = 4; // order 4
	controlP = MatrixXd::Zero(K + 3, points.cols());
	controlP.row(0) = points.row(0);
	controlP.row(K + 2) = points.row(K);

	// solve linear equation: AX = b
	VectorXd d0 = (1.0 / (params(0) - params(1)) + 1.0 / (params(0) - params(2)))*points.row(0)
		+ (1.0 / (params(1) - params(0)) + 1.0 / (params(2) - params(1)))*points.row(1)
		+ (1.0 / (params(2) - params(0)) + 1.0 / (params(2) - params(1)))*points.row(2);

	VectorXd dK = (1.0 / (params(K) - params(K - 1)) + 1.0 / (params(K) - params(K - 2)))*points.row(K)
		+ (1.0 / (params(K - 1) - params(K)) + 1.0 / (params(K - 2) - params(K - 1)))*points.row(K - 1)
		+ (1.0 / (params(K - 1) - params(K - 2)) - 1.0 / (params(K) - params(K - 2)))*points.row(K - 2);
	MatrixXd X(K + 1, points.cols()); //X_1,...,X_(K+1)
	MatrixXd A = MatrixXd::Zero(K + 1, K + 1);

	MatrixXd b = points;
	b.row(0) += (params(1) - params(0)) / 3 * d0;
	b.row(K) -= (params(K) - params(K - 1)) / 3 * dK;

	A(0, 0) = A(K, K) = 1.0;
	for (int i = 1; i <= K - 1; i++) {
		A(i, i - 1) = basis(i, 4, params(i), knots);
		A(i, i) = basis(i + 1, 4, params(i), knots);
		//cout << "iiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiii:" << i << endl;
		A(i, i + 1) = basis(i + 2, 4, params(i), knots);

	}
	//cout << "A:\n" << A << endl;
	X = A.inverse()*b;
	//cout << "X:\n" << X << endl;
	controlP.block(1, 0, K + 1, points.cols()) = X;
	controlPw = controlP;
}

bool NURBSCurve::insert(double t)
{
	/*cout << "before inserting:------------------" << endl;
	cout << "controlPw:\n" << controlPw << endl;
	cout << "knots:\n" << knots << endl;*/
	assert(t >= knots(0) && t <= knots(n + k));
	int L = find_ind(t);
	// befor insert: p_0, ..., p_(L-k+1), p_(L-k+2),... ,		 p_L,...
	// after insert: p_0, ..., p_(L-k+1), p'_(L-k+2),..., p'_L,  p'_(L+1)=p_L,...   
	int start = (L-k+1>=0)?L-k+1:0;
	int end = (L<=n)?L:n;
	MatrixXd new_controlPw(controlPw.rows()+1,controlPw.cols());
	for(int i=0;i<new_controlPw.rows();i++){
		if(i<=start){
			new_controlPw.row(i) = controlPw.row(i);
		}else if(i<=end){
			double factor = (t-knots(i))/(knots(i+k-1)-knots(i));
			new_controlPw.row(i)=factor*controlPw.row(i)+(1.0-factor)*controlPw.row(i-1);
			
		}else{
			new_controlPw.row(i)=controlPw.row(i-1);

		}
	}

	VectorXd new_knots(knots.size()+1);
	for(int i=0;i<new_knots.size();i++){
		if(i<=L){
			new_knots(i)=knots(i);
		}else if(i==L+1){
			new_knots(i)=t;
		}else{
			new_knots(i)=knots(i-1);
		}
	}

	controlPw = new_controlPw;
	knots = new_knots;
	n+=1;
	/*cout << "after inserting:------------------" << endl;
	cout << "controlPw:\n" << controlPw << endl;
	cout << "knots:\n" << knots << endl;*/

	return true;
}


// draw controlpolygon
void NURBSCurve::drawControlPolygon(igl::opengl::glfw::Viewer &viewer){
	
	// plot control points
	viewer.data().add_points(controlP, Eigen::RowVector3d(1, 1, 1));
	// plot control polygon
	for (int i = 0; i < n; i++){
		viewer.data().add_edges(
			controlP.row(i),
			controlP.row(i + 1),
			Eigen::RowVector3d(1, 1, 1));
	}
}

// draw NURBS surface
void NURBSCurve::drawSurface(igl::opengl::glfw::Viewer &viewer, double resolution){
	double left = knots(k - 1);
	double right = knots(n + 1);
	const int num = (right - left) / resolution;
	double new_resolution = (right - left) / num;
	// plot NURBSCurve
	for (int i = 0; i < num; i++)
	{
		Eigen::MatrixXd P1 = eval(left + 1.0*i* new_resolution);
		Eigen::MatrixXd P2 = eval(left + 1.0*(i + 1)*new_resolution);

		if (isRational) {
			viewer.data().add_edges(
				P1.rowwise().hnormalized(),
				P2.rowwise().hnormalized(),
				Eigen::RowVector3d(1, 0, 0));
		}
		else {
			viewer.data().add_edges(P1, P2, Eigen::RowVector3d(1, 0, 0));
		}

	}
}


// display by libigl
void NURBSCurve::draw(
	igl::opengl::glfw::Viewer& viewer, 
	bool showpolygon,bool showsurface,
	double resolution)
{
	if(isRational){
		controlP = controlPw.rowwise().hnormalized();
	}else{
		controlP = controlPw;
	}

	if(showpolygon){
		drawControlPolygon(viewer);
	}
	if(showsurface){
		drawSurface(viewer);
	}
	viewer.core.align_camera_center(controlP);
}





