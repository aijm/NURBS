// This file is part of NURBS, a simple NURBS library.
// github repo: https://github.com/aijm/NURBS
// Copyright (C) 2018 Jiaming Ai <aichangeworld@gmail.com>

#include "NURBSSurface.h"

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
NURBSSurface::NURBSSurface(
	VectorXi _order, 
	vector<MatrixXd> _controlP, 
	VectorXd _uknots, 
	VectorXd _vknots,
	bool _isRational)
{
	// pass and check the parameters of NURBS surface
	u_order = _order(0);
	v_order = _order(1);
	v_num = _controlP.size() - 1;
	assert(u_order >= 1 && v_order >= 1 && v_num >= 1);
	u_num = _controlP[0].rows() - 1;
	assert(_uknots.size() == u_num + u_order + 1 && _vknots.size() == v_num + v_order + 1);
	assert(u_num >= u_order - 1 && v_num >= v_order - 1);

	dimension = _controlP[0].cols(); // the dimension of control point 2 or 3 or 4(weighted)
	uknots = _uknots;
	vknots = _vknots;
	controlPw = _controlP;
	isRational = _isRational;

	if (isRational) { assert(dimension == 4); }
	else { assert(dimension == 3); }

}


// load
bool NURBSSurface::loadNURBS(string name){
	ifstream in(name);
	if(!in){
		return false;
	}
	char sep;
	dimension=3;
	in>> isRational;
	in>> u_num >> v_num;
	in>> u_order >> v_order;
	in>> dimension;
	controlPw = vector<MatrixXd>(v_num+1);
	for(int i=0;i<controlPw.size();i++){
		controlPw[i] = MatrixXd(u_num+1,dimension);
	}
	for(int i=0;i<controlPw.size();i++){
		for(int j=0;j<=u_num;j++){
			for(int k=0;k<dimension;k++){
				in>> controlPw[i](j,k);
			}
			in>>sep;
		}
	}

	uknots = VectorXd(u_num+u_order+1);
	vknots = VectorXd(v_num+v_order+1);

	for(int i=0;i<uknots.size();i++){
		in>>uknots(i);
	}
	
	for(int i=0;i<vknots.size();i++){
		in>>vknots(i);
	}
	return true;
}

// save 
bool NURBSSurface::saveNURBS(string name){
	if(controlPw.size()==0){
		cout<< "nothing to save!"<<endl;
		return false;
	}
	if(isRational){
		name+=".cptw";
	}else{
		name+=".cpt";
	}
	ofstream out(name);
	if(!out){
		return false;
	}
	IOFormat outputFmt(4, 0, " ", " ", "", ",");	
	out<< isRational<<endl;
	out<< u_num << " " <<v_num<<endl;
	out<< u_order << " "<< v_order<<endl;
	out<< controlPw[0].cols()<<endl;
	
	for(int i=0;i<controlPw.size();i++){
		for(int j=0;j<=u_num;j++){
			out<<controlPw[i].row(j).format(outputFmt);
		}
		out<<endl;
	}
	out<<uknots.transpose()<<endl;
	out<<vknots.transpose();
	return true;
}
	

// find the knot interval of t by binary searching
int NURBSSurface::find_ind(double t, int k, int n, const VectorXd& knots)
{
	if (t == knots(n + 1)) return n;
	int low = k - 1;
	int high = n + 1;
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
// calculate coordinate of curve point with parameter u & v
MatrixXd NURBSSurface::eval(double u, double v)
{
	//Calculating the Control Points of U-direction Isoparametric Line
	MatrixXd v_controlP(v_num + 1, dimension); 
	for (int i = 0; i <= v_num; i++)
	{
		v_controlP.row(i) = eval(u, controlPw[i], uknots);
	}
	return eval(v, v_controlP, vknots); // Calculating the coordinate of parameter v
}

MatrixXd NURBSSurface::eval(
	double t, 
	const MatrixXd &_controlP, 
	const VectorXd &knots)
{
	int n = _controlP.rows() - 1;
	int k = knots.size() - _controlP.rows();
	// find the knot interval of t by binary searching
	int L = find_ind(t, k, n, knots); //[t_L,t_(L+1)] 

	// P_(L-k+1),..,P_L control the interval [t_L,t_(L+1)]
	MatrixXd temp = _controlP.block(L - k + 1, 0, k, _controlP.cols());

	for (int r = 1; r <= k - 1; r++)
		for (int i = L - k + 1 + r; i <= L; i++) 
		{
			double factor = (t - knots(i)) / (knots(i + k - r) - knots(i));
			int start = i - (L - k + 1 + r);
			temp.row(start) = (1.0 - factor)*temp.row(start) + factor*temp.row(start + 1);
		}
	MatrixXd curvePoint = temp.row(0);
	return curvePoint;
}


// kont insertion
bool NURBSSurface::insert(double s, double t){
	return false;
}

// draw controlpolygon
void NURBSSurface::drawControlPolygon(igl::opengl::glfw::Viewer &viewer){
	
	//plot control points and control polygon
	
	//plot control points
	for (int i = 0; i < controlP.size(); i++)
	{
		viewer.data().add_points(
			controlP[i],
			Eigen::RowVector3d(1, 0, 0));
	}
	// plot control polygon
	for (int j = 0; j <= v_num; j++)
		for (int i = 0; i <= u_num; i++)
		{
			if (i != u_num) {
				viewer.data().add_edges(
					controlP[j].row(i),
					controlP[j].row(i + 1),
					Eigen::RowVector3d(1, 0, 0));
			}
			if (j != v_num) {
				viewer.data().add_edges(
					controlP[j].row(i),
					controlP[j + 1].row(i),
					Eigen::RowVector3d(1, 0, 0));
			}
		}
}

// draw NURBS surface
void NURBSSurface::drawSurface(igl::opengl::glfw::Viewer &viewer, double resolution){
	// cut apart the parameter domain
	double u_low = uknots(u_order - 1);
	double u_high = uknots(u_num + 1);
	const int uspan = (u_high - u_low) / resolution;
	double u_resolution = (u_high - u_low) / uspan;

	double v_low = vknots(v_order - 1);
	double v_high = vknots(v_num + 1);
	const int vspan = (v_high - v_low) / resolution;
	double v_resolution = (v_high - v_low) / vspan;

	mesh_V = MatrixXd((uspan + 1)*(vspan + 1), 3);
	mesh_F = MatrixXi(2 * uspan*vspan, 3);
	// discretize NURBS Surface into triangular mesh(V,F) in libigl mesh structure
	// calculate mesh_V
	for (int j = 0; j <= vspan; j++)
		for (int i = 0; i <= uspan; i++)
		{
			RowVectorXd curvePoint = eval(u_low + i*u_resolution, v_low + j*v_resolution).row(0);
			if (isRational) { mesh_V.row(j*(uspan + 1) + i) = curvePoint.hnormalized(); }
			else { mesh_V.row(j*(uspan + 1) + i) = curvePoint; }
		}

	for (int j = 0; j<vspan; j++)
		for (int i = 0; i < uspan; i++)
		{
			int V_index = j*(uspan + 1) + i;
			int F_index = 2 * j*uspan + 2 * i;
			mesh_F.row(F_index) << V_index, V_index + 1, V_index + uspan + 1;
			mesh_F.row(F_index + 1) << V_index + uspan + 1, V_index + 1, V_index + uspan + 2;
		}

	viewer.data().set_mesh(mesh_V, mesh_F);
}


void NURBSSurface::draw(
	igl::opengl::glfw::Viewer &viewer, 
	bool showpolygon,bool showsurface,
	double resolution){

	if(controlP.size()!=controlPw.size()){
		controlP = vector<MatrixXd>(controlPw.size());
	}
	if(isRational){
		for(int i=0;i<controlP.size();i++){
			controlP[i] = controlPw[i].rowwise().hnormalized();
		}
	}else{
		controlP = controlPw;
	}
	
	if(showpolygon){
		drawControlPolygon(viewer);
		viewer.core.align_camera_center(controlP[controlP.size()/2]);
	}
	if(showsurface){
		drawSurface(viewer,resolution);
		viewer.core.align_camera_center(mesh_V,mesh_F);
	}

}