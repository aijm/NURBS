#include <igl/opengl/glfw/Viewer.h>
#include "BezierCurve.h"
#include "NURBSCurve.h"
#include "NURBSSurface.h"
void testInterpolate(igl::opengl::glfw::Viewer &viewer)
{
	NURBSCurve nurbs;

	/*VectorXd knots(4);
	knots << 0, 0.25, 0.5, 0.75;
	cout << "value: " << nurbs.basis(knots, 0) << endl;
	cout << "value: " << nurbs.basis(knots, 0.3) << endl;
	cout << -1.5 + 12 * 0.3 - 16 * 0.3*0.3 << endl;*/

	MatrixXd points = (MatrixXd(9, 3) <<
		1.0, 0.0, 0.0,
		1.0, 1.0, 0.0,
		0.0, 1.0, 0.0,
		-1.0, 1.0, 0.0,
		-1.0, 0.0, 0.0,
		-1.0, -1.0, 0.0,
		0.0, -1.0, 0.0,
		1.0, -1.0, 0.0,
		1.0, 0.0, 0.0).finished();
	nurbs.interpolate(points);
	nurbs.show(viewer, 0.01);
	viewer.data().add_points(points,RowVector3d(0.0,1.0,0.0));
	viewer.core.align_camera_center(points);

	
}
void testCube(igl::opengl::glfw::Viewer &viewer)
{
	// Inline mesh of a cube
	const Eigen::MatrixXd V= (Eigen::MatrixXd(8,3)<<
	0.0,0.0,0.0,
	0.0,0.0,1.0,
	0.0,1.0,0.0,
	0.0,1.0,1.0,
	1.0,0.0,0.0,
	1.0,0.0,1.0,
	1.0,1.0,0.0,
	1.0,1.0,1.0).finished();
	const Eigen::MatrixXi F = (Eigen::MatrixXi(12,3)<<
	1,7,5,
	1,3,7,
	1,4,3,
	1,2,4,
	3,8,7,
	3,4,8,
	5,7,8,
	5,8,6,
	1,5,6,
	1,6,2,
	2,6,8,
	2,8,4).finished().array()-1;

	viewer.data().set_mesh(V, F);
	viewer.data().set_face_based(true);
}
void testBezierCurve(igl::opengl::glfw::Viewer &viewer)
{
	// test BezierCurve
	Eigen::MatrixXd controlP(5, 2);
	controlP <<
		0.0, 0.0,
		1.0, 1.0,
		2.0, -1.0,
		3.0, 0.0,
		4.0, 2.0;
	BezierCurve curve(controlP);
	curve.show(viewer);
	viewer.core.align_camera_center(controlP);
}

void testNURBSCurve(igl::opengl::glfw::Viewer &viewer)
{
	Eigen::MatrixXd controlP1(4, 4);
	controlP1 <<
		-1.0, 0.0, 0.0, 1.0,
		-1.0 / 3, 2.0 / 3, 0.0, 1.0 / 3,
		1.0 / 3, 2.0 / 3, 0.0, 1.0 / 3,
		1.0, 0.0, 0.0, 1.0;
	//cout << "controlP: " << controlP1 << endl;
	Eigen::VectorXd knots1(8);
	knots1 <<
		0.0,0.0,0.0,0.0,1.0,1.0, 1.0, 1.0;
	
	//cout << "hahahaha" << endl;
	NURBSCurve nurbs1(3, 4, controlP1, knots1,true);

	Eigen::MatrixXd controlP2(7, 3);
	controlP2 <<
		1.0, 0.0, 1.0,
		0.5, 0.5, 0.5,
		-0.5, 0.5, 0.5,
		-1.0, 0.0, 1.0,
		-0.5, -0.5, 0.5,
		0.5, -0.5, 0.5,
		1.0, 0.0, 1.0;
	

	Eigen::VectorXd knots2(10);
	knots2 <<
		0.0, 0.0, 0.0, 0.25, 0.5, 0.5, 0.75, 1.0, 1.0, 1.0;
	
	NURBSCurve nurbs2(6, 3, controlP2, knots2, true);


	nurbs2.show(viewer,0.01);
	viewer.core.align_camera_center(nurbs2.controlP);
}

void testCylindr(igl::opengl::glfw::Viewer &viewer)
{
	vector<MatrixXd> controlPolygon(2);
	controlPolygon[0] = (MatrixXd(9, 4) <<
		1.0, 0.0, 0.0, 1.0,
		0.7071, 0.7071, 0.0, 0.7071,
		0.0, 1.0, 0.0, 1.0,
		-0.7071, 0.7071, 0.0, 0.7071,
		-1.0, 0.0, 0.0, 1.0,
		-0.7071, -0.7071, 0.0, 0.7071,
		0.0, -1.0, 0.0, 1.0,
		0.7071, -0.7071, 0.0, 0.7071,
		1.0, 0.0, 0.0, 1.0).finished();

	controlPolygon[1] = (MatrixXd(9, 4) <<
		1.0, 0.0, 1.0, 1.0,
		0.7071, 0.7071, 0.7071, 0.7071,
		0.0, 1.0, 1.0, 1.0,
		-0.7071, 0.7071, 0.7071, 0.7071,
		-1.0, 0.0, 1.0, 1.0,
		-0.7071, -0.7071, 0.7071, 0.7071,
		0.0, -1.0, 1.0, 1.0,
		0.7071, -0.7071, 0.7071, 0.7071,
		1.0, 0.0, 1.0, 1.0).finished();
	
	Eigen::VectorXd uknots(12);
	Eigen::VectorXd vkonts(4);
	uknots << 0, 0, 0, 0.25, 0.25, 0.5, 0.5, 0.75, 0.75, 1, 1, 1;
	vkonts << 0, 0, 1, 1;

	VectorXi order(2);
	order << 3, 2;
	
	
	NURBSSurface nurbs(order, controlPolygon, uknots, vkonts, true);

	nurbs.show(viewer, 0.005);
	cout << "hahahahah" << endl;
	viewer.core.align_camera_center(nurbs.mesh_V, nurbs.mesh_F);
}

void testTorus(igl::opengl::glfw::Viewer &viewer)
{
	vector<MatrixXd> controlPolygon(9);

	{
		controlPolygon[0] = (MatrixXd(9, 4) <<
			-125.0, 0.0, 0.0, 1.0,
			-125.0, 0.0, 125.0, 1.0,
			0.0, 0.0, 125.0, 1.0,
			125.0, 0.0, 125.0, 1.0,
			125.0, 0.0, 0.0, 1.0,
			125.0, 0.0, -125.0, 1.0,
			0.0, 0.0, -125.0, 1.0,
			-125.0, 0.0, -125.0, 1.0,
			-125.0, 0.0, 0.0, 1.0).finished();

		controlPolygon[1] = (MatrixXd(9, 4) <<
			-88.375, 70.7, 0.0, 0.707,
			-88.375, 70.7, 88.375, 0.707,
			0.0, 70.7, 88.375, 0.707,
			88.375, 70.7, 88.375, 0.707,
			88.375, 70.7, 0.0, 0.707,
			88.375, 70.7, -88.375, 0.707,
			0.0, 70.7, -88.375, 0.707,
			-88.375, 70.7, -88.375, 0.707,
			-88.375, 70.7, 0.0, 0.707).finished();

		controlPolygon[2] = (MatrixXd(9, 4) <<
			-225.0, 100.0, 0.0, 1.0,
			-225.0, 100.0, 225.0, 1.0,
			0.0, 100.0, 225.0, 1.0,
			225.0, 100.0, 225.0, 1.0,
			225.0, 100.0, 0.0, 1.0,
			225.0, 100.0, -225.0, 1.0,
			0.0, 100.0, -225.0, 1.0,
			-225.0, 100.0, -225.0, 1.0,
			-225.0, 100.0, 0.0, 1.0).finished();

		controlPolygon[3] = (MatrixXd(9, 4) <<
			-229.775, 70.7, 0.0, 0.707,
			-229.775, 70.7, 229.775, 0.707,
			0.0, 70.7, 229.775, 0.707,
			229.775, 70.7, 229.775, 0.707,
			229.775, 70.7, 0.0, 0.707,
			229.775, 70.7, -229.775, 0.707,
			0.0, 70.7, -229.775, 0.707,
			-229.775, 70.7, -229.775, 0.707,
			-229.775, 70.7, 0.0, 0.707).finished();

		controlPolygon[4] = (MatrixXd(9, 4) <<
			-325.0, 0.0, 0.0, 1.0,
			-325.0, 0.0, 325.0, 1.0,
			0.0, 0.0, 325.0, 1.0,
			325.0, 0.0, 325.0, 1.0,
			325.0, 0.0, 0.0, 1.0,
			325.0, 0.0, -325.0, 1.0,
			0.0, 0.0, -325.0, 1.0,
			-325.0, 0.0, -325.0, 1.0,
			-325.0, 0.0, 0.0, 1.0).finished();

		controlPolygon[5] = (MatrixXd(9, 4) <<
			-229.775, -70.7, 0.0, 0.707,
			-229.775, -70.7, 229.775, 0.707,
			0.0, -70.7, 229.775, 0.707,
			229.775, -70.7, 229.775, 0.707,
			229.775, -70.7, 0.0, 0.707,
			229.775, -70.7, -229.775, 0.707,
			0.0, -70.7, -229.775, 0.707,
			-229.775, -70.7, -229.775, 0.707,
			-229.775, -70.7, 0.0, 0.707).finished();

		controlPolygon[6] = (MatrixXd(9, 4) <<
			-225.0, -100.0, 0.0, 1.0,
			-225.0, -100.0, 225.0, 1.0,
			0.0, -100.0, 225.0, 1.0,
			225.0, -100.0, 225.0, 1.0,
			225.0, -100.0, 0.0, 1.0,
			225.0, -100.0, -225.0, 1.0,
			0.0, -100.0, -225.0, 1.0,
			-225.0, -100.0, -225.0, 1.0,
			-225.0, -100.0, 0.0, 1.0).finished();

		controlPolygon[7] = (MatrixXd(9, 4) <<
			-88.375, -70.7, 0.0, 0.707,
			-88.375, -70.7, 88.375, 0.707,
			0.0, -70.7, 88.375, 0.707,
			88.375, -70.7, 88.375, 0.707,
			88.375, -70.7, 0.0, 0.707,
			88.375, -70.7, -88.375, 0.707,
			0.0, -70.7, -88.375, 0.707,
			-88.375, -70.7, -88.375, 0.707,
			-88.375, -70.7, 0.0, 0.707).finished();

		controlPolygon[8] = (MatrixXd(9, 4) <<
			-125.0, 0.0, 0.0, 1.0,
			-125.0, 0.0, 125.0, 1.0,
			0.0, 0.0, 125.0, 1.0,
			125.0, 0.0, 125.0, 1.0,
			125.0, 0.0, 0.0, 1.0,
			125.0, 0.0, -125.0, 1.0,
			0.0, 0.0, -125.0, 1.0,
			-125.0, 0.0, -125.0, 1.0,
			-125.0, 0.0, 0.0, 1.0).finished();
	}
	
	Eigen::VectorXd uknots(12);
	Eigen::VectorXd vkonts(12);
	uknots << 0, 0, 0, 0.25, 0.25, 0.5, 0.5, 0.75, 0.75, 1, 1, 1;
	vkonts << 0, 0, 0, 0.25, 0.25, 0.5, 0.5, 0.75, 0.75, 1, 1, 1;

	VectorXi order(2);
	order << 3, 3;


	NURBSSurface nurbs(order, controlPolygon, uknots, vkonts, true);

	nurbs.show(viewer, 0.005);
	// Use the z coordinate as a scalar field over the surface
	//Eigen::VectorXd Z = nurbs.mesh_V.col(2);
	//Eigen::MatrixXd C;
	//// Compute per-vertex colors
	//igl::jet(Z, true, C);

	//// Add per-vertex colors
	//viewer.data().set_colors(C);
	cout << "hahahahah" << endl;
	//viewer.core.align_camera_center(nurbs.mesh_V, nurbs.mesh_F);
}

void testSurface1(igl::opengl::glfw::Viewer &viewer)
{	
	vector<MatrixXd> controlPolygon(4);
	/*controlPolygon[0] = (MatrixXd(6, 3) <<
		-25, -25, -10,
		-25, -15, -5,
		-25, -5, 0,
		-25, 5, 0,
		-25, 15, -5,
		-25, 25, -10).finished();*/

	controlPolygon[0] = (MatrixXd(4, 3) <<
		//-15, -25, -8,
		-15, -15, -4,
		-15, -5, -4,
		-15, 5, -4,
		-15, 15, -4/*,
		-15, 25, -8*/).finished();

	controlPolygon[1] = (MatrixXd(4, 3) <<
		//-5, -25, -5,
		-5, -15, -3,
		-5, -5, -8, 
		-5, 5, -8, 
		-5, 15, -3/*,
		-5, 25, -5*/).finished();

	controlPolygon[2] = (MatrixXd(4, 3) <<
		//5, -25, -3,
		5, -15, -2,
		5, -5, -8,
		5, 5, -8,
		5, 15, -2/*,
		5, 25, -3*/).finished();

	controlPolygon[3] = (MatrixXd(4, 3) <<
		//15, -25, -8,
		15, -15, -4,
		15, -5, -4,
		15, 5, -4,
		15, 15, -4/*,
		15, 25, -8*/).finished();

	/*controlPolygon[5] = (MatrixXd(6, 3) <<
		25, -25, -10,
		25, -15, -5,
		25, -5, 2,
		25, 5, 2,
		25, 15, -5,
		25, 25, -10).finished();*/

	Eigen::VectorXd uknots(8);
	Eigen::VectorXd vkonts(8);
	uknots << -1.0, -1.0, 0.0, 1.0, 2.0, 3.0, 4.0, 4.0;
	vkonts << -1.0, -1.0, 0.0, 1.0, 2.0, 3.0, 4.0, 4.0;

	VectorXi order(2);
	order << 4, 4;


	NURBSSurface nurbs(order, controlPolygon, uknots, vkonts, false);

	nurbs.show(viewer, 1);
	cout << "hahahahah" << endl;
	viewer.core.align_camera_center(nurbs.mesh_V, nurbs.mesh_F);
	cout << "nurbs: \n" << nurbs.mesh_V << endl;
}

int main(int argc, char *argv[])
{
  

  // Plot the mesh
  igl::opengl::glfw::Viewer viewer;
  //testBezierCurve(viewer);
  //testCube(viewer);
  //testNURBSCurve(viewer);
  //testCylindr(viewer);
  //testTorus(viewer);
  //testSurface1(viewer);
  testInterpolate(viewer);
  viewer.launch();
}
