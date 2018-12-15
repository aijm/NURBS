#include <igl/opengl/glfw/Viewer.h>
#include "BezierCurve.h"
#include "NURBSCurve.h"
#include "NURBSSurface.h"
#include "test.h"


NURBSCurve nurbs;
//NURBSSurface nurbs;
double resolution = 0.01;

bool showpolygon = true;
bool showsurface = true;

void insert_loop(igl::opengl::glfw::Viewer &viewer) {
	double s = 0.0;
	double t = 0.0;
	while (true) {
		cout << "insert kont, format: s t" << endl;
		
		if (!(cin >> s /*>> t*/)) {
			cin.clear(); //clear the buffer
			cin.get();
			cout << "error! please use right format!" << endl;
			continue;
		}
		else {
			//cout << "insert : " << s << endl;
			nurbs.insert(s/*, t*/);
			//mesh.drawTmesh(viewer);
			//mesh.drawControlpolygon(viewer);
			//mesh.drawSurface(viewer);
			nurbs.draw(viewer, showpolygon, showsurface);
			return;
		}
	}
}
// This function is called every time a keyboard button is pressed
bool key_down(igl::opengl::glfw::Viewer& viewer, unsigned char key, int modifier) {
	//std::cout << "Key: " << key << " " << (unsigned int)key << std::endl;
	if (key == 32) {
		viewer.data().clear();
		insert_loop(viewer);
	}else if (key == 'P') {
		viewer.data().clear();
		showpolygon = !showpolygon;
		nurbs.draw(viewer, showpolygon, showsurface);
	}else if (key == 'S') {
		viewer.data().clear();
		showsurface = !showsurface;
		nurbs.draw(viewer, showpolygon, showsurface);
	}
	return false;
}
void testLampSkinning(igl::opengl::glfw::Viewer &viewer)
{
	vector<NURBSCurve> curves(8);
	curves[0].loadNURBS("../lamp1.cptw");
	curves[1].loadNURBS("../lamp2.cptw");
	curves[2].loadNURBS("../lamp3.cptw");
	curves[3].loadNURBS("../lamp4.cptw");
	curves[4].loadNURBS("../lamp5.cptw");
	curves[5].loadNURBS("../lamp6.cptw");
	curves[6].loadNURBS("../lamp7.cptw");
	curves[7].loadNURBS("../lamp8.cptw");

	curves[0].draw(viewer, true, true);
	curves[1].draw(viewer, true, true);
	curves[2].draw(viewer, true, true);
	curves[3].draw(viewer, true, true);
	curves[4].draw(viewer, true, true);
	curves[5].draw(viewer, true, true);
	curves[6].draw(viewer, true, true);
	curves[7].draw(viewer, true, true);
	//NURBSSurface nurbs;
	//vector<NURBSCurve> new_curves(8);
	//nurbs.skinning(curves, viewer);
	
	nurbs.draw(viewer,true,true);
	//nurbs.saveNURBS("lamp_skinning");
}
void testCircleSkinning(igl::opengl::glfw::Viewer &viewer)
{
	vector<NURBSCurve> circles(3);
	circles[0].loadNURBS("../circle.cptw");
	circles[1].loadNURBS("../circle1.cptw");
	circles[2].loadNURBS("../circle2.cptw");
	circles[0].draw(viewer,false,true);
	circles[1].draw(viewer,false,true);
	circles[2].draw(viewer,false,true);

	//nurbs.skinning(circles,viewer);
	nurbs.draw(viewer,true,true);

	
}
void testPiafit(igl::opengl::glfw::Viewer &viewer)
{
	MatrixXd points;
	loadpoints("../g.cur", points);
	
	//viewer.data().add_points(points, RowVector3d(0, 1, 0));
	nurbs.piafit(points, 100, 1e-5);
	//nurbs.piafit(nurbs_load.controlPw);
	
	nurbs.draw(viewer, false, true,0.0001);
}

void testLSPIAfit(igl::opengl::glfw::Viewer &viewer)
{
	MatrixXd points;
	loadpoints("../curve1.cpt", points);
	
	viewer.data().add_points(points, RowVector3d(0, 1, 0));
	//viewer.data().add_points(points, RowVector3d(0, 1, 0));
	nurbs.lspiafit(points, 5);
	//nurbs.piafit(nurbs_load.controlPw);

	nurbs.draw(viewer, true, true, 0.0001);
}

int main(int argc, char *argv[])
{
  

  // Plot the mesh
  igl::opengl::glfw::Viewer viewer;
  //testPiafit(viewer);
  testLSPIAfit(viewer);
  //testBezierCurve(viewer);
  //testCube(viewer);
//   testNURBSCurve(viewer);
//   testCylindr(viewer);
//   testTorus(viewer);
  //testBasisFunction();
   //testSurface1(viewer);
  //testInterpolate(viewer);
  //nurbs.loadNURBS("circle.cptw");
  //nurbs.draw(viewer);
  //testCircleSkinning(viewer);
  //testLampSkinning(viewer);

  //cout<< nurbs.loadNURBS("../circle1.cptw")<<endl;
  //viewer.data().add_label(Vector3d(1.0, 0.0, 1.0), "haha");
 /* viewer.callback_key_down = &key_down;
  nurbs.draw(viewer,showpolygon, showsurface);*/
  viewer.launch();
}
