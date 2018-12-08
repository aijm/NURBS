#include <igl/opengl/glfw/Viewer.h>
#include "BezierCurve.h"
#include "NURBSCurve.h"
#include "NURBSSurface.h"
#include "test.h"


NURBSCurve nurbs;
double resolution = 0.01;

bool showpolygon = true;
bool showsurface = true;

void insert_loop(igl::opengl::glfw::Viewer &viewer) {
	double s = 0.0;
	double t = 0.0;
	while (true) {
		cout << "insert kont, format: s t" << endl;
		
		if (!(cin >> s/* >> t*/)) {
			cin.clear(); //clear the buffer
			cin.get();
			cout << "error! please use right format!" << endl;
			continue;
		}
		else {
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




int main(int argc, char *argv[])
{
  

  // Plot the mesh
  igl::opengl::glfw::Viewer viewer;
  //testBezierCurve(viewer);
  //testCube(viewer);
//   testNURBSCurve(viewer);
//   testCylindr(viewer);
//   testTorus(viewer);
//   testSurface1(viewer);
  //testInterpolate(viewer);
  //nurbs.loadNURBS("circle.cptw");
  //nurbs.draw(viewer);
  cout<< nurbs.loadNURBS("curve1.cpt")<<endl;
  viewer.callback_key_down = &key_down;
  nurbs.draw(viewer,showpolygon, showsurface);
  viewer.launch();
}
