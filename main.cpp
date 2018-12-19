#include <igl/opengl/glfw/Viewer.h>
#include "BezierCurve.h"
#include "NURBSCurve.h"
#include "NURBSSurface.h"
#include "test.h"

#include <igl/opengl/glfw/imgui/ImGuiMenu.h>
#include <igl/opengl/glfw/imgui/ImGuiHelpers.h>
#include <imgui/imgui.h>


//NURBSCurve nurbs;
NURBSSurface nurbs;
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

	/*curves[0].draw(viewer, true, true);
	curves[1].draw(viewer, true, true);
	curves[2].draw(viewer, true, true);
	curves[3].draw(viewer, true, true);
	curves[4].draw(viewer, true, true);
	curves[5].draw(viewer, true, true);
	curves[6].draw(viewer, true, true);
	curves[7].draw(viewer, true, true);*/
	//NURBSSurface nurbs;
	//vector<NURBSCurve> new_curves(8);
	nurbs.skinning(curves, viewer);
	
	nurbs.draw(viewer,false,true);
	//nurbs.saveNURBS("lamp_skinning");
}
void testCircleSkinning(igl::opengl::glfw::Viewer &viewer)
{
	vector<NURBSCurve> circles(3);
	circles[0].loadNURBS("../curve.cptw");
	circles[1].loadNURBS("../curve1.cptw");
	circles[2].loadNURBS("../curve2.cptw");
	circles[0].draw(viewer,false,true);
	circles[1].draw(viewer, false,true);
	circles[2].draw(viewer, false,true);

	nurbs.skinning(circles,viewer);
	nurbs.draw(viewer,true,true);

	
}
void testPiafit(igl::opengl::glfw::Viewer &viewer)
{
	MatrixXd points;
	loadpoints("../g.cur", points);
	
	//viewer.data().add_points(points, RowVector3d(0, 1, 0));
	//nurbs.piafit(points, 100, 1e-5);
	//nurbs.piafit(nurbs_load.controlPw);
	
	nurbs.draw(viewer, false, true,0.0001);
}

void testLSPIAfit(igl::opengl::glfw::Viewer &viewer)
{
	MatrixXd points;
	loadpoints("../curve1.cpt", points);
	
	viewer.data().add_points(points, RowVector3d(0, 1, 0));
	//viewer.data().add_points(points, RowVector3d(0, 1, 0));
	//nurbs.lspiafit(points, 5);
	//nurbs.piafit(nurbs_load.controlPw);

	nurbs.draw(viewer, true, true, 0.0001);
}

int main(int argc, char *argv[])
{
  

  // Plot the mesh
  igl::opengl::glfw::Viewer viewer;
  //testPiafit(viewer);
  //testLSPIAfit(viewer);
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
  testLampSkinning(viewer);

  // Attach a menu plugin
  igl::opengl::glfw::imgui::ImGuiMenu menu;
  viewer.plugins.push_back(&menu);

  // Customize the menu
  double doubleVariable = 0.1f; // Shared between two menus

								// Add content to the default menu window
  menu.callback_draw_viewer_menu = [&]()
  {
	  // Draw parent menu content
	  menu.draw_viewer_menu();

	  // Add new group
	  if (ImGui::CollapsingHeader("New Group", ImGuiTreeNodeFlags_DefaultOpen))
	  {
		  // Expose variable directly ...
		  ImGui::InputDouble("double", &doubleVariable, 0, 0, "%.4f");

		  // ... or using a custom callback
		  static bool boolVariable = true;
		  if (ImGui::Checkbox("bool", &boolVariable))
		  {
			  // do something
			  std::cout << "boolVariable: " << std::boolalpha << boolVariable << std::endl;
		  }

		  // Expose an enumeration type
		  enum Orientation { Up = 0, Down, Left, Right };
		  static Orientation dir = Up;
		  ImGui::Combo("Direction", (int *)(&dir), "Up\0Down\0Left\0Right\0\0");

		  // We can also use a std::vector<std::string> defined dynamically
		  static int num_choices = 3;
		  static std::vector<std::string> choices;
		  static int idx_choice = 0;
		  if (ImGui::InputInt("Num letters", &num_choices))
		  {
			  num_choices = std::max(1, std::min(26, num_choices));
		  }
		  if (num_choices != (int)choices.size())
		  {
			  choices.resize(num_choices);
			  for (int i = 0; i < num_choices; ++i)
				  choices[i] = std::string(1, 'A' + i);
			  if (idx_choice >= num_choices)
				  idx_choice = num_choices - 1;
		  }
		  ImGui::Combo("Letter", &idx_choice, choices);

		  // Add a button
		  if (ImGui::Button("Print Hello", ImVec2(-1, 0)))
		  {
			  std::cout << "Hello\n";
		  }
	  }
  };

  // Draw additional windows
  menu.callback_draw_custom_window = [&]()
  {
	  // Define next window position + size
	  ImGui::SetNextWindowPos(ImVec2(180.f * menu.menu_scaling(), 10), ImGuiSetCond_FirstUseEver);
	  ImGui::SetNextWindowSize(ImVec2(200, 160), ImGuiSetCond_FirstUseEver);
	  ImGui::Begin(
		  "New Window", nullptr,
		  ImGuiWindowFlags_NoSavedSettings
	  );


	  // Expose the same variable directly ...
	  ImGui::PushItemWidth(-80);
	  ImGui::DragScalar("double", ImGuiDataType_Double, &doubleVariable, 0.1, 0, 0, "%.4f");
	  ImGui::PopItemWidth();

	  static std::string str = "bunny";
	  ImGui::InputText("Name", str);

	  ImGui::End();
  };

  viewer.callback_key_down = &key_down;
  //nurbs.draw(viewer, showpolygon, showsurface);
  viewer.launch();

  return 0;
}
