#include "line_cylinders.h"
#include "scene.h"
#include <igl/copyleft/tetgen/tetrahedralize.h>
#include <igl/opengl/glfw/Viewer.h>
#include <igl/opengl/glfw/imgui/ImGuiHelpers.h>
#include <igl/opengl/glfw/imgui/ImGuiMenu.h>
#include <imgui/imgui.h>
#include <iostream>

Eigen::MatrixXd V;
Eigen::MatrixXi F;
igl::opengl::glfw::Viewer mgpViewer;

float currTime = 0;

//initial values
float timeStep = 0.02;
float CRCoeff = 1.0;

double tolerance = 10e-3;
int maxIterations = 10000;

Scene scene;

Eigen::MatrixXd platV;
Eigen::MatrixXi platF;
Eigen::MatrixXi platT;
Eigen::RowVector3d platCOM;
Eigen::RowVector4d platOrientation;

void createPlatform()
{
    double platWidth = 100.0;
    platCOM << 0.0, -5.0, -0.0;
    platV.resize(8, 3);
    platF.resize(12, 3);
    platV << -platWidth, 0.0, -platWidth,
        -platWidth, 0.0, platWidth,
        platWidth, 0.0, platWidth,
        platWidth, 0.0, -platWidth,
        -platWidth, -platWidth / 10.0, -platWidth,
        -platWidth, -platWidth / 10.0, platWidth,
        platWidth, -platWidth / 10.0, platWidth,
        platWidth, -platWidth / 10.0, -platWidth;
    platF << 0, 1, 2,
        2, 3, 0,
        6, 5, 4,
        4, 7, 6,
        1, 0, 5,
        0, 4, 5,
        2, 1, 6,
        1, 5, 6,
        3, 2, 7,
        2, 6, 7,
        0, 3, 4,
        3, 7, 4;

    platOrientation << 1.0, 0.0, 0.0, 0.0;

    igl::copyleft::tetgen::tetrahedralize(platV, platF, "pq1.414Y", platV, platT, platF);

    MatrixXi tempF(platF.rows(), 3);
    tempF << platF.col(2), platF.col(1), platF.col(0);
    platF = tempF;
}

void updateMeshes(igl::opengl::glfw::Viewer &viewer)
{

    RowVector3d platColor;
    platColor << 0.8, 0.8, 0.8;
    RowVector3d meshColor;
    meshColor << 0.8, 0.2, 0.2;
    MatrixXd currV(scene.meshes[0].origPositions.size() / 3, 3);
    for (int j = 0; j < scene.meshes[0].origPositions.size() / 3; j++)
        currV.row(j) << scene.meshes[0].currPositions.segment(3 * j, 3).transpose();
    viewer.core.align_camera_center(currV);
    for (int i = 0; i < scene.meshes.size(); i++)
    {
        MatrixXd currV(scene.meshes[i].origPositions.size() / 3, 3);
        for (int j = 0; j < scene.meshes[i].origPositions.size() / 3; j++)
            currV.row(j) << scene.meshes[i].currPositions.segment(3 * j, 3).transpose();
        viewer.data_list[i].clear();
        viewer.data_list[i].set_mesh(currV, scene.meshes[i].F);
        viewer.data_list[i].set_face_based(true);
        viewer.data_list[i].set_colors(meshColor);
        viewer.data_list[i].show_lines = false;
    }
    viewer.data_list[0].show_lines = false;
    viewer.data_list[0].set_colors(platColor.replicate(scene.meshes[0].F.rows(), 1));
    viewer.data_list[0].set_face_based(true);
    //viewer.core.align_camera_center(scene.meshes[0].currV);

    //updating constraint viewing
    MatrixXi constF;
    MatrixXd constV, constC;

    MatrixXd P1(scene.userConstraints.size(), 3);
    MatrixXd P2(scene.userConstraints.size(), 3);
    for (int i = 0; i < scene.userConstraints.size(); i++)
    {
        cout << "scene.userConstraints[i].globalIndices:" << scene.userConstraints[i].globalIndices << endl;
        P1.row(i) = scene.globalPositions.segment(scene.userConstraints[i].globalIndices(0), 3).transpose();
        P2.row(i) = scene.globalPositions.segment(scene.userConstraints[i].globalIndices(3), 3).transpose();
    }

    MatrixXd cyndColors = RowVector3d(1.0, 1.0, 0.0).replicate(P1.size(), 1);

    double radius = 0.5;
    hedra::line_cylinders(P1, P2, radius, cyndColors, 8, constV, constF, constC);
    viewer.data_list[scene.meshes.size()].set_mesh(constV, constF);
    viewer.data_list[scene.meshes.size()].set_face_based(true);
    viewer.data_list[scene.meshes.size()].set_colors(constC);
    viewer.data_list[scene.meshes.size()].show_lines = false;
}

bool key_down(igl::opengl::glfw::Viewer &viewer, unsigned char key, int modifier)
{
    if (key == ' ')
    {
        viewer.core.is_animating = !viewer.core.is_animating;
        if (viewer.core.is_animating)
            cout << "Simulation running" << endl;
        else
            cout << "Simulation paused" << endl;
        return true;
    }

    if (key == 'S')
    {
        if (!viewer.core.is_animating)
        {
            scene.updateScene(timeStep, CRCoeff, tolerance, maxIterations);
            currTime += timeStep;
            updateMeshes(viewer);
            std::cout << "currTime: " << currTime << std::endl;
            return true;
        }
    }
    return false;
}

bool pre_draw(igl::opengl::glfw::Viewer &viewer)
{
    using namespace Eigen;
    using namespace std;

    if (viewer.core.is_animating)
    {
        scene.updateScene(timeStep, CRCoeff, tolerance, maxIterations);
        currTime += timeStep;
        //cout <<"currTime: "<<currTime<<endl;
        updateMeshes(viewer);
    }

    return false;
}

class CustomMenu : public igl::opengl::glfw::imgui::ImGuiMenu
{
	float optionalForce[3] = { 0, 0, 0 };
	float shrinkFactor = 11;
    virtual void draw_viewer_menu() override
    {
        // Draw parent menu
        ImGuiMenu::draw_viewer_menu();

        // Add new group
        if (ImGui::CollapsingHeader("Algorithm Options", ImGuiTreeNodeFlags_DefaultOpen))
        {
            ImGui::InputFloat("CR Coeff", &CRCoeff, 0, 0, 3);

			if (ImGui::InputFloat("Time Step", &timeStep))
			{
				mgpViewer.core.animation_max_fps = (((int)1.0 / timeStep));
			}
        }
		if (ImGui::CollapsingHeader("Interaction", ImGuiTreeNodeFlags_DefaultOpen))
		{
			//ImGui::InputFloat3("Force to squish", optionalForce, 2);
			ImGui::InputFloat("Shrink Factor", &shrinkFactor);
			if (ImGui::Button("Apply Shrink"))
			{
				scene.ShrinkMesh(shrinkFactor);
				//scene.setForceOptional((double)optionalForce[0], (double)optionalForce[1], (double)optionalForce[2]);
			}
		}
    }
};

int main(int argc, char *argv[])
{
    using namespace Eigen;
    using namespace std;

    // Load scene
    if (argc < 4)
    {
        cout << "Please provide path (argument 1), name of scene file (argument 2), and name of constraints file (argument 3)!" << endl;
        return 0;
    }
    cout << "scene file: " << std::string(argv[2]) << endl;
    //create platform
    createPlatform();

    scene.addMesh(platV, platF, platT, 10000.0, 0.5, 10000.0, true, platCOM, platOrientation);

    //load scene from file
    scene.loadScene(std::string(argv[1]), std::string(argv[2]), std::string(argv[3]));
    scene.initScene(timeStep, 0.1, 0.1);
    //scene.setPlatformBarriers(platV, CRCoeff);

    scene.updateScene(0.0, CRCoeff, tolerance, maxIterations);

    // Viewer Settings
    for (int i = 0; i < scene.meshes.size(); i++)
    {
        if (i != 0)
            mgpViewer.append_mesh();
        //mgpViewer.data_list[i].set_mesh(scene.meshes[i].currV, scene.meshes[i].F);
    }
    //mgpViewer.core.align_camera_center(scene.meshes[0].currV);

    //constraints mesh (for lines)
    mgpViewer.append_mesh();
    mgpViewer.callback_pre_draw = &pre_draw;
    mgpViewer.callback_key_down = &key_down;
    mgpViewer.core.is_animating = false;
    mgpViewer.core.animation_max_fps = 50.;
    updateMeshes(mgpViewer);
    CustomMenu menu;
    mgpViewer.plugins.push_back(&menu);

    cout << "Press [space] to toggle continuous simulation" << endl;
    cout << "Press 'S' to advance time step-by-step" << endl;

    mgpViewer.launch();
}
