#include<iostream>
#include<algorithm>
#include<chrono>

#include<Eigen/Dense>
#include<Eigen/Geometry>

#include <GLM/glm.hpp>
#include <GL/glew.h>
#include <GLFW/glfw3.h>
#include "opengl_util/Camera.h"
#include "opengl_util/shader.h"

#include "BVHACCEL.h"
#include "ModelIO.h"
#include "KDACCEL.h"

GLFWwindow* window = NULL;
Shader* companionShader = NULL, *gBoundsShader = NULL, *gDrawingWithPhoneShader = NULL;
// camera
glm::vec3 cameraPos = glm::vec3(0, 0.1, 0.5);
glm::vec3 cameraFront = glm::vec3(0.0f, 0.0f, -1.0f);
glm::vec3 cameraUp = glm::vec3(0.0f, 1.0f, 0.0f);
Camera camera = Camera(cameraPos);
glm::mat4 leftViewMatrix;
glm::mat4 leftProjectionMatrix;
glm::mat4 rightViewMatrix;
glm::mat4 rightProjectionMatrix;
glm::mat4 leftHMDMatrix;
glm::mat4 rightHMDMatrix;
bool firstMouse = true;
float lastX = 512.0f / 2.0;
float lastY = 512.0f / 2.0;
// timing
float deltaTime = 0.0f; // time between current frame and last frame
float lastFrame = 0.0f;

GLuint VAO, VBO, IBO;
GLuint BoxVAO, BoxVBO;

int boxNum = 1;

std::vector<shared_ptr<Mesh>> meshes;
vector<float> vertexes;

void mouse_callback(GLFWwindow* window, double xpos, double ypos) {
	if (firstMouse) {
		lastX = xpos;
		lastY = ypos;
		firstMouse = false;
	}

	float xoffset = xpos - lastX;
	float yoffset = lastY - ypos; // reversed since y-coordinates go from bottom to top

	lastX = xpos;
	lastY = ypos;

	camera.ProcessMouseMovement(xoffset, yoffset);
}

void scroll_callback(GLFWwindow* window, double xoffset, double yoffset) {
	camera.ProcessMouseScroll(yoffset);
}

bool setUpDll() {
	// glfw: initialize and configure
	// ------------------------------
	glfwInit();
	glfwWindowHint(GLFW_CONTEXT_VERSION_MAJOR, 4);
	glfwWindowHint(GLFW_CONTEXT_VERSION_MINOR, 3);
	glfwWindowHint(GLFW_OPENGL_PROFILE, GLFW_OPENGL_CORE_PROFILE);
	glfwWindowHint(GLFW_OPENGL_DEBUG_CONTEXT, GL_TRUE);

	// glfw window creation
	// --------------------
	window = glfwCreateWindow(512, 512, "Test Accel", NULL, NULL);
	if (window == NULL) {
		std::cout << "Failed to create GLFW window" << std::endl;
		glfwTerminate();
		return false;
	}
	glfwMakeContextCurrent(window);
	glfwSetCursorPosCallback(window, mouse_callback);
	glfwSetScrollCallback(window, scroll_callback);

	// tell GLFW to capture our mouse
	glfwSetInputMode(window, GLFW_CURSOR, GLFW_CURSOR_DISABLED);

	//init glew
	glewExperimental = GL_TRUE;
	GLenum nGlewError = glewInit();
	if (nGlewError != GLEW_OK) {
		printf("%s - Error initializing GLEW! %s\n", __FUNCTION__, glewGetErrorString(nGlewError));
		return false;
	}

	return true;
}

void setupShaders() {
	companionShader = new Shader("C:/Users/whatseven/Desktop/libstructure/shaders/simple.vs"
		, "C:/Users/whatseven/Desktop/libstructure/shaders/simple.fs");
}

bool setUpBuffer() {
	//vector<Eigen::Vector3f> t(20);
	//vertexes.resize(meshes.size()*20*3);
	//size_t meshIndex = 0;
	//for(const auto &mesh:meshes){
	//	//for (auto point : mesh->points) {
	//	//	vertexes.push_back(point[0]);
	//	//	vertexes.push_back(point[1]);
	//	//	vertexes.push_back(point[2]);
	//	//}
	//	t[0]=(mesh->getBounds().pMax);
	//	t[1]=(Eigen::Vector3f(mesh->getBounds().pMin[0], t[0][1], t[0][2]));
	//	t[2]=(mesh->getBounds().pMax);
	//	t[3]=(Eigen::Vector3f(t[0][0], mesh->getBounds().pMin[1], t[0][2]));
	//	t[4]=(mesh->getBounds().pMax);
	//	t[5]=(Eigen::Vector3f(t[0][0], t[0][1], mesh->getBounds().pMin[2]));
	//	t[6]=(mesh->getBounds().pMin);
	//	t[7]=(Eigen::Vector3f(t[6][0], t[6][1], t[0][2]));
	//	t[8]=(mesh->getBounds().pMin);
	//	t[9]=(Eigen::Vector3f(t[0][0], t[6][1], t[6][2]));
	//	t[10]=(mesh->getBounds().pMin);
	//	t[11]=(Eigen::Vector3f(t[6][0], t[0][1], t[6][2]));
	//	t[12]=(Eigen::Vector3f(t[6][0], t[0][1], t[6][2]));
	//	t[13]=(Eigen::Vector3f(t[0][0], t[0][1], t[6][2]));
	//	t[14]=(Eigen::Vector3f(t[6][0], t[0][1], t[6][2]));
	//	t[15]=(Eigen::Vector3f(t[6][0], t[0][1], t[0][2]));
	//	t[16]=(Eigen::Vector3f(t[0][0], t[6][1], t[0][2]));
	//	t[17]=(Eigen::Vector3f(t[0][0], t[6][1], t[6][2]));
	//	t[18]=(Eigen::Vector3f(t[0][0], t[6][1], t[0][2]));
	//	t[19]=(Eigen::Vector3f(t[6][0], t[6][1], t[0][2]));
	//	size_t pointIndex = 0;
	//	for (const auto &point : t) {
	//		vertexes[meshIndex * 20 * 3 + pointIndex * 3 + 0] = (point[0]);
	//		vertexes[meshIndex * 20 * 3 + pointIndex * 3 + 1]=(point[1]);
	//		vertexes[meshIndex * 20 * 3 + pointIndex * 3 + 2]=(point[2]);
	//		++pointIndex;
	//	}
	//	meshIndex++;
	//}

	unsigned int indices[] = {
		 0, 1,0,3,1,2,2,3
		,4,5,5,6,6,7,4,7
		,0,6,1,5,2,4,3,7
	};
	vector<unsigned int> indicesVector(vertexes.size());
	for(int i=0;i< indicesVector.size()/24;++i){
		for(int j=0;j<24;++j){
			indicesVector[i * 24 + j] = indices[j]+8*i;
		}
	}
	glGenVertexArrays(1, &VAO);
	glGenBuffers(1, &VBO);
	glGenBuffers(1, &IBO);
	glBindVertexArray(VAO);
	glBindBuffer(GL_ARRAY_BUFFER, VBO);
	glBufferData(GL_ARRAY_BUFFER, vertexes.size()*sizeof(float), &vertexes[0], GL_STATIC_DRAW);
	glBindBuffer(GL_ELEMENT_ARRAY_BUFFER, IBO);
	glBufferData(GL_ELEMENT_ARRAY_BUFFER, indicesVector.size() * sizeof(float), &indicesVector[0], GL_STATIC_DRAW);

	glEnableVertexAttribArray(0);
	glVertexAttribPointer(0, 3, GL_FLOAT, GL_FALSE, 3 * sizeof(float), (void*)0);
	glBindVertexArray(0);
	glBindBuffer(GL_ARRAY_BUFFER, 0);

	if (VAO == 0 || VBO == 0)
		return false;

	//glGenVertexArrays(1, &BoxVAO);
	//glGenBuffers(1, &BoxVBO);
	//glBindVertexArray(BoxVAO);
	//glBindBuffer(GL_ARRAY_BUFFER, BoxVBO);
	//glBufferData(GL_ARRAY_BUFFER, sizeof(companionWindowVertices), &companionWindowVertices[0], GL_STATIC_DRAW);

	//glEnableVertexAttribArray(0);
	//glVertexAttribPointer(0, 2, GL_FLOAT, GL_FALSE, 4 * sizeof(float), (void*)0);
	//glBindVertexArray(0);
	//glBindBuffer(GL_ARRAY_BUFFER, 0);

	//if (BoxVAO == 0 || BoxVBO == 0)
	//	return false;

	return true;
}

void handleInput(GLFWwindow* window) {
	if (glfwGetKey(window, GLFW_KEY_ESCAPE) == GLFW_PRESS)
		glfwSetWindowShouldClose(window, true);
	if (glfwGetKey(window, GLFW_KEY_W) == GLFW_PRESS)
		camera.ProcessKeyboard(FORWARD, deltaTime);
	if (glfwGetKey(window, GLFW_KEY_S) == GLFW_PRESS)
		camera.ProcessKeyboard(BACKWARD, deltaTime);
	if (glfwGetKey(window, GLFW_KEY_A) == GLFW_PRESS)
		camera.ProcessKeyboard(LEFT, deltaTime);
	if (glfwGetKey(window, GLFW_KEY_D) == GLFW_PRESS)
		camera.ProcessKeyboard(RIGHT, deltaTime);
	if (glfwGetKey(window, GLFW_KEY_EQUAL) == GLFW_PRESS){
		boxNum += 1;
		cout << "Now box num:" << boxNum << endl;
	}
	if (glfwGetKey(window, GLFW_KEY_MINUS) == GLFW_PRESS) {
		boxNum -= 1;
		cout << "Now box num:" << boxNum << endl;
	}

}

void Render() {
	while (!glfwWindowShouldClose(window)) {
		glClearColor(0.5f, 0.5f, 0.5f, 1.0f);
		glClear(GL_COLOR_BUFFER_BIT | GL_DEPTH_BUFFER_BIT);

		float currentFrame = glfwGetTime();
		deltaTime = currentFrame - lastFrame;
		lastFrame = currentFrame;
		handleInput(window);

		companionShader->use();
		glm::mat4 projection = glm::perspective(glm::radians(camera.Zoom), 1.0f, 0.1f, 10000.0f);
		glm::mat4 view = camera.GetViewMatrix();
		companionShader->setMat4("projection", projection);
		companionShader->setMat4("view", view);

		glBindVertexArray(VAO);
		glDrawElements(GL_LINES, boxNum*24, GL_UNSIGNED_INT, 0);

		glBindVertexArray(0);

		glfwSwapBuffers(window);
		glfwPollEvents();
	}
}

int main(){
	std::vector<Eigen::Vector3f> vertices, normals;
	std::vector<int> faces, texcoords;

	tie(vertices, normals, texcoords, faces)= loadPLY(
		"c:/Users/whatseven/Desktop/bencemark/NY-1/ProxyMesh/gesox-mesh_orig.ply"
		, true);
	 
	meshes= std::vector<shared_ptr<Mesh>>(faces.size()/3);

	auto startTimer = chrono::system_clock::now();
	std::vector<Eigen::Vector3f> meshPoints(3);
	for (int i = 0; i < meshes.size();++i) {
		meshPoints[0]=(vertices[faces[i * 3 + 0]]);
		meshPoints[1] = (vertices[faces[i * 3 + 1]]);
		meshPoints[2] = (vertices[faces[i * 3 + 2]]);
		meshes[i] = make_shared<Mesh>(Mesh(meshPoints));
	}
	auto endTimer = chrono::system_clock::now();
	auto duration = chrono::duration_cast<chrono::microseconds>(endTimer - startTimer);

	cout << "Prepare Construct" << endl;
	//getchar();
	BVHAccel bvhTree(meshes);
	KdTreeAccel kdAccel(meshes);

	LinearBVHNode* bNode = bvhTree.getLinearNodes();
	int totalNode = bvhTree.totalLinearNodes;
	int bvhIndex = 0;
	vertexes.resize(totalNode * 8 * 3);
	size_t vertexIndex = 0;
	while (bvhIndex< totalNode){
		for (auto val : bNode[bvhIndex].boundVertices)
			vertexes[vertexIndex++] = val;
		++bvhIndex;
	}


	if (!setUpDll()) {
		std::cout << "init dll error" << std::endl;
		return 0;
	}
	setUpBuffer();
	setupShaders();

	Render();

	//getchar();
    return 0;
}