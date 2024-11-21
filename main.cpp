#include <GLFW/glfw3.h>
#include <cmath>
#include <vector>
#include <iostream>

struct Point { // Define a struct for a 2D point
    float x; // X coordinate of the point
    float y; // Y coordinate of the point
};

struct Node : Point { // Define a struct for a node
    bool hasHandle1, hasHandle2; // Flags indicating whether the node has control handles/points
    Point handle1; // First control handle of the node
    Point handle2; // Second control handle of the node
};

std::vector<Node> nodes; // Vector to store all the nodes
Node* currentNode = nullptr; // Pointer to the currently selected node
Node* currentCP = nullptr; // Pointer to the currently selected control point

const float NODE_SIZE = 10.0f;
const float CONTROL_POINT_SIZE = 10.0f;
const float TOLERANCE = 10.0f;

// Function to calculate distance between two points
float distCalc(const Point& p1, const Point& p2) {
    float dx = p1.x - p2.x;
    float dy = p1.y - p2.y;
    return sqrt(dx * dx + dy * dy);
}

// Function to reset the program
void resetProgram() {
    nodes.clear();
    currentNode = nullptr;
    currentCP = nullptr;
}

// Function called when a key is pressed
void keyEPressed(GLFWwindow* window, int key, int scancode, int action, int mods) {
    // Check if the key pressed is "E" and if it's a key press event
    if (key == GLFW_KEY_E && action == GLFW_PRESS) {
        // If so, reset the program
        resetProgram();
    }
}

// Function called when a mouse button is pressed or released
void mouseLeftClick(GLFWwindow* window, int button, int action, int mods) {
    // Retrieve the framebuffer size to adjust mouse coordinates
    int W, H;
    glfwGetFramebufferSize(window, &W, &H);

    // Check if the left mouse button is pressed
    if (button == GLFW_MOUSE_BUTTON_LEFT && action == GLFW_PRESS) {
        double xcoor, ycoor;
        glfwGetCursorPos(window, &xcoor, &ycoor);
        ycoor = H - ycoor; // Invert y-axis to match OpenGL coordinate system

        bool existNode = false;
        bool existCP = false;

        // Check if clicked on an existing node or control point
        for (auto& n : nodes) {
            float dxNode = n.x - xcoor;
            float dyNode = n.y - ycoor;
            if (dxNode * dxNode + dyNode * dyNode < TOLERANCE * TOLERANCE) {
                currentNode = &n;
                existNode = true;
                break;
            }
            if (n.hasHandle1 || n.hasHandle2) {
                float dxCP1 = n.x + n.handle1.x - xcoor;
                float dyCP1 = n.y + n.handle1.y - ycoor;
                if (dxCP1 * dxCP1 + dyCP1 * dyCP1 < TOLERANCE * TOLERANCE) {
                    currentCP = &n;
                    existCP = true;
                    break;
                }
                if (n.hasHandle2) {
                    float dxCP2 = n.x + n.handle2.x - xcoor;
                    float dyCP2 = n.y + n.handle2.y - ycoor;
                    if (dxCP2 * dxCP2 + dyCP2 * dyCP2 < TOLERANCE * TOLERANCE) {
                        currentCP = &n;
                        existCP = true;
                        break;
                    }
                }
            }
        }

        // If not clicked on an existing node or control point, create a new node
        if (!existNode && !existCP) {
            Node newEndpoint;
            newEndpoint.x = xcoor;
            newEndpoint.y = ycoor;
            newEndpoint.hasHandle1 = true;
            newEndpoint.handle1.x = 0;
            newEndpoint.handle1.y = 50; // Set handle 1 to be 50 pixels above the node
            newEndpoint.hasHandle2 = false;

            // If there are no nodes or only 1 node, simply add the new node with its control point
            if (nodes.empty() || nodes.size() < 2) {
                nodes.push_back(newEndpoint);
            } else {
                Node& lastNode = nodes.back();
                float distToFirstNode = distCalc(newEndpoint, nodes.front());
                float distToLastNode = distCalc(newEndpoint, nodes.back());

                // Determine whether to insert the new node at the beginning or end
                if (distToFirstNode < distToLastNode) {
                    float currentHandle1Distance = sqrt(nodes.front().handle1.x * nodes.front().handle1.x + nodes.front().handle1.y * nodes.front().handle1.y);
                    float currentHandle1DirectionX = nodes.front().handle1.x / currentHandle1Distance;
                    float currentHandle1DirectionY = nodes.front().handle1.y / currentHandle1Distance;

                    nodes.front().hasHandle2 = true;
                    nodes.front().handle2.x = -currentHandle1DirectionX * currentHandle1Distance;
                    nodes.front().handle2.y = -currentHandle1DirectionY * currentHandle1Distance;
                    nodes.insert(nodes.begin(), newEndpoint);

                } else {
                    float currentHandle1Distance = sqrt(lastNode.handle1.x * lastNode.handle1.x + lastNode.handle1.y * lastNode.handle1.y);
                    float currentHandle1DirectionX = lastNode.handle1.x / currentHandle1Distance;
                    float currentHandle1DirectionY = lastNode.handle1.y / currentHandle1Distance;

                    lastNode.hasHandle2 = true;
                    lastNode.handle2.x = -currentHandle1DirectionX * currentHandle1Distance;
                    lastNode.handle2.y = -currentHandle1DirectionY * currentHandle1Distance;
                    nodes.push_back(newEndpoint);
                }
            }
        }

    } else if (button == GLFW_MOUSE_BUTTON_LEFT && action == GLFW_RELEASE) {
        // If the left mouse button is released, reset selected node and control point
        currentNode = nullptr;
        currentCP = nullptr;
    }
}

// Function called when the cursor position changes
void cursorMoving(GLFWwindow* window, double xcoor, double ycoor) {
    // Retrieve the framebuffer size to adjust mouse coordinates
    int W, H;
    glfwGetFramebufferSize(window, &W, &H);

    if (currentNode) {
        ycoor = H - ycoor; // Adjust y-axis to match OpenGL coordinate system
        currentNode->x = xcoor;
        currentNode->y = ycoor;

        if (currentNode->hasHandle1 && currentNode->hasHandle2) {
            float dx = xcoor - currentNode->x;
            float dy = ycoor - currentNode->y;
            currentNode->handle1.x += dx;
            currentNode->handle1.y += dy;
            currentNode->handle2.x -= dx;
            currentNode->handle2.y -= dy;
        }
    }
    if (currentCP) {
        ycoor = H - ycoor; // Invert y-axis
        if (currentCP->hasHandle1) {
            currentCP->handle1.x = xcoor - currentCP->x;
            currentCP->handle1.y = ycoor - currentCP->y;

            if (currentCP->hasHandle2) {
                float dx = currentCP->handle1.x;
                float dy = currentCP->handle1.y;
                currentCP->handle2.x = -dx;
                currentCP->handle2.y = -dy;
            }

        }
        if (currentCP->hasHandle2) {
            currentCP->handle2.x = xcoor - currentCP->x;
            currentCP->handle2.y = ycoor - currentCP->y;

            if (currentCP->hasHandle1) {
                float dx = currentCP->handle2.x;
                float dy = currentCP->handle2.y;
                currentCP->handle1.x = -dx;
                currentCP->handle1.y = -dy;
            }
        }
    }
}

// Function to render
void render() {
    glClear(GL_COLOR_BUFFER_BIT); // Clear the color buffer

    // Render spline
    glEnable(GL_LINE_SMOOTH);
    glEnable(GL_BLEND);
    glBlendFunc(GL_SRC_ALPHA, GL_ONE_MINUS_SRC_ALPHA);
    glColor3f(1.0f, 1.0f, 1.0f);
    glLineWidth(3.0f);
    glBegin(GL_LINE_STRIP);
    for (size_t i = 0; i < nodes.size(); ++i) { // Iterate over all nodes
        const Node& n = nodes[i]; // Get the current node
        if (i == 0) { // If it's the first node
            glVertex2f(n.x, n.y); // Draw the vertex
        }else if (i==1){ // If it's the second node
            const Node& prev = nodes[0];
            const float segments = 200.0f; // Number of segments for bezier curve
                for (float t = 0; t <= 1; t += 1 / segments) { // Iterate over the curve
                    float x = (1 - t) * (1 - t) * (1 - t) * prev.x + 3 * (1 - t) * (1 - t) * t * (prev.x + prev.handle1.x) + 3 * (1 - t) * t * t * (n.x + n.handle1.x) + t * t * t * n.x;
                    float y = (1 - t) * (1 - t) * (1 - t) * prev.y + 3 * (1 - t) * (1 - t) * t * (prev.y + prev.handle1.y) + 3 * (1 - t) * t * t * (n.y + n.handle1.y) + t * t * t * n.y;
                    glVertex2f(x, y);
                }

        } else { // For intermediate nodes
            const Node& prev = nodes[i - 1];
            // Cubic Bezier curve calculation
            const float segments = 200.0f; // Number of segments for bezier curve
            for (float t = 0; t <= 1; t += 1 / segments) { // Iterate over the curve
                float x = (1 - t) * (1 - t) * (1 - t) * prev.x + 3 * (1 - t) * (1 - t) * t * (prev.x + prev.handle2.x) + 3 * (1 - t) * t * t * (n.x + n.handle1.x) + t * t * t * n.x;
                float y = (1 - t) * (1 - t) * (1 - t) * prev.y + 3 * (1 - t) * (1 - t) * t * (prev.y + prev.handle2.y) + 3 * (1 - t) * t * t * (n.y + n.handle1.y) + t * t * t * n.y;
                glVertex2f(x, y);

            }
        }
    }
    glEnd();

    // Render nodes
    glColor3f(0.0f, 0.0f, 1.0f); // Set color to blue
    glPointSize(NODE_SIZE);
    //glBegin(GL_POINTS);
    for (const auto& n : nodes) {
        glBegin(GL_QUADS); // Begin drawing quads
        float halfNode = NODE_SIZE / 2.0f;
        glVertex2f(n.x - halfNode, n.y - halfNode); // Bottom-left corner
        glVertex2f(n.x + halfNode, n.y - halfNode); // Bottom-right corner
        glVertex2f(n.x + halfNode, n.y + halfNode); // Top-right corner
        glVertex2f(n.x - halfNode, n.y + halfNode); // Top-left corner
        glEnd();
    }
    //glEnd();

    // Render control points
    glColor3f(1.0f, 1.0f, 1.0f); // White color
    glEnable(GL_POINT_SMOOTH);
    glEnable(GL_BLEND);
    glBlendFunc(GL_SRC_ALPHA, GL_ONE_MINUS_SRC_ALPHA);
    glPointSize(CONTROL_POINT_SIZE);
    glBegin(GL_POINTS);
    for (const auto& n : nodes) {
        if (n.hasHandle1) { // If node has handle 1
            glVertex2f(n.x + n.handle1.x, n.y + n.handle1.y); // Draw control point 1
        }
        if (n.hasHandle2) { // If node has handle 2
            glVertex2f(n.x + n.handle2.x, n.y + n.handle2.y); // Draw control point 2
        }
    }
    glEnd();

    // Render dotted lines connecting control points to nodes
    glEnable(GL_LINE_STIPPLE);
    glLineStipple(1, 0xAAAA);
    glColor3f(0.6f, 0.8f, 0.9f); // Light blue color
    glBegin(GL_LINES);
    for (const auto& n : nodes) {
        if (n.hasHandle1) { // If node has handle 1
            // Draw line from node to control point 1
            glVertex2f(n.x, n.y);
            glVertex2f(n.x + n.handle1.x, n.y + n.handle1.y);
        }
        if (n.hasHandle2) { // If node has handle 2
            // Draw line from node to control point 2
            glVertex2f(n.x, n.y);
            glVertex2f(n.x + n.handle2.x, n.y + n.handle2.y);
        }
    }
    glEnd();

    glDisable(GL_LINE_STIPPLE);
}

// Main function
int main(int argc, char* argv[]) {
    // Check if screen width and height arguments are provided
    if (argc < 3) {
        std::cerr << "Usage: " << argv[0] << " <screen_width> <screen_height>" << std::endl;
        return -1;
    }

    // Parse screen width and height from command-line arguments
    int W = std::atoi(argv[1]);
    int H = std::atoi(argv[2]);

    // Initialize GLFW
    if (!glfwInit()) {
        return -1;
    }


    // Set GLFW window hints; Enable 4 times multisampling
    glfwWindowHint(GLFW_SAMPLES, 4);

    // Create a windowed mode window and its OpenGL context
    GLFWwindow* window = glfwCreateWindow(W, H, "Cubic Bezier Splines", NULL, NULL);
    if (!window) {
        glfwTerminate();
        return -1;
    }

    // Make the window's context current
    glfwMakeContextCurrent(window);

    // Set up viewport
    glViewport(0, 0, W, H);
    glMatrixMode(GL_PROJECTION);
    glLoadIdentity();
    glOrtho(0.0f, W, 0.0f, H, -1.0f, 1.0f);
    glMatrixMode(GL_MODELVIEW);

    glClearColor(0.0f, 0.0f, 0.0f, 1.0f);

    // Set GLFW callback functions
    glfwSetMouseButtonCallback(window, mouseLeftClick);
    glfwSetCursorPosCallback(window, cursorMoving);
    glfwSetKeyCallback(window, keyEPressed);

    // Loop until the user closes the window
    while (!glfwWindowShouldClose(window)) {
        // Render here
        render();

        // Swap and back buffers
        glfwSwapBuffers(window);

        // Poll for and process events
        glfwPollEvents();
    }

    glfwTerminate();
    return 0;
}