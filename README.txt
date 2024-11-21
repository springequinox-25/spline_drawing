Cubic Bezier Splines Drawing Program

DESCRIPTION

This program allows you to interactively create and manipulate cubic Bezier splines using OpenGL and GLFW.
The program provides a graphical interface where you can create and manipulate nodes that define the shape of the cubic Bezier spline. Each node can have control handles that affect the curvature of the spline segments connected to it.


FEATURES

- Interactive creation and manipulation of cubic Bezier splines.
- Control handles for adjusting the curvature of spline segments.
- Reset functionality to clear all nodes and start anew.


SPECIAL LIBRARIES

This program relies on the following libraries (other than general c++ libraries):

- GLFW (OpenGL Framework)
- OpenGL


COMPILATION and RUN

- The program has been successfully run on CodeBlocks. (*The command line below is not being tested on the terminal)
- To run this program, ensure you have the necessary libraries installed on the system. Then, compile the program using a C++ compiler. For example, you can compile it using g++:
	g++ -o cubic_bezier_splines main.cpp -lOpenGl -lglfw

- After compiling the program, execute the binary providing the screen width and height as command-line arguments. For example:
	./cubic_bezier_splines ScreenWidth ScreenHeight
	./cubic_bezier_splines 800 600

- This will open a window where you can interactively create and manipulate cubic Bezier splines.


CONTROLS/INTERACTIONS

- Left Mouse Click (Press): Create or select a node or control point.
- Left Mouse Click (Release): Deselect the current node or control point.
- Mouse Movement: Move the selected node or control point.
- Key E: Reset the program, clearing all nodes.