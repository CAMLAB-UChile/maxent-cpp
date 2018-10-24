# maxent-cpp
A C++ implementation of the maximum-entropy basis functions

# Author
<a href="https://github.com/rsilvv">Rodrigo Silva-Valenzuela</a>, M.Sc. student, Department of Mechanical Engineering, Universidad de Chile.

# Supervisor
<a href="https://github.com/aaortizb">Alejandro Ortiz-Bernardin</a>, Assistant Professor, Department of Mechanical Engineering, Universidad de Chile.

# Instructions
This program is controlled by the "main.cpp" file located in the foler "src". 
Use "main.cpp" to define a test example. The default "main.cpp" file defines a 
a 3D cloud of 8 nodes and solves the maxent approximation at a given sampling 
point located in this cloud. Other examples are available in the same folder 
(see files "3nodes1d.cpp", "4nodes2d.cpp", "6nodes3d.cpp").

This program relies on g++ compiler and eigen library.

# Compilation
- On the top directory of maxent-cpp create a folder named "bin"
- Open a terminal and, at the top directory of maxent-cpp, type make to build an
  executable named "maxent-cpp" that will appear in the folder "bin"
- Go inside the folder "bin" and execute "maxent-cpp"

# License
This project is licensed under the GPL3 License. This program is free software; it can be redistributed or modified under the terms of the GNU General Public License 3 as published by the Free Software Foundation. 
