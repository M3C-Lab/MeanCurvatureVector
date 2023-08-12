// Gmsh project created on Fri Apr 14 19:47:47 2023
SetFactory("OpenCASCADE");
//+
Mesh.MshFileVersion = 2;
//+
Sphere(1) = {0, 0, 0, 4, -Pi/2, Pi/2, 2*Pi};
//+
Physical Surface("surf") = {1};
//+
Physical Volume("vol") = {1};
