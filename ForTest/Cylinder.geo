// Gmsh project created on Thu Aug 10 16:31:53 2023
SetFactory("OpenCASCADE");
//+
Mesh.MshFileVersion = 2;
//+
Circle(1) = {0, 0, 0, 1, 0, 2*Pi};
//+
Extrude {0, 0, 2} {
  Curve{1}; 
}
//+
Physical Curve("inlet", 1) = {3};
//+
Physical Curve("outlet", 2) = {1};
//+
Curve Loop(2) = {3};
//+
Surface(2) = {2};
//+
Curve Loop(4) = {1};
//+
Surface(3) = {4};
//+
Surface Loop(1) = {2, 1, 3};
//+
Volume(1) = {1};
//+
Physical Surface("surf", 3) = {1};
//+
Physical Volume("vol", 4) = {1};
