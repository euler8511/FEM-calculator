// Gmsh project created on Sun Sep 21 19:37:26 2025
//+
Point(1) = {0, 0, 0, 1.0};
//+
Point(2) = {0, 1, 0, 1.0};
//+
Point(3) = {0.7, 1, 0, 1.0};
//+
Point(4) = {0.7, 0, 0, 1.0};
//+
Line(1) = {1, 2};
//+
Line(2) = {2, 3};
//+
Line(3) = {3, 4};
//+
Recursive Delete {
  Curve{2}; 
}
//+
Point(5) = {0.35, 1, 0, 1};
//+
Line(4) = {2, 5};
//+
Line(5) = {5, 3};
//+
Physical Point("fix", 6) = {1, 4};
//+
Physical Point("load_y", 7) = {5};
//+
Physical Curve("beam", 8) = {1, 4, 5, 3};
//+
Physical Curve("l_section", 9) = {1, 3};
//+
Physical Curve("c_section", 10) = {4, 5};
//+
Physical Curve(" beam", 8) -= {4, 5, 1, 3};
