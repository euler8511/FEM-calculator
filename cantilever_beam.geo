// Gmsh project created on Sun Sep 21 22:50:23 2025
//+
Point(1) = {0, 0, 0, 1};
//+
Point(2) = {2, 0, 0, 1};
//+
Line(1) = {1, 2};
//+
Physical Point("fix", 2) = {1, 2};
//+
Physical Point(" fix", 2) -= {2};
//+
Physical Point("load_y", 3) = {2};
//+
Physical Curve("beam", 4) = {1};
