
Include "plate.geo";

// Remove previous surfaces
Delete Physicals;
Delete {
  Surface {1};
  Line{4};
}

Point(5) = {a,0.0,0.0,lc2};
Point(6) = {2*a,0.0,0.0,lc1};

Line(4) = {4,5};
Line(5) = {5,6};
Line(6) = {6,1};

Line Loop(13) = {1,2,3,4,5,6};
Plane Surface(1) = {13};
Physical Surface(1) = {1};
