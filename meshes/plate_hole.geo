Include "plate.geo";

// Remove previous surfaces
Delete Physicals;
Delete {
  Surface {1};
  Line {3,4};
}

// Domain construction
Point(5) = {0.0,L,0.0,lc1};
Point(6) = {0.0,a,0.0,lc2};
Point(7) = {a,0.0,0.0,lc2};

Line(3) = {3,5};
Line(4) = {5,6};
Circle(5) = {6,4,7};
Line(6) = {7,1};

Line Loop(12) = {1,2,3,4,5,6};
Plane Surface(1) = {12};
Physical Surface(1) = {1};
