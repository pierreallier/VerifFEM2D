
// Geometric parameters
H = 10;  //  semiheight of plate
L = H/2;  //  semiwidth of plate
a = 1; // semilength of crack (center of crack is at x1=x2=0)

// Discretization parameters
lc1 = 5; // element size at the border
lc2 = .5; // element size at the crack tip

// Domain construction
Point(1) = {L,0.0,0.0,lc1};
Point(2) = {L,H,0.0,lc1};
Point(3) = {0.0,H,0.0,lc1};
Point(4) = {0.0,0.0,0.0,lc1};

Line(1) = {1,2};
Line(2) = {2,3};
Line(3) = {3,4};
Line(4) = {4,1};

Line Loop(11) = {1,2,3,4};
Plane Surface(1) = {11};

// To return only 2D element in msh file
Physical Surface(1) = {1};
