////////////////////////////////////////////////////////////
// GEOMETRY FOR THE IRON YOKE OF SUPERCONDUCTING BENDING  //
// MAGNET USED IN THE HEAVY ION SYNCHROTRON.              //
//                                                        //
// AUTHOR: Ye Yang (QST)                                  //
// DATE  : 2020.05.23                                     //
////////////////////////////////////////////////////////////
SetFactory("OpenCASCADE");

// UNITS
mm = 1e-3;
cm = 1e-2;
m  = 1.00;
deg= Pi/180.0;

// SIZE
size1 = 2.2*mm;
size2 = 150*mm;
size3 = 1.5*mm;

// CONSTANTS
Rc     = 1.890*m;
Xbound = 550*mm;
Ybound = 450*mm;

// PARAMETERS
Xin    = 125*mm;
Yin    = 109*mm;
Xout1  = 80*mm; 
Xout2  = 310*mm;
Yout1  = 85*mm;
Yout2  = 230*mm;
Ahole  = 200*mm;
Rcirc1 = 5*mm;
Rcirc2 = 5*mm;
Rcirc3 = 5*mm;
Theta1 = 10*deg; 
Theta2 = 40*deg; 
Theta3 = 60*deg; 

// ELLIPTIC PARAMETERS
Bhole = Yin/Xin*Ahole;
Ehole = Sqrt(Ahole*Ahole - Bhole*Bhole);
Eta   = 0.5*Log((1+Bhole/Ahole)/(1-Bhole/Ahole));

// SETUP POINTS 
Point( 1) = {        Rc,    0.0, 0.0, size3};
Point( 2) = {-Xbound+Rc,    0.0, 0.0, size2};
Point( 3) = { Xbound+Rc,    0.0, 0.0, size2};
Point( 4) = {-Xbound+Rc, Ybound, 0.0, size2};
Point( 5) = { Xbound+Rc, Ybound, 0.0, size2};
Point( 6) = {        Rc,    Yin, 0.0, size1};
Point( 7) = {   -Xin+Rc,    0.0, 0.0, size1};
Point( 8) = {    Xin+Rc,    0.0, 0.0, size1};
Point( 9) = { -Xout2+Rc,    0.0, 0.0, size1};
Point(10) = {  Xout2+Rc,    0.0, 0.0, size1};
Point(11) = { -Xout2+Rc,  Yout1, 0.0, size1};
Point(12) = {  Xout2+Rc,  Yout1, 0.0, size1};
Point(13) = { -Xout1+Rc,  Yout2, 0.0, size1};
Point(14) = {  Xout1+Rc,  Yout2, 0.0, size1};

// SETUP HOLE FOR FIELD ADJUSTMENT
Circle(1) = {Rc+Ehole*Cosh(Eta)*Cos(Theta1), Ehole*Sinh(Eta)*Sin(Theta1), 0, Rcirc1, 0, 2*Pi};
Circle(2) = {Rc+Ehole*Cosh(Eta)*Cos(Theta2), Ehole*Sinh(Eta)*Sin(Theta2), 0, Rcirc2, 0, 2*Pi};
Circle(3) = {Rc+Ehole*Cosh(Eta)*Cos(Theta3), Ehole*Sinh(Eta)*Sin(Theta3), 0, Rcirc3, 0, 2*Pi};
Circle(4) = {Rc-Ehole*Cosh(Eta)*Cos(Theta1), Ehole*Sinh(Eta)*Sin(Theta1), 0, Rcirc1, 0, 2*Pi};
Circle(5) = {Rc-Ehole*Cosh(Eta)*Cos(Theta2), Ehole*Sinh(Eta)*Sin(Theta2), 0, Rcirc2, 0, 2*Pi};
Circle(6) = {Rc-Ehole*Cosh(Eta)*Cos(Theta3), Ehole*Sinh(Eta)*Sin(Theta3), 0, Rcirc3, 0, 2*Pi};

// ELLIPTICAL ARC
Ellipse(7) = {8, 1, 7, 6};
Ellipse(8) = {7, 1, 8, 6};
Line( 9) = {7, 1};
Line(10) = {1, 8};

// IRON YOKE
Line(11) = {8, 10};
Line(12) = {10, 12};
Line(13) = {12, 14};
Line(14) = {14, 13};
Line(15) = {13, 11};
Line(16) = {11, 9};
Line(17) = {9, 7};

// SIMULATION BOUNDARY
Line(18) = {10, 3};
Line(19) = {3, 5};
Line(20) = {5, 4};
Line(21) = {4, 2};
Line(22) = {2, 9};

// SURFACE::ELLIPTICAL ARC::AIR
Curve Loop(1) = {8, -7, -10, -9};
Surface(1) = {1};
// SURFACE::HOLE::AIR
Curve Loop(3) = {4};
Surface(2) = {3};
Curve Loop(5) = {5};
Surface(3) = {5};
Curve Loop(7) = {6};
Surface(4) = {7};
Curve Loop(9) = {3};
Surface(5) = {9};
Curve Loop(11) = {2};
Surface(6) = {11};
Curve Loop(13) = {1};
Surface(7) = {13};
// SURFACE::IRON YOKE::IRON
Curve Loop(19) = {8, -7, 11, 12, 13, 14, 15, 16, 17};
Curve Loop(20) = {1};
Curve Loop(21) = {2};
Curve Loop(22) = {3};
Curve Loop(23) = {6};
Curve Loop(24) = {5};
Curve Loop(25) = {4};
Plane Surface(10) = {19, 20, 21, 22, 23, 24, 25};
// SURFACE::AIR::AIR
Curve Loop(17) = {15, 16, -22, -21, -20, -19, -18, 12, 13, 14};
Surface(9) = {17};

// BOUNDARY AND INTERFACE CONDITIONS
Physical Curve("BC_NORM"    , 23) = {22, 17, 9, 10, 11, 18};
Physical Curve("BC_TANG"    , 24) = {21, 20, 19};
Physical Curve("IC_ARC_IN"  , 25) = {8, 7};
Physical Curve("IC_YOKE_OUT", 26) = {16, 15, 14, 13, 12};
Physical Curve("IC_HOLES"   , 27) = {6, 5, 4, 3, 2, 1};

// SURFACE OF DOMAINS 
Physical Surface("SURF_AIR_OUT", 28) = {9};
Physical Surface("SURF_YOKE"   , 29) = {10};
Physical Surface("SURF_HOLES"  , 30) = {4, 3, 2, 7, 6, 5};
Physical Surface("SURF_AIR_IN" , 31) = {1};

// MESHING
// HOLE
Transfinite Curve {5, 6, 4, 3, 2, 1} = 14 Using Progression 1;
// AIR INNER BORE
Transfinite Curve {8, 7 } = 75 Using Progression 1;
Transfinite Curve {9, 10} = 70 Using Progression 1;
// IRON YOKE
Transfinite Curve {16, 12} = 50 Using Progression 1;
Transfinite Curve {15, 13} = 40 Using Progression 1;
Transfinite Curve {17, 11, 14} = 40 Using Progression 1;

// Meshing
Smoother Surface {2,3,4,5,6,7} = 15;
Mesh 1;
Mesh 2;
SetOrder 1;

// write out gmsh file
Save "opt_SCBM_IronYokeShape.msh";

