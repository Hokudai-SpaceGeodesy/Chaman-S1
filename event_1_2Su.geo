//left-lateral for 20160513 and 20160710 on Chaman fault

// Gmsh project created on Sat Mar 26 12:42:40 2011
//it took a long time to make the residual from -8 to -6
//-3~2,0.6,0.35
//Strike 36

Point(1)={246850,3380900,-50,750};
Point(2)={255550,3395000,-50,750};
Point(6)={258950,3403000,-50,750};
//Point(7)={264200,3419500,-50,1500};
Point(9)={265250,3423600,-50,750};

Point(11)={246450,3381000,-20000,8000};
Point(12)={254850,3393500,-20000,8000};
Point(16)={258800,3403000,-20000,8000};
//Point(17)={264400,3419500,-20000,5000};
Point(19)={265150,3423600,-20000,8000};


Spline(1)={1,2,6,9};
Spline(2)={9,19};
Spline(3)={19,16,12,11};
Spline(4)={11,1};

//Spline(1)={1,2,6,7,9};
//Spline(2)={9,19};
//Spline(3)={19,17,16,12,11};
//Spline(4)={11,1};

Line Loop(5)={1,2,3,4};
Ruled Surface(6)={5};
Physical Surface(7)={6};

Coherence;
