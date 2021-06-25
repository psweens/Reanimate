//
// Created by Paul Sweeney on 14/02/2021.
//

#include "MicroCell.hpp"

using namespace reanimate;

void MicroCell::hexCell2D() {

    nnod = 12;
    nnodbc = 6;
    nseg = 12;

    setup_microCellArrays();
    for (int iseg = 0; iseg < nseg; iseg++)   {segname(iseg) = iseg;}
    for (int inod = 0; inod < nnod; inod++)   {nodname(inod) = inod;}

    double l_mean = mean(lengthDistrib);
    double l_SD = stddev(lengthDistrib);
    int n = 1;
    double Lm = sum(lognormal(l_mean, l_SD, n));

    double L1 = Lm/(2*sin(M_PI/6));
    double L2 = Lm/(2*tan(M_PI/6));
    double L3 = Lm*cos(M_PI/3);
    double L4 = Lm*sin(M_PI/3);

    double cap_mean = mean(diamDistrib);
    double cap_SD = stddev(diamDistrib);

    // Randomly assign diameters based on lognormal distribution
    for (int iseg = 0; iseg < nseg; iseg++)   {diam(iseg) = sum(lognormal(cap_mean, cap_SD, n));}

    // Along horizontal
    cnode(0,0) = L1;
    cnode(0,1) = 3*L1/2;
    cnode(0,2) = -L1;
    cnode(0,3) = -3*L1/2;

    // y > 0
    cnode(0,4) = L1/2;
    cnode(1,4) = L2;
    cnode(0,5) = cnode(0,4) + L3/2;
    cnode(1,5) = cnode(1,4) + L4/2;
    cnode(0,6) = -L1/2;
    cnode(1,6) = L2;
    cnode(0,7) = -(cnode(0,4) + L3/2);
    cnode(1,7) = cnode(1,4) + L4/2;


    // y < 0
    cnode(0,8) = L1/2;
    cnode(1,8) = -L2;
    cnode(0,9) = cnode(0,4) + L3*0.5;
    cnode(1,9) = -(cnode(1,4) + L4*0.5);
    cnode(0,10) = -L1/2;
    cnode(1,10) = -L2;
    cnode(0,11) = -(cnode(0,4) + L3*0.5);
    cnode(1,11) = -(cnode(1,4) + L4*0.5);

    ista(0) = 0;
    iend(0) = 1;
    ista(1) = 0;
    iend(1) = 4;
    ista(2) = 4;
    iend(2) = 5;
    ista(3) = 6;
    iend(3) = 4;
    ista(4) = 6;
    iend(4) = 7;
    ista(5) = 6;
    iend(5) = 2;
    ista(6) = 3;
    iend(6) = 2;
    ista(7) = 2;
    iend(7) = 10;
    ista(8) = 10;
    iend(8) = 11;
    ista(9) = 10;
    iend(9) = 8;
    ista(10) = 8;
    iend(10) = 9;
    ista(11) = 0;
    iend(11) = 8;
    indexNodeConnectivity();

    // Organise cell dimensions
    findBoundingBox();
    alz = diam.max();

    // Compute lseg vector
    findLengths();

    uvec idx = find(nodtyp == 1);
    nnodbc = (int) idx.n_elem;
    bcnodname = nodname(idx);
    bcnod = zeros<ivec>(nnodbc);
    Bin = zeros<ivec>(nnodbc);
    Bout = zeros<ivec>(nnodbc);
    indexBCconnectivity();

    // Pairing of boundary nodes
    bcPairs = zeros<ivec>(nnodbc);
    bcPairs(0) = 1;
    bcPairs(1) = 1;
    bcPairs(2) = 2;
    bcPairs(3) = 3;
    bcPairs(4) = 2;
    bcPairs(5) = 3;

    // Defining in/out boundary nodes
    for (int inodbc = 0; inodbc < nnodbc; inodbc++)   {
        if (cnode(0,bcnod(inodbc)) == 0.)   {Bin(inodbc) = 1;}
        else if (cnode(0,bcnod(inodbc)) == alx)    {Bout(inodbc) = 1;}
        if (cnode(1,bcnod(inodbc)) == 0.)  {Bin(inodbc) = 1;}
        else if (cnode(1,bcnod(inodbc)) == aly) {Bout(inodbc) = 1;}
    }

    // Print micro-cell image
    pictureNetwork("MicroCellDiameters.ps", diam);

}


void MicroCell::crossCell2D() {

    cell2D = true;

    nnod = 5;
    nnodbc = 4;
    nseg = 4;

    setup_microCellArrays();
    for (int iseg = 0; iseg < nseg; iseg++)   {segname(iseg) = iseg;}
    for (int inod = 0; inod < nnod; inod++)   {nodname(inod) = inod;}

    // Randomly assign lengths from lognormal distribution
    double l_mean = mean(lengthDistrib);
    double l_SD = stddev(lengthDistrib);

    int n = 1;
    double L1 = sum(lognormal(l_mean, l_SD, n));
    double L2 = sum(lognormal(l_mean, l_SD, n));

    // Randomly assign diameters based on lognormal distribution
    double d_mean = mean(diamDistrib);
    double d_SD = stddev(diamDistrib);
    for (int iseg = 0; iseg < nseg; iseg++)   {diam(iseg) = sum(lognormal(d_mean, d_SD, n));}

    rseg = diam / 2.0;

    // Setup nodal coordinates
    cnode(1,0) = L2;

    cnode(0,1) = L1;
    cnode(1,1) = L2;

    cnode(0,2) = 2*L1;
    cnode(1,2) = L2;

    cnode(0,3) = L1;

    cnode(0,4) = L1;
    cnode(1,4) = 2*L2;
    if (abs(rotationAngle) > 0.)    {
        mat rotate = zeros<mat>(3,3);
        rotate(0,0) = cos(rotationAngle);
        rotate(0,1) = -sin(rotationAngle);
        rotate(1,0) = sin(rotationAngle);
        rotate(1,1) = cos(rotationAngle);
        rotate(2,2) = 1;
        for (int inod = 0; inod < nnod; inod++)   {cnode.col(inod) = rotate * cnode.col(inod);}
    }


    // Define node connections
    ista(0) = 0;
    iend(0) = 1;
    ista(1) = 1;
    iend(1) = 2;
    ista(2) = 1;
    iend(2) = 3;
    ista(3) = 4;
    iend(3) = 1;

    indexNodeConnectivity();

    // Organise cell dimensions
    findBoundingBox();
    alz = max(diam);

    // Compute lseg vector
    findLengths();

    // BC node names
    uvec idx = find(nodtyp == 1);
    nnodbc = (int) idx.n_elem;
    bcnodname = nodname(idx);
    bcnod = zeros<ivec>(nnodbc);
    Bin = zeros<ivec>(nnodbc);
    Bout = zeros<ivec>(nnodbc);
    indexBCconnectivity();

    // Pairing of boundary nodes
    bcPairs = zeros<ivec>(nnodbc);
    bcPairs(0) = 1;
    bcPairs(1) = 1;
    bcPairs(2) = 2;
    bcPairs(3) = 2;

    // Defining in/out boundary nodes
    for (int inodbc = 0; inodbc < nnodbc; inodbc++)   {
        // Planes perp. to x-axis
        if (cnode(0,bcnod(inodbc)) == 0.)   {Bin(inodbc) = 1;}
        else if (cnode(0,bcnod(inodbc)) == alx)    {Bout(inodbc) = 1;}
        // Planes perp. to y-axis
        if (cnode(1,bcnod(inodbc)) == 0.)  {Bin(inodbc) = 1;}
        else if (cnode(1,bcnod(inodbc)) == aly) {Bout(inodbc) = 1;}
    }

    // Print micro-cell image
    pictureNetwork("MicroCellDiameters.ps", diam);

}


void MicroCell::crossCell3D() {

    nnod = 7;
    nnodbc = 6;
    nseg = 6;

    setup_microCellArrays();
    for (int iseg = 0; iseg < nseg; iseg++)   {segname(iseg) = iseg;}
    for (int inod = 0; inod < nnod; inod++)   {nodname(inod) = inod;}

    // Randomly assign lengths from lognormal distribution
    double l_mean = mean(lengthDistrib);
    double l_SD = stddev(lengthDistrib);
    //mat lengths = branch_angles();
    int n = 1;
    double L1 = sum(lognormal(l_mean, l_SD, n));//sum(lognrnd(lengths(0,0), lengths(0,1), 1));
    double L2 = sum(lognormal(l_mean, l_SD, n));//sum(lognrnd(lengths(1,0), lengths(1,1), 1));
    double L3 = sum(lognormal(l_mean, l_SD, n));//sum(lognrnd(lengths(2,0), lengths(2,1), 1));

    // Randomly assign diameters based on lognormal distribution
    double d_mean = mean(diamDistrib);
    double d_SD = stddev(diamDistrib);
    for (int iseg = 0; iseg < nseg; iseg++)   {diam(iseg) = sum(lognormal(d_mean, d_SD, n));}
    rseg = diam / 2.0;

    // Setup nodal coordinates
    cnode(1,0) = L2;
    cnode(2,0) = L3;

    cnode(0,1) = L1;
    cnode(1,1) = L2;
    cnode(2,1) = L3;

    cnode(0,2) = 2*L1;
    cnode(1,2) = L2;
    cnode(2,2) = L3;

    cnode(0,3) = L1;
    cnode(2,3) = L3;

    cnode(0,4) = L1;
    cnode(1,4) = 2*L2;
    cnode(2,4) = L3;

    cnode(0,5) = L1;
    cnode(1,5) = L2;

    cnode(0,6) = L1;
    cnode(1,6) = L2;
    cnode(2,6) = 2*L3;

    if (abs(rotationAngle) > 0.) {
        mat rotate = zeros<mat>(3, 3);
        rotate(0, 0) = cos(rotationAngle);
        rotate(0, 1) = -sin(rotationAngle);
        rotate(1, 0) = sin(rotationAngle);
        rotate(1, 1) = cos(rotationAngle);
        rotate(2, 2) = 1;
        for (int inod = 0; inod < nnod; inod++) { cnode.col(inod) = rotate * cnode.col(inod); }
    }


    // Define node connections
    ista(0) = 0;
    iend(0) = 1;
    ista(1) = 1;
    iend(1) = 2;
    ista(2) = 1;
    iend(2) = 3;
    ista(3) = 4;
    iend(3) = 1;
    ista(4) = 1;
    iend(4) = 5;
    ista(5) = 6;
    iend(5) = 1;
    indexNodeConnectivity();

    // Organise cell dimensions
    findBoundingBox();

    // Compute lseg vector
    findLengths();

    // BC node names
    uvec idx = find(nodtyp == 1);
    nnodbc = (int) idx.n_elem;
    bcnodname = nodname(idx);
    bcnod = zeros<ivec>(nnodbc);
    Bin = zeros<ivec>(nnodbc);
    Bout = zeros<ivec>(nnodbc);
    indexBCconnectivity();

    // Pairing of boundary nodes
    bcPairs = zeros<ivec>(nnodbc);
    bcPairs(0) = 1;
    bcPairs(1) = 1;
    bcPairs(2) = 2;
    bcPairs(3) = 2;
    bcPairs(4) = 3;
    bcPairs(5) = 3;

    // Defining in/out boundary nodes
    for (int inodbc = 0; inodbc < nnodbc; inodbc++)   {
        // Planes perp. to x-axis
        if (cnode(0,bcnod(inodbc)) == 0.)   {Bin(inodbc) = 1;}
        else if (cnode(0,bcnod(inodbc)) == alx)    {Bout(inodbc) = 1;}
        // Planes perp. to y-axis
        if (cnode(1,bcnod(inodbc)) == 0.)  {Bin(inodbc) = 1;}
        else if (cnode(1,bcnod(inodbc)) == aly) {Bout(inodbc) = 1;}
        // Planes perp. to z-axis
        if (cnode(2,bcnod(inodbc)) == 0.)  {Bin(inodbc) = 1;}
        else if (cnode(2,bcnod(inodbc)) == alz) {Bout(inodbc) = 1;}
    }

    // Print micro-cell image
    pictureNetwork("MicroCellDiameters.ps", diam);

}
