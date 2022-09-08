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

    double l_mean = mean(eucLengths);
    double l_SD = stddev(eucLengths);
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
    l_mean = mean(lengthDistrib);
    l_SD = stddev(lengthDistrib);
    double length{},seg{};
    for (int iseg = 0; iseg < nseg; iseg++) {
        seg = lseg(iseg);
        do {length = sum(lognormal(l_mean, l_SD, n));}
        while (length < seg);
        lseg(iseg) = length;
    }

    uvec idx = find(nodtyp == 1);
    nnodbc = (int) idx.n_elem;
    bcnodname = nodname(idx);
    bcnod = zeros<uvec>(nnodbc);
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
    double l_mean = mean(eucLengths);
    double l_SD = stddev(eucLengths);

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
    l_mean = mean(lengthDistrib);
    l_SD = stddev(lengthDistrib);
    double length{},seg{};
    for (int iseg = 0; iseg < nseg; iseg++) {
        seg = lseg(iseg);
        while (length < seg) {length = sum(lognormal(l_mean, l_SD, n));}
        lseg(iseg) = length;
    }

    // BC node names
    uvec idx = find(nodtyp == 1);
    nnodbc = (int) idx.n_elem;
    bcnodname = nodname(idx);
    bcnod = zeros<uvec>(nnodbc);
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


/*void MicroCell::crossCell3D() {

    nnod = 7;
    nnodbc = 6;
    nseg = 6;

    setup_microCellArrays();
    for (int iseg = 0; iseg < nseg; iseg++)   {segname(iseg) = iseg;}
    for (int inod = 0; inod < nnod; inod++)   {nodname(inod) = inod;}

    // Randomly assign lengths from lognormal distribution
    double l_mean = mean(eucLengths);
    double l_SD = stddev(eucLengths);

    int n = 1;
    double L1 = sum(lognormal(l_mean, l_SD, n));
    double L2 = sum(lognormal(l_mean, l_SD, n));
    double L3 = sum(lognormal(l_mean, l_SD, n));

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
        for (int inod = 0; inod < nnod; inod++) {cnode.col(inod) = rotate * cnode.col(inod); }
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
    l_mean = mean(lengthDistrib);
    l_SD = stddev(lengthDistrib);
    double length{},seg{};
    for (int iseg = 0; iseg < nseg; iseg++) {
        seg = lseg(iseg);
        while (length < seg) {length = sum(lognormal(l_mean, l_SD, n));}
        lseg(iseg) = length;
    }

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

}*/


void MicroCell::crossCell3D() {

    // No. of tessellations in each coordinate
    int xGrid = 10;
    int yGrid = 10;
    int zGrid = 10;
    gridSpacing = zeros<vec>(3);
    gridSpacing(0) = xGrid;
    gridSpacing(1) = yGrid;
    gridSpacing(2) = zGrid;

    // Initial setup
    nnod = 7;
    nseg = 6;
    nodsegm = 6;
    cnode = zeros<mat>(3,nnod);
    ista = zeros<ivec>(nseg);
    iend = zeros<ivec>(nseg);

    // Initial length
    double initLength = 0.5;

    // Setup nodal coordinates of initial grid
    cnode(1,0) = initLength;
    cnode(2,0) = initLength;

    cnode(0,1) = initLength;
    cnode(1,1) = initLength;
    cnode(2,1) = initLength;

    cnode(0,2) = 2*initLength;
    cnode(1,2) = initLength;
    cnode(2,2) = initLength;

    cnode(0,3) = initLength;
    cnode(2,3) = initLength;

    cnode(0,4) = initLength;
    cnode(1,4) = 2*initLength;
    cnode(2,4) = initLength;

    cnode(0,5) = initLength;
    cnode(1,5) = initLength;

    cnode(0,6) = initLength;
    cnode(1,6) = initLength;
    cnode(2,6) = 2*initLength;

    // Define node connections
    ista(0) = 0;
    iend(0) = 1;
    ista(1) = 1;
    iend(1) = 2;
    ista(2) = 3;
    iend(2) = 1;
    ista(3) = 1;
    iend(3) = 4;
    ista(4) = 5;
    iend(4) = 1;
    ista(5) = 1;
    iend(5) = 6;

    // Tessellate 3D grid
    ivec istaCopy = ista;
    ivec iendCopy = iend;
    mat cnodeCopy = cnode;
    int maxidx{};
    ivec tmpVec;
    mat tmpMat;
    bool skip = true;
    for (int i = 0; i < xGrid; i++)    {
        for (int j = 0; j < yGrid; j++)    {
            for (int k = 0; k < zGrid; k++)    {
                if (!skip)    {
                    tmpMat = cnodeCopy;
                    tmpMat.row(0) += i;
                    tmpMat.row(1) += j;
                    tmpMat.row(2) += k;
                    cnode.insert_cols(cnode.n_cols, tmpMat);

                    maxidx = max(ista.max(), iend.max());

                    tmpVec = istaCopy;
                    tmpVec += (maxidx + 1);
                    ista.insert_rows(ista.n_rows, tmpVec);

                    tmpVec = iendCopy;
                    tmpVec += (maxidx + 1);
                    iend.insert_rows(iend.n_rows, tmpVec);
                }
                else {skip = false;}
            }
        }
    }

    // Find spatially overlapping nodes and connect tessellations
    nseg = (int) ista.n_elem;
    nnod = (int) cnode.n_cols;
    indexNodeConnectivity();
    vec nod{};
    ivec nodflag = zeros<ivec>(nnod);
    for (int inod = 0; inod < nnod; inod++) {
        nod = cnode.col(inod);
        for (int jnod = inod + 1; jnod < nnod; jnod++) {
            if (nod(0) == cnode(0, jnod) && nod(1) == cnode(1, jnod) && nod(2) == cnode(2, jnod))   {
                nodflag(jnod) = 1;
                for (int iseg = 0; iseg < nseg; iseg++) {
                    if (ista(iseg) == jnod) {
                        ista(iseg) = inod;
                        iseg = nseg;
                    }
                    else if (iend(iseg) == jnod)    {
                        iend(iseg) = inod;
                        iseg = nseg;
                    }
                }
            }
        }
    }

    // Generate seg & nodal names. Remove redundant nodes
    nodname = zeros<ivec>(nnod);
    segname = zeros<ivec>(nseg);
    for (int inod = 0; inod < nnod; inod++)   {nodname(inod) = inod;}
    for(int iseg = 0; iseg < nseg; iseg++)  {segname(iseg) = iseg;}
    segnodname = zeros<imat>(2,nseg);
    for (int iseg = 0; iseg < nseg; iseg++) {
        segnodname(0,iseg) = nodname(ista(iseg));
        segnodname(1,iseg)= nodname(iend(iseg));
    }
    nodname.shed_rows(find(nodflag == 1));
    cnode.shed_cols(find(nodflag == 1));
    nnod = (int) cnode.n_cols;
    indexSegmentConnectivity();
    indexNodeConnectivity();

    // Remove nodes of nodtyp == 2
    nodflag = zeros<ivec>(nnod);
    int seg1{},seg2{},nod2{};
    ivec segflag = zeros<ivec>(nseg);
    for (int inod = 0; inod < nnod; inod++) {
        if (nodtyp(inod) == 2)  {
            seg1 = nodseg(0, inod);
            seg2 = nodseg(1, inod);
            nod2 = nodnod(1,inod);
            if (segnodname(0, seg1) == nodname(inod))   {segnodname(0, seg1) = nodname(nod2);}
            else if (segnodname(1, seg1) == nodname(inod))   {segnodname(1, seg1) = nodname(nod2);}
            nodflag(inod) = 1;
            segflag(seg2) = 1;
        }
    }
    segname.shed_rows(find(segflag == 1));
    segnodname.shed_cols(find(segflag == 1));
    nodname.shed_rows(find(nodflag == 1));
    nodtyp.shed_rows(find(nodflag == 1));
    cnode.shed_cols(find(nodflag == 1));
    nseg = (int) segname.n_elem;
    nnod = (int) cnode.n_cols;
    indexSegmentConnectivity();
    indexNodeConnectivity();

    // Setup remaining arrays
    uvec idx = find(nodtyp == 1);
    nnodbc = (int) idx.n_elem;
    bcnodname = nodname(idx);
    setup_microCellArrays(true);
    indexBCconnectivity();

    findBoundingBox();
    findBCpairs();
    int bnod{},knodbc{};
    int minRemoval = round(0.4*nnodbc);
    int maxRemoval = round(0.45*nnodbc);
    ivec flagSeg = zeros<ivec>(nseg);
    ivec removalShuffle = linspace<ivec>(minRemoval, maxRemoval, (maxRemoval-minRemoval)+1);
    ivec bcidx = linspace<ivec>(0,nnodbc-1,nnodbc);
    removalShuffle = shuffle(removalShuffle);
    bcidx = shuffle(bcidx);
/*    for (int inodbc = 0; inodbc < removalShuffle(0); inodbc++) {
        knodbc = bcidx(inodbc);
        for (int iseg = 0; iseg < nseg; iseg++) {
            if (ista(iseg) == bcnod(knodbc) || iend(iseg) == bcnod(knodbc)) {
                flagSeg(iseg) = 1;
                iseg = nseg;
            }
        }
        uvec oppBC = find(bcPairs(knodbc) == bcPairs);
        for (int jnodbc = 0; jnodbc < 2; jnodbc++)  {
            bnod = oppBC(jnodbc);
            for (int iseg = 0; iseg < nseg; iseg++) {
                if (ista(iseg) == bcnod(bnod) || iend(iseg) == bcnod(bnod)) {
                    flagSeg(iseg) = 1;
                    iseg = nseg;
                }
            }
        }
    }
    bctyp = zeros<ivec>(nnodbc);
    bcprfl = zeros<vec>(nnodbc);
    bchd = zeros<vec>(nnodbc);
    subNetwork(flagSeg, true, false);*/

    // Randomly remove segments to reduce nodal connectivity
    int nodsum{},segsum{};
    uvec jdx;
    flagSeg = zeros<ivec>(nseg);
    for (int inodbc = 0; inodbc < nnodbc; inodbc++) {
        flagSeg(nodseg(0,bcnod(inodbc))) = -1;
    }

    ivec copytyp = nodtyp;
    int nod3{},seg{};
    ivec segshuffle = linspace<ivec>(0,nseg-1,nseg);
    segshuffle = shuffle(segshuffle);
    for (int iseg = 0; iseg < nseg; iseg++) {
        seg = segshuffle(iseg);
        if (flagSeg(seg) == 0) {
            nod3 = ista(seg);
            nod2 = iend(seg);
            if (copytyp(nod3) == 1 && copytyp(nod2) > 3) {
                flagSeg(seg) = 1;
                copytyp(nod2) -= 1;
                copytyp(nod3) -= 1;
            }
            else if (copytyp(nod3) > 3 && copytyp(nod2) == 1) {
                flagSeg(seg) = 1;
                copytyp(nod2) -= 1;
                copytyp(nod3) -= 1;
            }
        }
    }
    segshuffle = shuffle(segshuffle);
    for (int iseg = 0; iseg < nseg; iseg++) {
        seg = segshuffle(iseg);
        if (flagSeg(seg) == 0) {
            nod3 = ista(seg);
            nod2 = iend(seg);
            if (copytyp(nod3) > 3 && copytyp(nod2) > 3) {
                flagSeg(seg) = 1;
                copytyp(nod2) -= 1;
                copytyp(nod3) -= 1;
            }
        }
    }
    segshuffle = shuffle(segshuffle);
    for (int iseg = 0; iseg < nseg; iseg++) {
        seg = segshuffle(iseg);
        if (flagSeg(seg) == 0) {
            nod3 = ista(seg);
            nod2 = iend(seg);
            if (copytyp(nod3) > 3 && copytyp(nod2) > 2)    {
                flagSeg(seg) = 1;
                copytyp(nod2) -= 1;
                copytyp(nod3) -= 1;
            }
            else if (copytyp(nod3) > 2 && copytyp(nod2) > 3)    {
                flagSeg(seg) = 1;
                copytyp(nod2) -= 1;
                copytyp(nod3) -= 1;
            }
        }
    }

    flagSeg(find(flagSeg == -1)).fill(0);
    bctyp = zeros<ivec>(nnodbc);
    bcprfl = zeros<vec>(nnodbc);
    bchd = zeros<vec>(nnodbc);
    subNetwork(flagSeg, true, false);

    // Setup microcell arrays
    setup_microCellArrays(true);
    indexBCconnectivity();
    findBoundingBox();
    findBCpairs();

    // Organise cell dimensions
    findBoundingBox();

    // Classify vessels aligned with axes
    int nod1{};
    for (int iseg = 0; iseg < nseg; iseg++) {
        nod1 = ista(iseg);
        nod2 = iend(iseg);
        if (cnode(1,nod1) == cnode(1,nod2) && cnode(2,nod1) == cnode(2,nod2)) {
            vessOrient(iseg) = 1; // Parallel to x-axis
        }
        else if (cnode(0,nod1) == cnode(0,nod2) && cnode(2,nod1) == cnode(2,nod2))   {
            vessOrient(iseg) = 2; // Parallel to y-axis
        }
        else if (cnode(0,nod1) == cnode(0,nod2) && cnode(1,nod1) == cnode(1,nod2))   {
            vessOrient(iseg) = 3; // Parallel to z-axis
        }
        else {
            vessOrient(iseg) = 4; // Diagonal vessels
        }
    }

}

void MicroCell::assignGeometry()   {

    // To remove any artificial enlargement
    findBoundingBox();

    // Normalise cell size
    cnode.row(0) /= alx;
    cnode.row(1) /= aly;
    cnode.row(2) /= alz;
    alx = 1.;
    aly = 1.;
    alz = 1.;

    // Randomly assign diameters based on lognormal distribution
    int n{1};
    double sample{};
    double d_mean = mean(diamDistrib);
    double d_SD = stddev(diamDistrib);
    double dmin = min(diamDistrib);
    for (int iseg = 0; iseg < nseg; iseg++)   {
        sample = 0.;
        while (sample < dmin)    {sample = sum(lognormal(d_mean, d_SD, n));}
        diam(iseg) = sample;
    }
    //diam.fill(5.);
    rseg = diam / 2.0;

    // Calculate euclidean lengths
    double l_mean{},l_SD{};
    l_mean = mean(eucLengths);
    l_SD = stddev(eucLengths);

    vec gridLength = zeros<vec>(3);
    for (int i = 0; i < 3; i++) {
        sample = sum(lognormal(l_mean, l_SD, n));
        gridLength(i) = sample;
        cnode.row(i) *= (sample  * gridSpacing(i));
    }

    l_mean = 5.104;//mean(lengthDistrib);
    l_SD = 4.579;//stddev(lengthDistrib);
    int nitmax{100},cntr{};
    for (int iseg = 0; iseg < nseg; iseg++)   {
        sample = 0.;
        cntr = 0;
        while (sample * diam(iseg) < min(gridLength) && cntr < nitmax){
            sample = sum(lognormal(l_mean, l_SD, n));
            cntr += 1;
        }
        lseg(iseg) = sample * diam(iseg);
    }

    // Ensure diameter periodicity at boundaries
    double bcDiam{};
/*    for (int inodbc = 0; inodbc < nnodbc; inodbc++) {
        for (int iseg = 0; iseg < nseg; iseg++) {
            if (ista(iseg) == bcnod(inodbc) || iend(iseg) == bcnod(inodbc)) {
                bcDiam = diam(iseg);
                lseg(iseg) /= 2.;
                iseg = nseg;
            }
        }
        for (int jnodbc = 0; jnodbc < nnodbc; jnodbc++) {
            if (bcPairs(inodbc) == bcPairs(jnodbc)) {
                for (int iseg = 0; iseg < nseg; iseg++) {
                    if (ista(iseg) == bcnod(jnodbc) || iend(iseg) == bcnod(jnodbc)) {
                        diam(iseg) = bcDiam;
                        rseg(iseg) = 0.5 * bcDiam;
                        lseg(iseg) /= 2.;
                        iseg = nseg;
                        jnodbc = nnodbc;
                    }
                }
            }
        }
    }*/

    // Artificially inflate tissue box due to "half-sized" boundary segments
    findBoundingBox();
    alx += gridLength(0);
    aly += gridLength(1);
    alz += gridLength(2);

    // Network volume
    netVol = sum(lseg % pow(rseg, 2))*M_PI;

    // Vascular density
    vascDens = 100*netVol/(alx*aly*alz);

    // Print micro-cell image
/*    int cntr = 0;
    vec flagPairs = zeros<vec>(nseg);
    for (int inodbc = 0; inodbc < nnodbc; inodbc++) {
        cntr += 1;
        for (int iseg = 0; iseg < nseg; iseg++) {
            if (ista(iseg) == bcnod(inodbc) || iend(iseg) == bcnod(inodbc)) {
                flagPairs(iseg) = cntr;
                iseg = nseg;
            }
        }
        for (int jnodbc = 0; jnodbc < nnodbc; jnodbc++) {
            if (bcPairs(inodbc) == bcPairs(jnodbc)) {
                for (int iseg = 0; iseg < nseg; iseg++) {
                    if (ista(iseg) == bcnod(jnodbc) || iend(iseg) == bcnod(jnodbc)) {
                        flagPairs(iseg) = cntr;
                        iseg = nseg;
                        jnodbc = nnodbc;
                    }
                }
            }
        }
    }*/
    pictureNetwork("MicroCellDiameters.ps", diam);

}

void MicroCell::findBCpairs()   {

    // Pairing of boundary nodes
    int cntr{},nod1{},nod2{},pairIdx{1};
    vec gridDim = zeros<vec>(3);
    gridDim(0) = alx;
    gridDim(1) = aly;
    gridDim(2) = alz;
    for (int i = 0; i < 3; i++) {
        for (int inodbc = 0; inodbc < nnodbc; inodbc++) {
            cntr = 0;
            nod1 = bcnod(inodbc);
            if (cnode(i,nod1) == 0.) {
                for (int jnodbc = 0; jnodbc < nnodbc; jnodbc++) {
                    nod2 = bcnod(jnodbc);
                    if (inodbc != jnodbc && cnode(i,nod2) == gridDim(i))   {
                        for (int j = 0; j < 3; j++) {
                            if (i != j && cnode(j,nod1) == cnode(j,nod2)) {
                                cntr += 1;
                            }
                        }
                        if (cntr == 2)  {
                            bcPairs(inodbc) = pairIdx;
                            Bin(inodbc) = 1;
                            bcPairs(jnodbc) = pairIdx;
                            Bout(jnodbc) = 1;
                            pairIdx += 1;
                            jnodbc = nnodbc;
                        }
                        else {cntr = 0;}
                    }
                }
            }
        }
    }

}
