#include "Network.hpp"

using namespace reanimate;

void Network::setup_networkArrays() {

    ista = zeros<ivec>(nseg);
    iend = zeros<ivec>(nseg);
    nodout = zeros<ivec>(nnod);
    nodrank = zeros<ivec>(nnod);
    nk = zeros<ivec>(nnod);
    sgraphTag = zeros<ivec>(nseg);
    ngraphTag = zeros<ivec>(nnod);

    nodtyp = zeros<ivec>(nnod);
    bcnod = zeros<ivec>(nnodbc);
    articPnt = zeros<ivec>(nnod);

    nodnod = zeros<imat>(nodsegm,nnod);
    nodseg = zeros<imat>(nodsegm,nnod);

}

void Network::setup_flowArrays(bool popMatrices)    {

    c = zeros<vec>(nseg);
    conductance = zeros<vec>(nseg);
    segpress = zeros<vec>(nseg);
    tau = zeros<vec>(nseg);
    deadends = zeros<ivec>(nseg);
    subGraphs = zeros<ivec>(nseg);
    BCflow = zeros<vec>(nnodbc);
    BCpress = zeros<vec>(nnodbc);
    Qo = zeros<vec>(nnod);

    noflow = zeros<ivec>(nseg);

    M = zeros<sp_mat>(nseg,nnod);
    L = zeros<sp_mat>(nnod,nseg);

    qold = q;
    hdold = hd;

    if (!unknownBCs && popMatrices)    {

        for (int inod = 0; inod < nnod; inod++)    {
            int countIndx = 0;
            for (int iseg = 0; iseg < nseg; iseg++)    {
                if (inod == ista(iseg))   {
                    L(inod,iseg) = -1.;
                    countIndx += 1;
                }
                else if (inod == iend(iseg))  {
                    L(inod,iseg) = 1.;
                    countIndx += 1;
                }
                if (nodtyp(inod) == countIndx)  {
                    iseg = nseg;
                }
            }
        }
    }


}

void Network::setup_estimationArrays()  {

    // Creating an index to flag boundary nodes with unknown boundary conditions
    unknownnod = zeros<ivec>(nnod);
    for (int inodbc = 0; inodbc < nnodbc; inodbc++) {
        if (bctyp(inodbc) == 3)    {unknownnod(bcnod(inodbc)) = 1;}
    }
    unknownnod_idx = find(unknownnod == 0);
    bcpress_idx = find(bctyp == 0);
    nunknown = (int) accu(unknownnod);
    printText("Total unknown conditions = "+to_string(nunknown)+" ("+to_string(100*(float(nunknown)/float(nnodbc)))+"%)");
    nIBnod = nnod - nunknown; // No. of internal and known boundary nodes
    estimationarraysize = nnod + nIBnod;

    // Flow reversal storage
    storeBC = bcprfl;
    storeBCtyp = bctyp;
    storeBChd = bchd;
    storeHD = hd;

    // Store previous flow iterations
    oldFlowsign = zeros<vec>(nseg);
    oldTau = zeros<vec>(nseg);
    oldHd = zeros<vec>(nseg);
    oldq = zeros<vec>(nseg);
    oldNodpress = zeros<vec>(nnod);

    flowsign = zeros<vec>(nseg);
    tau0 = zeros<vec>(nseg);
    p0 = zeros<vec>(nnod);
    Qo = zeros<vec>(nIBnod);
    B = zeros<vec>(estimationarraysize);

    A = zeros<sp_mat>(estimationarraysize,estimationarraysize);
    H1 = zeros<sp_mat>(nnod,nseg);
    H2 = zeros<sp_mat>(nseg,nnod);
    W = zeros<sp_mat>(nnod,nnod);
    L = zeros<sp_mat>(nIBnod,nseg);

    // Check haematocrit
    if (accu(hd) == 0.) {hd.fill(consthd);}
    if (accu(bchd) == 0.)   {hd.fill(consthd);}

    // Assign target pressure
    //targPress = 31.;
    targPress = mean(bcprfl(find(bctyp == 0))); // Set target pressure as the mean of the assign boundary pressure conditions
    printNum("Target Pressure = ",targPress,"mmHg");
    p0.fill(targPress * alpha);

    // Target shear stress - intially set with random directions unless network flow is known
    targStress = 15.;
    printNum("Target Wall Shear Stress = ",targStress,"dyn/cm2");
    targStress *= beta;
    int ran{};
    srand((u_int) time(0));
    for (int iseg = 0; iseg < nseg; iseg++)    {
        ran = rand()%2;
        if (ran == 0)   {
            ran = -1;
        }
        tau0(iseg) = targStress*ran;
    }
    if (accu(q) > 0.0)  {
        tau0 = abs(tau0);
        tau0 %= sign(q);
        printText("Wall shear stress initialised using network file flow direction",2,0);
    }


    // Constructing vector Qo
    int cntr{};
    for (int inod = 0; inod < nnod; inod++) {
        if (unknownnod(inod) == 0)  {
            for (int inodbc = 0; inodbc < nnodbc; inodbc++) {
                if (inod == bcnod(inodbc) && bctyp(inodbc) == 0) {
                    Qo(inod - cntr) = -bcprfl(inodbc)*alpha;
                }
                else if (inod == bcnod(inodbc) && (bctyp(inodbc) == 1 || bctyp(inodbc) == 2))  {
                    Qo(inod - cntr) = bcprfl(inodbc)*gamma;
                }
            }
        }
        else {
            cntr += 1;
        }
    }


    // Construct W matrix
    for (int iseg = 0; iseg < nseg; iseg++) {
        long tista = ista(iseg);
        long tiend = iend(iseg);
        double tlseg = lseg(iseg);
        W(tista,tista) += (0.5*kp*tlseg);
        W(tiend,tiend) += (0.5*kp*tlseg);
    }


    // Construct L matrix
    cntr = 0;
    for (int inod = 0; inod < nnod; inod++)    {
        if (unknownnod(inod) == 0)  {
            int idx = 0;
            for (int iseg = 0; iseg < nseg; iseg++)    {
                if (inod == ista(iseg))   {
                    L(inod-cntr,iseg) = -1.;
                    idx += 1;
                }
                else if (inod == iend(iseg))  {
                    L(inod-cntr,iseg) = 1.;
                    idx += 1;
                }
                if (nodtyp(inod) == idx)  {
                    iseg = nseg;
                }
            }
        }
        else    {
            cntr += 1;
        }
    }

    // Partially construct B vector
    B(span(0,nIBnod-1)) = -Qo;

}