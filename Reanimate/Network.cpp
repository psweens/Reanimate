//
//  Network.cpp
//  TEEPEE
//
//  Created by Paul Sweeney on 16/06/2020.
//  Copyright © 2020 Paul Sweeney. All rights reserved.
//

#include "Network.hpp"
#include "misc_func.h"

#include <string>

using namespace std;
using namespace reanimate;

Network::Network() {

    gamma = 1./(1.e3*60);
    alpha = 0.1333;
    beta = 1.e-4;
    xi = 0.001*1e-3;

    consthd = 0.45;

    kp = 0.1;
    ktau = 1.e-4;

    unknownBCs=false;

}
Network::~Network() = default;

void Network::setup_networkArrays() {

    ista = zeros<ivec>(nseg);
    iend = zeros<ivec>(nseg);
    nodout = zeros<ivec>(nnod);
    nodrank = zeros<ivec>(nnod);
    nk = zeros<ivec>(nnod);

    rseg = zeros<vec>(nseg);
    nodtyp = zeros<ivec>(nnod);
    bcnod = zeros<ivec>(nnodbc);

    nodnod = zeros<imat>(nodsegm,nnod);
    nodseg = zeros<imat>(nodsegm,nnod);

}

void Network::setup_flowArrays()    {

    c = zeros<vec>(nseg);
    conductance = zeros<vec>(nseg);
    segpress = zeros<vec>(nseg);
    BCflow = zeros<vec>(nnodbc);
    BCpress = zeros<vec>(nnodbc);
    Qo = zeros<vec>(nnod);

    zeroflow = zeros<ivec>(nseg);

    M = zeros<sp_mat>(nseg,nnod);
    L = zeros<sp_mat>(nnod,nseg);


    if (!unknownBCs)    {

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
        if (bctyp(inodbc) == 3)    {
            unknownnod(bcnod(inodbc)) = 1;
        }
    }
    nunknown = (int) accu(unknownnod);
    cout<<"Total unknown conditions = "<<nunknown<<" ("<<100*(float(nunknown)/float(nnodbc))<<"%)"<<endl;
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
    targPress = 31.;//mean(bcprfl(find(bctyp == 0))); // Set target pressure as the mean of the assign boundary pressure conditions
    cout<<"Target Pressure = "<<targPress<<" mmHg"<<endl;
    p0.fill(targPress * alpha);

    // Target shear stress - intially set with random directions unless network flow is known
    targStress = 15.;
    cout<<"Target Wall Shear Stress = "<<targStress<<" dyn/cm2"<<endl;
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
        cout<<"Wall shear stress initialised using network file flow direction ..."<<endl;
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

void Network::loadNetwork(const string& filepath)   {
    
    int max=200;
    char bb[200];

    networkPath = filepath;

    FILE *ifp;

    cout<<"\nImporting network data..."<<endl;
    ifp = fopen(filepath.c_str(),"r");

    fgets(bb,max,ifp);
    printf("%s",bb);
    
    networkName = bb;
    
    fscanf(ifp, "%f %f %f", &alx,&aly,&alz); fgets(bb,max,ifp);
    fscanf(ifp, "%i %i %i", &mxx,&myy,&mzz); fgets(bb,max,ifp);
    fscanf(ifp, "%f", &lb); fgets(bb,max,ifp);
    fscanf(ifp, "%f", &maxl); fgets(bb,max,ifp);
    fscanf(ifp,"%i", &nodsegm);
    fgets(bb,max,ifp);
    
    // Number of segments in vessel network
    fscanf(ifp,"%i", &nseg); fgets(bb,max,ifp);
    fgets(bb,max,ifp);

    // Segment properties: name type nodename(start), nodename(end), diameter, flow, hematocrit
    segname = zeros<ivec>(nseg);
    vesstyp = zeros<ivec>(nseg);
    segnodname = zeros<imat>(2,nseg);
    diam = zeros<vec>(nseg);
    lseg = zeros<vec>(nseg);
    q = zeros<vec>(nseg);
    hd = zeros<vec>(nseg);

    int num = detect_col(ifp);
    if (num == 7)   {
        for(int iseg = 0; iseg < nseg; iseg++){
            fscanf(ifp, "%lli %lli %lli %lli %lf %lf %lf\n",
                   &segname(iseg),&vesstyp(iseg),&segnodname(0,iseg),&segnodname(1,iseg),&diam(iseg),&q(iseg),&hd(iseg));
        }
        computeLseg = 1;
    }
    else if (num == 8)  {
        for(int iseg = 0; iseg < nseg; iseg++){
            fscanf(ifp, "%lli %lli %lli %lli %lf %lf %lf %lf\n",
                   &segname(iseg),&vesstyp(iseg),&segnodname(0,iseg),&segnodname(1,iseg),&diam(iseg),&lseg(iseg),&q(iseg),&hd(iseg));
        }
    }
    else    {
        printf("*** Error in Network File: Invalid Segment Format ***");
    }

    
    // Number of nodes in vessel network
    fscanf(ifp,"%i", &nnod);
    fgets(bb,max,ifp);
    fgets(bb,max,ifp);
    
    // Coordinates of nodes
    nodname = zeros<ivec>(nnod);
    cnode = zeros<mat>(3,nnod);
    nodpress = zeros<vec>(nnod);
    
    num = detect_col(ifp);
    if (num == 4)   {
        for(int inod = 0; inod < nnod; inod++)  {
            fscanf(ifp, "%lli %lf %lf %lf\n", &nodname(inod),&cnode(0,inod),&cnode(1,inod),&cnode(2,inod));
        }
    }
    else if (num == 5)  {
        for(int inod = 0; inod < nnod; inod++)  {
            fscanf(ifp, "%lli %lf %lf %lf %lf\n", &nodname(inod),&cnode(0,inod),&cnode(1,inod),&cnode(2,inod),&nodpress(inod));
        }
    }
    else    {
        printf("*** Error in Network File: Invalid Node Format ***");
    }
    
    // Boundary nodes
    fscanf(ifp,"%i", &nnodbc);
    fgets(bb,max,ifp);
    fgets(bb,max,ifp);
    
    bcnodname = zeros<ivec>(nnodbc);
    bctyp = zeros<ivec>(nnodbc);
    bcprfl = zeros<vec>(nnodbc);
    bchd = zeros<vec>(nnodbc);
    
    nsol = detect_col(ifp);
    if (nsol == 4)   {
        for(int inodbc = 0; inodbc < nnodbc; inodbc++){
            fscanf(ifp,"%lli %lli %lf %lf\n", &bcnodname(inodbc),&bctyp(inodbc),&bcprfl(inodbc),&bchd(inodbc));
        }
        nsol -= 3;
        bcp = zeros<mat>(nnodbc,nsol);
    }
    else if (nsol > 4)  {
        for(int inodbc = 0; inodbc < nnodbc; inodbc++){
            nsol -= 4; // Count extra solutes (exc. PO2)
            bcp = zeros<mat>(nnodbc,nsol);
            fscanf(ifp,"%lli %lli %lf %lf", &bcnodname(inodbc),&bctyp(inodbc),&bcprfl(inodbc),&bchd(inodbc));
            for (int i = 0; i < nsol; i++)  {
                fscanf(ifp," %lf",&bcp(inodbc,i));
            }
            fscanf(ifp,"\n");
        }
    }
    else    {
        printf("*** Error in Network File: Invalid Boundary Node Format ***");
    }
    
    fclose(ifp);
    
    nodsegm += 1; // Armadillo indexing starts at zero

    setup_networkArrays();
    analyse_network();
}

void Network::analyse_network()   {

    printf("No. of segments = %i\n",nseg);
    printf("No. of nodes = %i\n",nnod);
    printf("No. of boundary nodes = %i\n",nnodbc);

    for (int iseg = 0; iseg < nseg; iseg++)	{
        //Search for nodes corresponding to this segment
        for (int i = 0; i < 2; i++) {
            for (int inod = 0; inod < nnod; inod++) {
                if (nodname(inod) == segnodname(i,iseg))    {
                    if (i == 0) {
                        ista(iseg) = inod;
                        goto foundit;
                    }
                    else if (i == 1)    {
                        iend(iseg) = inod;
                        goto foundit;
                    }
                }
            }
            printf("*** Error: No matching node found for segname %lli\n", segname(iseg));
            foundit:;
        }
    }


    // Setup nodtyp, nodseg and nodnod
    // 'nodseg' -> for each nodes, store the corresponding segment index
    // 'nodnod' -> for each node, store the nodal index of the opposite side of the segment
    for (int iseg = 0; iseg < nseg; iseg++) {
        int inod1 = (int) ista(iseg);
        int inod2 = (int) iend(iseg);
        nodtyp(inod1) += 1;
        nodtyp(inod2) += 1;
        if (nodtyp(inod1) > nodsegm) {
            printf("*** Error: Too many segments connected to node %i\n", inod1);
        }
        if (nodtyp(inod2) > nodsegm) {
            printf("*** Error: Too many segments connected to node %i\n", inod2);
        }
        // "-1" due to indexing starting at zero
        nodseg(nodtyp(inod1) - 1,inod1) = iseg;
        nodseg(nodtyp(inod2) - 1,inod2) = iseg;
        nodnod(nodtyp(inod1) - 1,inod1) = inod2;
        nodnod(nodtyp(inod2) - 1,inod2) = inod1;
    }

    for (int inodbc = 0; inodbc < nnodbc; inodbc++){
        // Search for node corresponding to this node name
        for (int inod = 0; inod < nnod; inod++) {
            if(nodname(inod) == bcnodname(inodbc))  {
                bcnod(inodbc) = inod;
                if(nodtyp(inod) != 1)   {
                    printf("*** Error: Boundary node %lli is not a 1-segment node\n", nodname(inod));
                }
                goto foundit2;
            }
        }
        printf("*** Error: No matching node found for nodname %lli, inodbc %i\n", bcnodname(inodbc),inodbc);
        foundit2:;
    }


    // Start(k,iseg) = coordinates of starting point of segment iseg
    // End(k,iseg) = coordinates of ending point of segment iseg
    computeLseg = 1;
    qq = abs(q);
    vec ss = zeros<vec>(3);
    mat Start = zeros<mat>(3,nseg);
    mat End = zeros<mat>(3,nseg);
    for (int iseg = 0; iseg < nseg; iseg++) {
        rseg(iseg) = diam(iseg) / 2.0;
        qq(iseg) = abs(q(iseg));
        for (int k = 0; k < 3; k++){
            Start(k,iseg) = cnode(k,ista(iseg));
            End(k,iseg) = cnode(k,iend(iseg));
            ss(k) = End(k,iseg) - Start(k,iseg);
        }
        if (computeLseg == 1)  {
            lseg(iseg) = sqrt(pow(ss(0),2) + pow(ss(1),2) + pow(ss(2),2));
        }
    }

    BCgeo = zeros<ivec>(nnodbc);
    for (int inodbc = 0; inodbc < nnodbc; inodbc++) {
        for (int iseg = 0; iseg < nseg; iseg++) {
            if (bcnod(inodbc) == ista(iseg) || bcnod(inodbc) == iend(iseg)) {
                if (vesstyp(iseg) == 1)    {
                    BCgeo(inodbc) = 1;
                }
                else if (vesstyp(iseg) == 2)   {
                    BCgeo(inodbc) = 2;
                }
                else {
                    BCgeo(inodbc) = 3;
                }
            }
        }
    }

}

void Network::subNetwork(ivec &index) {
    
    cout<<"Creating subnetwork..."<<endl;
    for (int iseg = 0; iseg < nseg; iseg++) {
        if (index(iseg) == 1)   {
            segname.shed_row(iseg);
            vesstyp.shed_row(iseg);
            segnodname.shed_col(iseg);
            diam.shed_row(iseg);
            lseg.shed_row(iseg);
            q.shed_row(iseg);
            hd.shed_row(iseg);
            index.shed_row(iseg);
            nseg -= 1;
            iseg = 0;
        }
    }
    
    ista = zeros<ivec>(nseg);
    iend = zeros<ivec>(nseg);
    for (int iseg = 0; iseg < nseg; iseg++)    {
        //Search for nodes corresponding to this segment
        for (int i = 0; i < 2; i++) {
            for (int inod = 0; inod < nnod; inod++) {
                if (nodname(inod) == segnodname(i,iseg))    {
                    if (i == 0) {
                        ista(iseg) = inod;
                        goto foundit;
                    }
                    else if (i == 1)    {
                        iend(iseg) = inod;
                        goto foundit;
                    }
                }
            }
            printf("*** Error: No matching node found for segname %lli\n", segname(iseg));
        foundit:;
        }
    }
    
    
    //Setup nodtyp, nodseg and nodnod
    nodtyp.zeros();
    for (int iseg = 0; iseg < nseg; iseg++) {
        int inod1 = (int) ista(iseg);
        int inod2 = (int) iend(iseg);
        nodtyp(inod1) += 1;
        nodtyp(inod2) += 1;
        if(nodtyp(inod1) > nodsegm) {
            printf("*** Error: Too many segments connected to node %i\n", inod1);
        }
        if(nodtyp(inod2) > nodsegm) {
            printf("*** Error: Too many segments connected to node %i\n", inod2);
        }
    }
    
    for (int inod = 0; inod < nnod; inod++) {
        if (nodtyp(inod) == 0)  {
            nodname.shed_row(inod);
            nodtyp.shed_row(inod);
            cnode.shed_col(inod);
            inod = 0;
            nnod -= 1;
        }
    }
    
    ivec copyBCnodname = bcnodname;
    ivec copyBCtyp = bctyp;
    vec copyBCprfl = bcprfl;
    
    nnodbc = 0;
    for (int inod = 0; inod < nnod; inod++) {
        if (nodtyp(inod) == 1)  {
            nnodbc += 1;
        }
    }
    
    bcnodname = zeros<ivec>(nnodbc);
    bctyp = zeros<ivec>(nnodbc);
    bcprfl = zeros<vec>(nnodbc);
    bchd = zeros<vec>(nnodbc);
    bchd.fill(0.4);
    
    int jnodbc = 0;
    for (int inod = 0; inod < nnod; inod++) {
        if (nodtyp(inod) == 1)  {
            bcnodname(jnodbc) = nodname(inod);
            for (int knodbc = 0; knodbc < copyBCnodname.n_elem; knodbc++)   {
                if (bcnodname(jnodbc) == copyBCnodname(knodbc)) {
                    bctyp(jnodbc) = copyBCtyp(knodbc);
                    bcprfl(jnodbc) = copyBCprfl(knodbc);
                    goto here;
                }
                else {
                    bctyp(jnodbc) = 3;
                }
            }
        here:;
            jnodbc += 1;
        }
    }
    
    printf("Subnetwork: No. of segments = %i\n",nseg);
    printf("Subnetwork: No. of nodes = %i\n",nnod);
    printf("Subnetwork: No. of boundary nodes = %i\n",nnodbc);
}
