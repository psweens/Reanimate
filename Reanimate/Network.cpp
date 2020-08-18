#include "Network.hpp"

#include <string>
#include <dirent.h>

using namespace std;
using namespace reanimate;

void Network::setBuildPath() {

    printText("Setting build directory");
    DIR* dir = opendir(buildPath.c_str());
    if (dir) {
        printText("Directory exits, removing contents",2,0);
        string contents = buildPath + "*";
        system(("rm " + contents).c_str());
        closedir(dir);
    } else if (ENOENT == errno) {
        printText("Directory does not exist, creating folder",2,0);
        system(("mkdir " + buildPath).c_str());
    }

}

Network::Network() {

    gamma = 1./(1.e3*60);
    alpha = 0.1333;
    beta = 1.e-4;
    xi = 0.001*1e-3;

    consthd = 0.45;

    kp = 0.1;
    ktau = 1.e-4;

    unknownBCs=false;

    rLog = "Reanimate_Log.txt";
    initLog();

}
Network::~Network() = default;


void Network::loadNetwork(const string &filename, const bool directFromAmira)   {
    
    int max=200;
    char bb[200];

    if (!directFromAmira)   {networkPath = loadPath + filename;}
    else {networkPath = buildPath + filename;}

    FILE *ifp;
    ifp = fopen(networkPath.c_str(),"r");

    fgets(bb,max,ifp);
    bb[strcspn(bb, "\n")] = 0;
    networkName = bb;

    printText(networkName, 3);
    printText("Importing network data");
    
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
        printText("Network File -> Invalid Segment Format",4);
    }
    qq = abs(q);


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
        printText("Network file -> Invalid Node Format",4);
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
        nsol -= 4; // Count extra solutes
        bcp = zeros<mat>(nnodbc,nsol);
        for(int inodbc = 0; inodbc < nnodbc; inodbc++){
            fscanf(ifp,"%lli %lli %lf %lf", &bcnodname(inodbc),&bctyp(inodbc),&bcprfl(inodbc),&bchd(inodbc));
            for (int i = 0; i < nsol; i++)  {
                fscanf(ifp," %lf",&bcp(inodbc,i));
            }
            fscanf(ifp,"\n");
        }
    }
    else    {
        printText("Network file -> Invalid Boundary Node Format",4);
    }
    
    fclose(ifp);
    
    nodsegm += 1; // Armadillo indexing starts at zero

    setup_networkArrays();
    analyse_network();
    edgeNetwork();
    printNetwork("loadedNetworkFile.txt");

    field<string> headers(2,1);
    headers(0,0) = "Diameter";
    headers(1,0) = "Length";
    mat data = zeros<mat>(nseg,2);
    data.col(0) = diam;
    data.col(1) = lseg;
    printHistogram(buildPath + "DiamLength_HistogramData.txt", data, headers);

}

void Network::analyse_network(bool graph)   {

    if (graph)  {
        printNum("No. of edges =", nseg);
        printNum("No. of vertices =", nnod);
        printNum("No. of leaf vertices =", nnodbc);
    }
    else {
        printNum("No. of segments =", nseg);
        printNum("No. of nodes =", nnod);
        printNum("No. of boundary nodes =", nnodbc);
    }


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
            printText( "No matching node found for segname " + to_string(segname(iseg)),4);
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
            printText( "Too many segments connected to node " + to_string(inod1),4);
        }
        if (nodtyp(inod2) > nodsegm) {
            printText( "Too many segments connected to node " + to_string(inod2),4);
        }
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
                    printText( "Boundary node " + to_string(nodname(inod)) + " is not a 1-segment node",4);
                }
                goto foundit2;
            }
        }
        printText("No matching node found for nodname " + to_string(bcnodname(inodbc)) + " , " + to_string(inodbc), 4);
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

void Network::subNetwork(ivec &index, bool graph) {

    printText("Creating subnetwork");
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
            printText( "No matching node found for segname " + to_string(segname(iseg)),4);
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
            printText( "Too many segments connected to node " + to_string(inod1),4);
        }
        if(nodtyp(inod2) > nodsegm) {
            printText( "Too many segments connected to node " + to_string(inod2),4);
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
            for (int knodbc = 0; knodbc < (int) copyBCnodname.n_elem; knodbc++)   {
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

    (*this).setup_networkArrays();
    (*this).analyse_network(graph);

}
