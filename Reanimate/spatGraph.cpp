#include "spatGraph.hpp"
#include <fstream>

using namespace reanimate;
using namespace std;

spatGraph::spatGraph() {

    Dmin = 16.5;
    Dr = 1.25;

}
spatGraph::~spatGraph() = default;

void spatGraph::setup_graphArrays() {

    //q = zeros<vec>(nseg);
    rseg = zeros<vec>(nseg);
    //vesstyp = zeros<ivec>(nseg);
    ista = zeros<ivec>(nseg);
    iend = zeros<ivec>(nseg);
    geometry = zeros<ivec>(nseg);
    isBridge = zeros<ivec>(nseg);
    isBridgehead = zeros<ivec>(nseg);
    subGraphs = zeros<ivec>(nseg);

    nodout = zeros<ivec>(nnod);
    nodrank = zeros<ivec>(nnod);
    nk = zeros<ivec>(nnod);
    nodpress = zeros<vec>(nnod);
    nodtyp = zeros<ivec>(nnod);
    articPnt = zeros<ivec>(nnod);

    bcnod = zeros<ivec>(nnodbc);

    nodnod = zeros<imat>(nodsegm,nnod);
    nodseg = zeros<imat>(nodsegm,nnod);

    segnod = zeros<sp_mat>(nnod,nnod);

}


//template <class CallNetwork>
void spatGraph::generate(Network &network, bool print)  {

    buildPath = network.buildPath;
    loadPath = network.loadPath;

    printText("Generating spatial graph");

    nodtyp = zeros<ivec>(network.getNnod());

    uvec segIdx = find_unique(network.edgeLabels,false);
    nseg = segIdx.n_elem;
    segname = network.edgeLabels(segIdx);
    vesstyp = network.vesstyp(segIdx);
    diam = network.ediam(segIdx);
    lseg = network.elseg(segIdx);
    q = network.q(segIdx);
    hd = network.hd(segIdx);


    // Finding edge vertices
    ivec flagLoop = zeros<ivec>(nseg);
    segnodname = zeros<imat>(2,nseg);
    ngraphTag = zeros<ivec>(nnod);
    for (int jseg = 0; jseg < nseg; jseg++) {
        nodtyp.zeros();
        for (int iseg = 0; iseg < network.getNseg(); iseg++) {
            if (segname(jseg) == network.edgeLabels(iseg)) {
                nodtyp(network.ista(iseg)) += 1;
                nodtyp(network.iend(iseg)) += 1;
            }
        }
        if (min(nodtyp(find(nodtyp))) == 1)   {
            uvec idx = find(nodtyp == 1);
            segnodname(0, jseg) = network.nodname(idx(1));
            segnodname(1, jseg) = network.nodname(idx(0));
        }
    }


    // List nodes & coordinates
    ivec nodIdx = unique(segnodname);
    nnod = nodIdx.n_elem;
    nodname = zeros<ivec>(nnod);
    cnode = zeros<mat>(3,nnod);
    for (int inod = 0; inod < nnod; inod++) {
        for (int jnod = 0; jnod < (int) network.nodname.n_elem; jnod++) {
            if (nodIdx(inod) == network.nodname(jnod))  {
                nodname(inod) = network.nodname(jnod);
                cnode.col(inod) = network.cnode.col(jnod);
                jnod = network.nodname.n_elem;
            }
        }
    }


    bcnodname = network.bcnodname;
    bctyp = network.bctyp;
    bcprfl = network.bcprfl;
    bchd = network.bchd;
    nnodbc = bcnodname.n_elem;

    nodsegm = network.nodsegm;
    ista = zeros<ivec>(nseg);
    iend = zeros<ivec>(nseg);
    indexSegmentConnectivity();

    alx = network.alx;
    aly = network.aly;
    alz = network.alz;
    mxx = network.mxx;
    myy = network.myy;
    mzz = network.mzz;
    lb = network.lb;
    maxl = network.maxl;


    setup_graphArrays();
    setup_networkArrays();
    analyse_network(true, print);
    printNetwork("spatialGraph.txt");

    int cntr = 0;
    for (int inodbc = 0; inodbc < network.getNnodbc(); inodbc++)    {
        if (network.bcprfl(inodbc) == 0.0 && network.bctyp(inodbc) == 1)    {
            cntr += 1;
        }
    }
    cntr -= network.getNnodbc();
    if (nnodbc != abs(cntr)) {printText("Incorrect number of boundary nodes. Target no. = "+to_string(abs(cntr)),5);}
    for (int inodbc = 0; inodbc < nnodbc; inodbc++) {
        uvec idx = find(bcnodname(inodbc) == network.bcnodname);
        if (idx.n_elem == 0)    {printText("New boundary node detected",5);}
    }

    for (int iseg = 0; iseg < nseg; iseg++) {
        int inod1 = (int) ista(iseg);
        int inod2 = (int) iend(iseg);
        segnod(inod1,inod2) = iseg;
        segnod(inod2,inod1) = iseg;
    }

}


void spatGraph::linkEdges()  {

    for (int inod = 0; inod < nnod; inod++) {
        if (nodtyp(inod) == 2)  {
            // Find connected edges between two bifurcations
            ivec tag = -ones<ivec>(nnod);
            ivec segFlag = zeros<ivec>(nseg);
            dfsBasic(inod,1,tag, true, 2);

            for (int iseg = 0; iseg < nseg; iseg++) {
                if (tag(ista(iseg)) == 1 || tag(iend(iseg)) == 1) {segFlag(iseg) = 1; }
            }

            diam(find(segFlag == 1)).fill(mean(diam(find(segFlag == 1))));
            lseg(find(segFlag == 1)).fill(accu(lseg(find(segFlag == 1))));
        }
    }

}


void spatGraph::defineTrunk()   {

    cout<<"Running network topology classifier ..."<<endl;

    nTrees = 0;
    cout<<"How many feeding/draining vessels does the network contain? ";
    cin >> nTrees;

    int classify = 0;
    int nodInput = 0;
    InOutlets = zeros<imat>(nTrees,4);
    for (int ives = 0; ives < nTrees; ives++) {

        nodnameRewind:;
        cout<<"Please enter node name: ";
        cin>>nodInput;
        // Node name checker
        int nodCheck = 0;
        int nodindx = 0;
        for (int inod = 0; inod < nnod; inod++) {
            if (nodname(inod) == nodInput)  {
                nodCheck += 1;
                nodindx = inod;
            }
        }
        if (nodCheck == 0)  {
            cout<<"*** Error: Non-existent nodname name ***"<<endl;
            goto nodnameRewind;
        }
        cout<<"Feeding (1) or draining (2) vessel? ";
        cin>>classify;
        if (classify == 1)  {
            narterioles += 1;
        }
        else if (classify == 2) {
            nvenules += 1;
        }
        for (int iseg = 0; iseg < nseg; iseg++)    {
            if (nodname(ista(iseg)) == nodInput) {
                InOutlets(ives,0) = iseg;
                InOutlets(ives,1) = classify;
                InOutlets(ives,2) = diam(iseg);
                InOutlets(ives,3) = nodindx;
            }
            else if (nodname(iend(iseg)) == nodInput)   {
                InOutlets(ives,0) = iseg;
                InOutlets(ives,1) = classify;
                InOutlets(ives,2) = diam(iseg);
                InOutlets(ives,3) = nodindx;
            }
        }
    }

    classifyNetwork(InOutlets, geometry);

}

void spatGraph::analyseTopology(imat predefinedInput)   {

    // predefinedInput: col(0) nodname, col(1) arteriole/venule; rows() trunks
    nTrees = predefinedInput.n_rows;
    InOutlets = zeros<imat>(nTrees,4);
    for (int i = 0; i < nTrees; i++)    {
        for (int iseg = 0; iseg < nseg; iseg++)    {
            if (nodname(ista(iseg)) == predefinedInput(i,0)) {
                InOutlets(i,0) = iseg;
                InOutlets(i,1) = predefinedInput(i,1);
                InOutlets(i,2) = diam(iseg);
                InOutlets(i,3) = ista(iseg);
            }
            else if (nodname(iend(iseg)) == predefinedInput(i,0))   {
                InOutlets(i,0) = iseg;
                InOutlets(i,1) = predefinedInput(i,1);
                InOutlets(i,2) = diam(iseg);
                InOutlets(i,3) = iend(iseg);
            }
        }
    }


    classifyNetwork(InOutlets, geometry);

}

void spatGraph::loadTrunks(const string &filepath) {

    cout<<"Loading network trunks ..."<<endl;

    int n{}, max{200};
    char bb[200];

    FILE *data;
    data = fopen(filepath.c_str(),"r");
    fscanf(data,"%i", &n); fgets(bb,max,data);
    fgets(bb,max,data);

    imat array = zeros<imat>(n,2);
    for (int i = 0; i < n; i++) {
        fscanf(data, "%lli %lli", &array(i,0), &array(i,1));
    }

    fclose(data);

    analyseTopology(array);

}