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

    int nod1{},nod2{};
    for (int iseg = 0; iseg < nseg; iseg++) {
        nod1 = (int) ista(iseg);
        nod2 = (int) iend(iseg);
        segnod(nod1,nod2) = iseg;
        segnod(nod2,nod1) = iseg;
    }

    // Store Vertex indices in full network graph & index of a single segment from edge
    uvec idx,node1,node2;
    edgeSta = zeros<ivec>(nseg);
    edgeEnd = zeros<ivec>(nseg);
    edgeSeg = zeros<ivec>(nseg);
    int segmax{};
    for (int iseg = 0; iseg < nseg; iseg++) {
        idx = find(segname(iseg) == network.edgeLabels);
        if ((int) idx.n_elem > segmax)    {segmax = idx.n_elem;}
        if (idx.n_elem == 0)  {printText("No segments found in full network graph",5);}
        else {edgeSeg(iseg) = idx(0);}
        node1 = find(segnodname(0, iseg) == network.nodname);
        node2 = find(segnodname(1, iseg) == network.nodname);
        if (node1.n_elem > 1 || node2.n_elem > 1)  {printText("Multiple nodes found in full network graph for vertex",5);}
        else {
            edgeSta(iseg) = node1(0);
            edgeEnd(iseg) = node2(0);
        }
    }
    edgeseg = zeros<umat>(segmax,nseg);
    edgetyp = zeros<ivec>(nseg);
    for (int iseg = 0; iseg < nseg; iseg++) {
        idx = find(segname(iseg) == network.edgeLabels);
        edgetyp(iseg) = (int) idx.n_elem;
        for (int jseg = 0 ; jseg < edgetyp(iseg); jseg++) {
            edgeseg(jseg,iseg) = idx(jseg);
        }
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

    printText("Vessel Classification Module", 3);

    nTrees = 0;
    printText("How many feeding/draining vessels does the network contain?",0,1);
    cin >> nTrees;
    printNum("No. of trees =", nTrees);

    int classify = 0;
    int nodInput = 0;
    InOutlets = zeros<imat>(nTrees,5);
    for (int ives = 0; ives < nTrees; ives++) {

        nodnameRewind:;
        printText("Please enter node name",0,0);
        cin>>nodInput;
        printNum("Node name =", nodInput);
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
            printText("Non-existent nodname name",4);
            goto nodnameRewind;
        }
        printText("Feeding (1) or draining (2) vessel?",0,0);
        cin>>classify;
        printNum("Type =", classify);
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
                InOutlets(ives,4) = nodInput;
            }
            else if (nodname(iend(iseg)) == nodInput)   {
                InOutlets(ives,0) = iseg;
                InOutlets(ives,1) = classify;
                InOutlets(ives,2) = diam(iseg);
                InOutlets(ives,3) = nodindx;
                InOutlets(ives,4) = nodInput;
            }
        }
    }

    classifyNetwork(InOutlets, geometry);
    uvec art = find(geometry == 1);
    uvec cap = find(geometry == 2);
    uvec ven = find(geometry == 3);
    printNum("Arterioles (%) =", 100*double(art.n_elem) / double(nseg));
    printNum("Capillaries (%) =", 100*double(cap.n_elem) / double(nseg));
    printNum("Venules (%) =", 100*double(ven.n_elem) / double(nseg));

}

void spatGraph::analyseTopology(imat predefinedInput)   {

    printText("Vessel Classification Module", 3);
    printText("Loading network trunks");

    // predefinedInput: col(0) nodname, col(1) arteriole/venule; rows() trunks
    nTrees = predefinedInput.n_rows;
    InOutlets = zeros<imat>(nTrees,5);
    for (int i = 0; i < nTrees; i++)    {
        for (int iseg = 0; iseg < nseg; iseg++)    {
            if (nodname(ista(iseg)) == predefinedInput(i,0)) {
                InOutlets(i,0) = iseg;
                InOutlets(i,1) = predefinedInput(i,1);
                InOutlets(i,2) = diam(iseg);
                InOutlets(i,3) = ista(iseg);
                InOutlets(i,4) = predefinedInput(i,0);
            }
            else if (nodname(iend(iseg)) == predefinedInput(i,0))   {
                InOutlets(i,0) = iseg;
                InOutlets(i,1) = predefinedInput(i,1);
                InOutlets(i,2) = diam(iseg);
                InOutlets(i,3) = iend(iseg);
                InOutlets(i,4) = predefinedInput(i,0);
            }
        }
    }

    classifyNetwork(InOutlets, geometry);

    uvec art = find(geometry == 1);
    uvec cap = find(geometry == 2);
    uvec ven = find(geometry == 3);
    printNum("Arterioles (%) =", 100*double(art.n_elem) / double(nseg));
    printNum("Capillaries (%) =", 100*double(cap.n_elem) / double(nseg));
    printNum("Venules (%) =", 100*double(ven.n_elem) / double(nseg));
    vesstyp = geometry;

}

void spatGraph::loadTrunks(const string &filepath) {

    printText("Vessel Classification Module", 3);
    printText("Loading network trunks");

    int n{}, max{200};
    char bb[200];

    FILE *data;
    data = fopen((loadPath + filepath).c_str(),"r");
    fscanf(data,"%i", &n); fgets(bb,max,data);
    fgets(bb,max,data);

    imat array = zeros<imat>(n,2);
    for (int i = 0; i < n; i++) {
        fscanf(data, "%lli %lli", &array(i,0), &array(i,1));
    }

    fclose(data);

    analyseTopology(array);

}