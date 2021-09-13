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
    uvec idx;
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
            idx = find(nodtyp == 1);
            segnodname(0, jseg) = network.nodname(idx(1));
            segnodname(1, jseg) = network.nodname(idx(0));
        }
        else {
            printText("Loop detected for edge "+to_string(segname(jseg)),5);
            idx = find(nodtyp > 0);
            for (int inod = 0; inod < (int) idx.n_elem; inod++)    {
                if (nodtyp(idx(inod)) != network.nodtyp(idx(inod))) {
                    segnodname(0, jseg) = network.nodname(idx(inod));
                    segnodname(1, jseg) = network.nodname(idx(inod));
                }
            }
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

    // Find and store vertex indices in full network
    indexVertices(network);

    // Find points per edge and label vertices
    indexPoints(network);

}


void spatGraph::indexVertices(Network &network)    {

    int nod1{},nod2{},found{};
    segvert = zeros<imat>(2, nseg);
    for (int iseg = 0; iseg < nseg; iseg++) {
        nod1 = segnodname(0,iseg);
        nod2 = segnodname(1,iseg);
        if (nod1 == nod2)   {
            // Loops
            for (int inod = 0; inod < network.getNnod(); inod++)    {
                if (nod1 == network.nodname(inod))  {
                    segvert(0,iseg) = inod;
                    segvert(1,iseg) = inod;
                }
            }
        }
        else {
            found = 0;
            for (int inod = 0; inod < network.getNnod(); inod++)    {
                if (nod1 == network.nodname(inod))   {
                    segvert(0,iseg) = inod;
                    found += 1;
                    if (found == 2) {inod = nnod;}
                }
                else if (nod2 == network.nodname(inod))  {
                    segvert(1,iseg) = inod;
                    found += 1;
                    if (found == 2) {inod = nnod;}
                }
            }
            if (found == 0) {printText("Vertex not found for edge " + to_string(segname(iseg)),4);}
        }
    }

}


void spatGraph::indexPoints(Network &network)    {

    uvec idx;
    edgePnts = zeros<ivec>(nseg);
    for (int iseg = 0; iseg < nseg; iseg++)  {
        idx = find(segname(iseg) == network.edgeLabels);
        edgePnts(iseg) = (int) idx.n_elem + 1;
    }
    network.npoint = accu(edgePnts);

    // Log full network node indices in order from start to end vertices
    int nod1,nod2{},enod{},pnt{},npnts{},maxpnts=(int) max(edgePnts);
    ivec tag{};
    segpoints = zeros<imat>(maxpnts, nseg);
    for (int iseg = 0; iseg < nseg; iseg++) {
        npnts = edgePnts(iseg);
        if (npnts == 2)    {
            // Loops
            segpoints(0,iseg) = segvert(0, iseg);
            segpoints(1, iseg) = segvert(1, iseg);
        }
        else {
            nod2 = -1;
            idx = find(segname(iseg) == network.edgeLabels);
            tag = zeros<ivec>((int) idx.n_elem);

            pnt = 0;
            nod1 = segvert(0, iseg);
            enod = segvert(1, iseg);
            segpoints(pnt, iseg) = nod1;
            pnt += 1;
            while (nod2 != enod) {
                for (int jseg = 0; jseg < (int) idx.n_elem; jseg++) {
                    if (network.ista(idx(jseg)) == nod1 && tag(jseg) == 0) {
                        nod2 = network.iend(idx(jseg));
                        nod1 = nod2;
                        tag(jseg) = 1;
                        segpoints(pnt, iseg) = nod1;
                        pnt += 1;
                        jseg = idx.n_elem;
                    }
                    else if (network.iend(idx(jseg)) == nod1 && tag(jseg) == 0) {
                        nod2 = network.ista(idx(jseg));
                        nod1 = nod2;
                        tag(jseg) = 1;
                        segpoints(pnt, iseg) = nod1;
                        pnt += 1;
                        jseg = idx.n_elem;
                    }
                }
            }
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

void spatGraph::analyseTopology(imat predefinedInput, Network &network)   {

    printText("Vessel Classification Module", 3);
    printText("Loading network trunks");

    popInletsOutlets(predefinedInput);

    classifyNetwork(InOutlets, geometry);

    uvec art = find(geometry == 1);
    uvec cap = find(geometry == 2);
    uvec ven = find(geometry == 3);
    printNum("Arterioles (%) =", 100*double(art.n_elem) / double(nseg));
    printNum("Capillaries (%) =", 100*double(cap.n_elem) / double(nseg));
    printNum("Venules (%) =", 100*double(ven.n_elem) / double(nseg));
    vesstyp = geometry;

    for (int iseg = 0; iseg < network.getNseg(); iseg++)    {
        for (int jseg = 0; jseg < nseg; jseg++) {
            if (network.edgeLabels(iseg) == segname(jseg))  {
                network.vesstyp(jseg) = vesstyp(jseg);
                jseg = nseg;
            }
        }
    }

}

void spatGraph::loadTrunks(const string &filepath, Network &network) {

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

    analyseTopology(array, network);

}

void spatGraph::findTree(imat input)  {

    popInletsOutlets(input);

    ivec nodvtyp = zeros<ivec>(nnod);
    for (int iseg = 0; iseg < nseg; iseg++) {
        if (vesstyp(iseg) == 1 || vesstyp(iseg) == 3)   {
            nodvtyp(ista(iseg)) = -1;
            nodvtyp(iend(iseg)) = -1;

        }
    }

    ivec oldnodtyp = nodvtyp;
    flagTree = zeros<ivec>(nseg);
    for (int i = 0; i < (int) input.n_rows; i++)    {
        for (int inodbc = 0; inodbc < nnodbc; inodbc++) {
            if (bcnodname(inodbc) == input(i,0))    {
                dfsBasic(inodbc, i+1, nodvtyp);
                inodbc = nnodbc;
            }
        }
        for (int iseg = 0; iseg < nseg; iseg++) {
            if (nodvtyp(ista(iseg)) == i+1 && nodvtyp(iend(iseg)) == i+1)   {
                flagTree(iseg) = i+1;
            }
        }
        nodvtyp = oldnodtyp;
    }

    geometry = vesstyp;

}

void spatGraph::popInletsOutlets(imat &predefinedInput)    {

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

}