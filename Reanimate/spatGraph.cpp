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

    q = zeros<vec>(nseg);
    rseg = zeros<vec>(nseg);
    vesstyp = zeros<ivec>(nseg);
    ista = zeros<ivec>(nseg);
    iend = zeros<ivec>(nseg);
    geometry = zeros<ivec>(nseg);

    nodout = zeros<ivec>(nnod);
    nodrank = zeros<ivec>(nnod);
    nk = zeros<ivec>(nnod);
    nodpress = zeros<vec>(nnod);
    nodtyp = zeros<ivec>(nnod);

    bcnod = zeros<ivec>(nnodbc);

    nodnod = zeros<imat>(nodsegm,nnod);
    nodseg = zeros<imat>(nodsegm,nnod);

}

void spatGraph::generate(Network &network)  {

    cout<<"\nGenerating spatial graph ..."<<endl;

    segname = unique(network.edgeLabels);
    nseg = segname.n_elem;

    uvec segIdx = find_unique(network.edgeLabels,false);
    diam = network.ediam(segIdx);
    lseg = network.elseg(segIdx);
    hd = network.hd;

    // Calculate no. of nodes
    nnod = network.getNnodbc();
    for (int inod = 0; inod < network.getNnod(); inod++)    {
        if (network.nodtyp(inod) > 2)   { // Counting vertices
            nnod += 1;
        }
    }

    nodtyp = zeros<ivec>(network.getNnod());
    segnodname = zeros<imat>(2,nseg);
    for (int jseg = 0; jseg < nseg; jseg++) {
        nodtyp.zeros();
        for (int iseg = 0; iseg < network.getNseg(); iseg++) {
            if (segname(jseg) == network.edgeLabels(iseg)) {
                nodtyp(network.ista(iseg)) += 1;
                nodtyp(network.iend(iseg)) += 1;
            }
        }
        // Find vertices
        uvec idx = find(nodtyp == 1 || nodtyp > 2);
        segnodname(0,jseg) = network.nodname(idx(0));
        segnodname(1,jseg) = network.nodname(idx(1));
    }


    // List nodes & coordinates
    ivec nodIdx = unique(segnodname);
    nodname = network.nodname(find(network.nodname == nodIdx));
    cnode = network.cnode.cols(find(network.nodname == nodIdx));

    nnodbc = network.getNnodbc();
    bcnodname = network.bcnodname;
    bctyp = network.bctyp;
    bcprfl = network.bcprfl;
    bchd = network.bchd;

    alx = network.alx;
    aly = network.aly;
    alz = network.alz;
    mxx = network.mxx;
    myy = network.myy;
    mzz = network.mzz;
    lb = network.lb;
    maxl = network.maxl;
    nodsegm = network.nodsegm;

    setup_graphArrays();
    analyse_network();

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
        fscanf(data, "%i %i", &array(i,0), &array(i,1));
    }

    fclose(data);

    analyseTopology(array);

}