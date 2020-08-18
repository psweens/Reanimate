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
    isBridge = zeros<ivec>(nseg);
    subGraphs = zeros<ivec>(nseg);

    nodout = zeros<ivec>(nnod);
    nodrank = zeros<ivec>(nnod);
    nk = zeros<ivec>(nnod);
    nodpress = zeros<vec>(nnod);
    nodtyp = zeros<ivec>(nnod);

    bcnod = zeros<ivec>(nnodbc);

    nodnod = zeros<imat>(nodsegm,nnod);
    nodseg = zeros<imat>(nodsegm,nnod);
    segnod = zeros<imat>(nnod,nnod);

}

//template <class CallNetwork>
void spatGraph::generate(Network network, bool setNoflow)  {

    printText("Generating spatial graph");
    buildPath = network.buildPath;
    loadPath = network.loadPath;

    nodtyp = zeros<ivec>(network.getNnod());
    if (setNoflow)  {

        ivec fnodtyp = zeros<ivec>(network.getNnod());
        for (int i = 1; i <= max(network.edgeLabels); i++)  {
            nodtyp.zeros();
            for (int iseg = 0; iseg < network.getNseg(); iseg++) {
                int nod1 = network.ista(iseg);
                int nod2 = network.iend(iseg);
                if (network.edgeLabels(iseg) == i) {
                    nodtyp(nod1) += 1;
                    nodtyp(nod2) += 1;
                }
            }
            uvec idx = find(nodtyp == 1);
            uvec eidx = find(network.edgeLabels == i);
            if (idx.n_elem == 0)    {
                network.noflow(eidx).fill(1);
            }
            else if (network.nodpress(idx(0)) == network.nodpress(idx(1)))  {
                network.noflow(eidx).fill(1);
            }

        }

        for (int iseg = 0; iseg < network.getNseg(); iseg++) {
            int nod1 = network.ista(iseg);
            int nod2 = network.iend(iseg);
            if (network.noflow(iseg) == 1)  {
                fnodtyp(nod1) += 1;
                fnodtyp(nod2) += 1;
            }
            if (network.noflow(iseg) == 1)  {
                network.noflow(find(network.edgeLabels == network.edgeLabels(iseg))).fill(1);
            }
        }
        for (int iseg = 0; iseg < network.getNseg(); iseg++)    {
            int nod1 = network.ista(iseg);
            int nod2 = network.iend(iseg);
            if (fnodtyp(nod1) + 1 == nodtyp(nod1) && fnodtyp(nod2) + 1 == nodtyp(nod2)) {
                network.noflow(find(network.edgeLabels == network.edgeLabels(iseg))).fill(1);
            }
        }

    }


    //segname = unique(network.edgeLabels);
    uvec segIdx = find_unique(network.edgeLabels,false);
    if (setNoflow)  {segIdx = segIdx(find(network.noflow(segIdx) == 0));}
    nseg = segIdx.n_elem;
    q = network.q(segIdx);
    hd = network.hd(segIdx);
    segname = network.edgeLabels(segIdx);
    diam = network.ediam(segIdx);
    lseg = network.elseg(segIdx);
    vesstyp = network.vesstyp(segIdx);


    // Finding edge ends
    ivec ring = zeros<ivec>(nseg);
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
        uvec idx = find(nodtyp == 1);
        if (idx.n_elem == 0)    { // Ring detected -> remove
            ring(jseg) = 1;
        }
        else {
            segnodname(0, jseg) = network.nodname(idx(1));
            segnodname(1, jseg) = network.nodname(idx(0));
        }
    }
    if (accu(ring) > 0) {
        for (int iseg = 0; iseg < nseg; iseg++) {
            if (ring(iseg) == 1)    {
                segname.shed_row(iseg);
                vesstyp.shed_row(iseg);
                segnodname.shed_col(iseg);
                diam.shed_row(iseg);
                lseg.shed_row(iseg);
                q.shed_row(iseg);
                hd.shed_row(iseg);
                ring.shed_row(iseg);
                nseg -= 1;
                iseg = 0;
            }
        }
    }


    // List nodes & coordinates
    ivec nodIdx = unique(segnodname);
    nnod = nodIdx.n_elem;
    nodname = zeros<ivec>(nnod);
    cnode = zeros<mat>(3,nnod);
    for (int inod = 0; inod < nnod; inod++) {
        for (int jnod = 0; jnod < network.nodname.n_elem; jnod++) {
            if (nodIdx(inod) == network.nodname(jnod))  {
                nodname(inod) = network.nodname(jnod);
                cnode.col(inod) = network.cnode.col(jnod);
                jnod = network.nodname.n_elem;
            }
        }
    }

    ivec noflowBC = zeros<ivec>(network.getNnodbc());
    if (setNoflow)  {

        if (accu(network.noflow) > 0)   {
            for (int inodbc = 0; inodbc < network.getNnodbc(); inodbc++)    {
                for (int iseg = 0; iseg < network.getNseg(); iseg++)    {
                    if (network.noflow(iseg) == 1 && (network.bcnod(inodbc) == network.ista(iseg) || network.bcnod(inodbc) == network.iend(iseg))) {
                        noflowBC(inodbc) = 1;
                        iseg = network.getNseg();
                    }
                }
            }
        }

    }


    bcnodname = network.bcnodname(find(noflowBC == 0));
    nnodbc = bcnodname.n_elem;
    bctyp = network.bctyp(find(noflowBC == 0));
    bcprfl = network.bcprfl(find(noflowBC == 0));
    bchd = network.bchd(find(noflowBC == 0));

    ista = zeros<ivec>(nseg);
    iend = zeros<ivec>(nseg);
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
            printText( "No matching node found for segname "+to_string(segname(iseg)),4);
            foundit:;
        }
    }

    nodtyp = zeros<ivec>(nnod);
    for (int iseg = 0; iseg < nseg; iseg++) {
        int inod1 = (int) ista(iseg);
        int inod2 = (int) iend(iseg);
        nodtyp(inod1) += 1;
        nodtyp(inod2) += 1;
    }
    for (int inod = 0; inod < nnod; inod++) {
        uvec idx = find(nodname(inod) == bcnodname);
        if (nodtyp(inod) == 1 && idx.n_elem == 0)  {
            bcnodname.insert_rows(nnodbc,1);
            bcnodname(nnodbc) = nodname(inod);
            bctyp.insert_rows(nnodbc,1);
            bctyp(nnodbc) = 1;
            bcprfl.insert_rows(nnodbc,1);
            bcprfl(nnodbc) = 0.0;
            bchd.insert_rows(nnodbc,1);
            for (int iseg = 0; iseg < network.getNseg(); iseg++) {
                if (network.segnodname(0,iseg) == nodname(inod) || network.segnodname(1,iseg) == nodname(inod)) {
                    bchd(nnodbc) = network.hd(iseg);
                    iseg = network.getNseg();
                }
            }
            //bchd(nnodbc) = consthd;
            nnodbc += 1;
        }
    }

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
    analyse_network(true);

    for (int iseg = 0; iseg < nseg; iseg++) {
        int inod1 = (int) ista(iseg);
        int inod2 = (int) iend(iseg);
        segnod(inod1,inod2) = iseg;
        segnod(inod2,inod1) = iseg;
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
        fscanf(data, "%i %i", &array(i,0), &array(i,1));
    }

    fclose(data);

    analyseTopology(array);

}