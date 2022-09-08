#include "Network.hpp"
#include "omp.h"

using namespace reanimate;

void Network::analyse_network(bool graph, bool print)   {

    if (print)  {
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
    }

    // Populate ista/iend indices
    indexSegmentConnectivity();

    // Setup nodtyp, nodseg and nodnod
    indexNodeConnectivity();
    indexBCconnectivity();

    // Start(k,iseg) = coordinates of starting point of segment iseg
    // End(k,iseg) = coordinates of ending point of segment iseg
    computeLseg = 1;
    qq = abs(q);
    rseg = diam / 2.0;
    if (computeLseg == 1)   {findLengths();}

    BCgeo = zeros<ivec>(nnodbc);
    analyseBoundaryType();

/*    // Network volume
    double netVol = sum(lseg % pow(rseg, 2))*M_PI;

    // Vascular density
    double vascDens = 100*netVol/(alx*aly*alz);

    cout<<"Vascular Density: "<<vascDens<<endl;*/

}


void Network::indexSegmentConnectivity() {

    int seg{};
    bool found = false;
    ista = zeros<ivec>(nseg);
    iend = zeros<ivec>(nseg);
    #pragma omp parallel for schedule(auto) default(none) shared(nseg, nnod, segnodname, nodname, segname, ista, iend) private(seg, found)
    for (int iseg = 0; iseg < nseg; iseg++)	{
        //Search for nodes corresponding to this segment
        for (int i = 0; i < 2; i++) {
            seg = segnodname(i,iseg);
            for (int inod = 0; inod < nnod; inod++) {
                if (nodname(inod) == seg)    {
                    if (i == 0) {
                        ista(iseg) = inod;
                        inod = nnod;
                        found = true;
                    }
                    else if (i == 1)    {
                        iend(iseg) = inod;
                        inod = nnod;
                        found = true;
                    }
                }
            }
            if (!found) {printText("No matching node found for segname " + to_string(segname(iseg)),4);}
            else {found = false;}
        }
    }
}


void Network::indexNodeConnectivity()   {

    // Setup nodtyp, nodseg and nodnod
    // 'nodseg' -> for each nodes, store the corresponding segment index
    // 'nodnod' -> for each node, store the nodal index of the opposite side of the segment
    int inod1{},inod2{};
    nodtyp = zeros<ivec>(nnod);
    nodnod = zeros<imat>(nodsegm,nnod);
    nodseg = zeros<imat>(nodsegm,nnod);
    for (int iseg = 0; iseg < nseg; iseg++) {
        inod1 = (int) ista(iseg);
        inod2 = (int) iend(iseg);
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

}

void Network::indexBCconnectivity() {

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
        printText("No matching boundary node found for nodname " + to_string(bcnodname(inodbc)) + ", " + to_string(inodbc), 4);
        foundit2:;
    }

}

void Network::findLengths() {

    vec ss = zeros<vec>(3);
    mat Start = zeros<mat>(3,nseg);
    mat End = zeros<mat>(3,nseg);
    for (int iseg = 0; iseg < nseg; iseg++) {
        for (int k = 0; k < 3; k++){
            Start(k,iseg) = cnode(k,ista(iseg));
            End(k,iseg) = cnode(k,iend(iseg));
            ss(k) = End(k,iseg) - Start(k,iseg);
        }
        lseg(iseg) = sqrt(pow(ss(0),2) + pow(ss(1),2) + pow(ss(2),2));
    }

}


void Network::analyseBoundaryType()  {

    for (int inodbc = 0; inodbc < nnodbc; inodbc++) {
        if (vesstyp(nodseg(0,bcnod(inodbc))) == 1)  {BCgeo(inodbc) = 1;}
        else if (vesstyp(nodseg(0,bcnod(inodbc))) == 2)  {BCgeo(inodbc) = 2;}
        else if (vesstyp(nodseg(0,bcnod(inodbc))) == 3)  {BCgeo(inodbc) = 3;}
        //else {printText( "Boundary node "+to_string(bcnodname(inodbc))+" not classified",4);}
    }

}


void Network::scaleNetwork(double sfactor)    {

    diam *= sfactor;
    rseg *= sfactor;
    lseg *= sfactor;
    cnode *= sfactor;

}


void Network::edgeNetwork() {

    printText( "Analysing network edges");

    int edgIdx = 1;
    ivec flag = zeros<ivec>(nseg);
    edgeLabels = -ones<ivec>(nseg);
    elseg = zeros<vec>(nseg);
    ediam = zeros<vec>(nseg);

    int ntyp{};
    uvec bnod = find(nodtyp != 2); // Find index of boundary or bifurcating nodes
    ivec trackNode = -ones<ivec>(nnod);
    for (int inod = 0; inod < nnod; inod++) {
        if (nodtyp(inod) != 2)  {
            bool runme = true;
            if (runme)    {
                ntyp = (int) nodtyp(inod);
                // Cycle through connections to node using dfs
                for (int jnod = 0; jnod < ntyp; jnod++) {
                    if (trackNode(nodnod(jnod, inod)) == -1)    {
                        ivec tag = -ones<ivec>(nnod);
                        tag(inod) = 1;
                        dfsBasic(nodnod(jnod, inod),1,tag, true, 2);
                        tag(inod) = 0;
                        trackNode(find(tag == 1)).fill(1);
                        for (int iseg = 0; iseg < nseg; iseg++) {
                            if (tag(ista(iseg)) == 1 || tag(iend(iseg)) == 1) {
                                if (flag(iseg) == 0) {
                                    edgeLabels(iseg) = edgIdx;
                                    flag(iseg) = 1;
                                }
                            }
                        }
                        edgIdx += 1;
                    }
                }
            }
        }
    }
    // Find edges with a single segment connecting bifurcations
    for (int iseg = 0; iseg < nseg; iseg++) {
        if (flag(iseg) == 0)    {
            if (nodtyp(ista(iseg)) != 2 && nodtyp(iend(iseg)) != 2) {
                edgeLabels(iseg) = edgIdx;
                flag(iseg) = 1;
                trackNode(ista(iseg)) = 1;
                trackNode(iend(iseg)) = 1;
                edgIdx += 1;
            }
        }
    }


    for (int i = 0; i < (int) edgeLabels.n_elem; i++) {
        if (edgeLabels(i) <= 0) {
            printText( "Segment(s) not assigned an edge",4);
            i = edgeLabels.n_elem;
        }
    }

    ivec edgeIdx = unique(edgeLabels);
    for (int iseg = 0; iseg < (int) edgeIdx.n_elem; iseg++)  {
        uvec segs = find(edgeLabels == edgeIdx(iseg));
        vec test = diam(segs);
        double d = mean(diam(segs));
        double l = sum(lseg(segs));
        ediam(segs).fill(d);
        elseg(segs).fill(l);
    }
    nedge = edgeIdx.n_elem;
    uvec newIdx = find_unique(edgeLabels);
    printNum("No. of edges =", nedge);
    printStat("Edge lengths =", elseg(newIdx), "um");
    printStat("Edge diameters =", ediam(newIdx), "um");
    printStat("Length / diameter ratio =", elseg(newIdx)/ediam(newIdx), "um");

    nedge = elseg.n_rows;

}

double Network::nodeAverage(const int &inod, const vec &param)    {

    double var{};
    for (int jnod = 0; jnod < nodtyp(inod); jnod++) {
        var += param(nodseg(jnod,inod));
    }

    return var /= nodtyp(inod);

}
