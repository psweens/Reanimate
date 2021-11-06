#include "spatGraph.hpp"

using namespace reanimate;
using namespace std;

void spatGraph::dfsBridge(int v, int p) {
    // Depth first search
    visited(v) = 1;
    tin(v) = timer++;
    low(v) = tin(v);
    int to{};
    int n = nodtyp(v);
    for (int inod = 0; inod < n; inod++)    {
        to = nodnod(inod,v);
        if (to == p) {continue;}
        if (visited(to) == 1) {low(v) = min(low(v), tin(to));}
        else {
            dfsBridge(to, v);
            low(v) = min(low(v), low(to));
            if (low(to) > tin(v))   {isBridge(segnod(v,to)) = 1;}
        }
    }

}

void spatGraph::findBridges() {

    timer = 0;
    visited = zeros<ivec>(nnod);
    tin = -ones<ivec>(nnod);
    low = -ones<ivec>(nnod);
    for (int i = 0; i < nnod; i++) {
        if (visited(i) == 0)
            dfsBridge(i);
    }
    printNum("No. of bridges = ", accu(isBridge));

}


void spatGraph::findBridgeheads()   {

    printText("Finding graph bridgeheads");
    findBridges();

    // Traverse all boundary nodes to find bridgeheads (no flow boundaries removed during graph generation)
    ivec notBridgehead = zeros<ivec>(nnod);
    for (int iseg = 0; iseg < nseg; iseg++) {
        if (isBridge(iseg) == 1)    {
            if (nodtyp(ista(iseg)) <= 2) {notBridgehead(ista(iseg)) = -1;}
            if (nodtyp(iend(iseg)) <= 2) {notBridgehead(iend(iseg)) = -1;}
        }
    }

    // Link bridges to edges with a nodtyp=2 connection to bridge (edges in series)
    for (int inod = 0; inod < nnod; inod++) {
        if (notBridgehead(inod) == -1)  {
            ivec linkBridge = -ones<ivec>(nnod);
            dfsBasic(inod,1,linkBridge, true, 2);
            uvec idx  = find(linkBridge == 1);
            notBridgehead(idx).fill(-1);
        }
    }

    for (int iseg = 0; iseg < nseg; iseg++) {
        if (notBridgehead(ista(iseg)) == -1 || notBridgehead(iend(iseg)) == -1)  {
            isBridge(iseg) = 1;
        }
    }

    // Remove boundary edges from bridge list
    for (int inodbc = 0; inodbc < nnodbc; inodbc++) {dfsBasic(bcnod(inodbc),1,notBridgehead, true, 2);}
    isBridgehead = isBridge;
    for (int iseg = 0; iseg < nseg; iseg++) {
        if (notBridgehead(ista(iseg)) == 1 || notBridgehead(iend(iseg)) == 1)   {
            isBridgehead(iseg) = 0;
        }
    }

    mat ED = zeros<mat>(nseg, 1);
    ED.col(0) = conv_to<vec>::from(isBridgehead);
    printAmira("amiraBridgeheads4.am", ED);

    printNum("No. of bridgeheads = ", accu(isBridgehead));

}


ivec spatGraph::findParallelEdges() {

    ivec idx = zeros<ivec>(nseg);
    for (int iseg = 0; iseg < nseg; iseg++) {
        for (int jseg = iseg + 1; jseg < nseg; jseg++)  {
            if ((ista(iseg) == ista(jseg) && iend(iseg) == iend(jseg)) || (ista(iseg) == iend(jseg) && iend(iseg) == ista(jseg)))   {
                idx(iseg) = 1;
                idx(jseg) = 1;
            }
        }
    }
    return idx;

}

ivec spatGraph::findParallelLoops(ivec &x)    {

    ivec idx = zeros<ivec>(nseg);
    for (int iseg = 0; iseg < nseg; iseg++) {
        if (x(iseg) == 1)    {
            if (nodtyp(ista(iseg)) == 2 || nodtyp(iend(iseg)) == 2)    {idx(iseg) = 1;}
        }
    }
    return idx;

}