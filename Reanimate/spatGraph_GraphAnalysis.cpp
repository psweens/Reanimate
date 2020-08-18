// Based on algorithms from https://cp-algorithms.com/graph/

#include "spatGraph.hpp"

using namespace reanimate;
using namespace std;

void spatGraph::dfsBasic(int v, int tag) {

    int to{};
    int n = nodtyp(v);
    visited(v) = tag;
    for (int inod = 0; inod < n; inod++)    {
        to = nodnod(inod,v);
        if (visited(to) == -1)   {dfsBasic(to,tag);}
    }

}

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
        if (visited(to) == 1) {
            low(v) = min(low(v), tin(to));
        }
        else {
            dfsBridge(to, v);
            low(v) = min(low(v), low(to));
            if (low(to) > tin(v))   {
                isBridge(segnod(v,to)) = 1;
            }
        }
    }

}

void spatGraph::findBridges() {

    timer = 0;
    visited = zeros<ivec>(nnod);
    visited(find(nodtyp == 1)).fill(1);
    for (int inodbc = 0; inodbc < nnodbc; inodbc++) {
        if (bcprfl(inodbc) == 0. && bctyp(inodbc) == 1) {
            visited(bcnod(inodbc) ) = 0;
        }
    }
    tin = -ones<ivec>(nnod);
    low = -ones<ivec>(nnod);
    for (int i = 0; i < nnod; i++) {
        if (visited(i) == 0)
            dfsBridge(i);
    }

    cout<<"No. of bridges = "<<accu(isBridge)<<endl;

}

void spatGraph::findBridgeheads()   {

    printText("Finding graph bridgeheads");
    findBridges();

    // Find bridgeheads by removing boundary bridges
    isBridgehead = isBridge;
    for (int iseg = 0; iseg < nseg; iseg++) {
        if (isBridge(iseg) == 1)    {
            if (nodtyp(ista(iseg)) == 1 || nodtyp(iend(iseg)) == 1) {
                isBridgehead(iseg) = 0;
            }
        }
    }
    printNum("No. of bridgeheads =", accu(isBridgehead));


    ivec copyBridgehead = isBridgehead;
    spatGraph allGraphs = *this;
    allGraphs.subNetwork(copyBridgehead, true); // Remove bridges
    allGraphs.traverseGraph();

    for (int iseg = 0; iseg < nseg; iseg++) {
        if (isBridgehead(iseg) == 1)    {
            subGraphs(iseg) = -3;
        }
        else {
            int seg = (int) segname(iseg);
            for (int jseg = 0; jseg < allGraphs.getNseg(); jseg++) {
                if (seg == allGraphs.segname(jseg)) {
                    subGraphs(iseg) = allGraphs.subGraphs(jseg);
                }
            }
        }
    }

}

void spatGraph::traverseGraph()  {

    int ngraphs = 0; // No. of well-connected (i.e. subgraphs with flowing boundary vertices)
    visited = -ones<ivec>(nnod);
    for (int inodbc = 0; inodbc < nnodbc; inodbc++) {
        if (visited(bcnod(inodbc)) == -1) {
            if (bctyp(inodbc) != 1) {
                ngraphs += 1;
                dfsBasic(bcnod(inodbc), inodbc);
            }
        }
    }
    int disconnectedGraphs = 0;
    for (int inod = 0; inod < nnod; inod++) {
        if (visited(inod) == -1)    {
            disconnectedGraphs += 1;
            dfsBasic(inod, -2);
        }
    }

    printNum("No. of graphs with flowing leaf vertices =", ngraphs);
    printNum("No. of disconnected graphs =", disconnectedGraphs);

    for (int iseg = 0; iseg < nseg; iseg++) {
        if (visited(ista(iseg)) == visited(iend(iseg))) {
            subGraphs(iseg) = visited(ista(iseg));
        }
        else {printText("Edge connected to two subgraphs",4);}
    }

}