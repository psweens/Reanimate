#include "Network.hpp"

using namespace std;
using namespace reanimate;

void Network::dfsArtic(int v, int p) {

    visited(v) = 1;
    tin(v) = timer++;
    low(v) = tin(v);
    int children{};
    int to{};
    int n = nodtyp(v);
    for (int inod = 0; inod < n; inod++)    {
        to = nodnod(inod,v);
        if (to == p) {continue;}
        if (visited(to) == 1) {low(v) = min(low(v), tin(to));}
        else {
            dfsArtic(to, v);
            low(v) = min(low(v), low(to));
            if (low(to) >= tin(v) && p != -1)   {articPnt(v) = 1;}
            ++children;
        }
    }
    if(p == -1 && children > 1) {articPnt(v) = 1;}

}

void Network::findArticulationPoints() {

    printText("Finding articulation points");
    timer = 0;
    visited = zeros<ivec>(nnod);
    tin = -ones<ivec>(nnod);
    low = -ones<ivec>(nnod);
    for (int i = 0; i < nnod; i++) {
        if (visited(i) == 0)
            dfsArtic(i);
    }
    printNum("No. of articulation points =", accu(articPnt));
    printNum("% of total nodes =", 100*accu(articPnt) / float(nnod));

}

ivec Network::findDeadends() {

    findArticulationPoints();
    articPnt(find(articPnt == 1 && nodtyp < 2)).fill(0); // Target aritculation points w/ nodtyp > 2

    printText("Finding dead ends");
    ivec nodeTags = zeros<ivec>(nnod);
    ivec deadEnds = zeros<ivec>(nseg);
    ivec tag = -ones<ivec>(nnod);
    uvec idx;
    int cntr{};
    for (int inod = 0; inod < nnod; inod++) {
        if (articPnt(inod) == 1 && nodeTags(inod) == 0)    {
            // Cycle through connections to node using dfs
            for (int jnod = 0; jnod < (int) nodtyp(inod); jnod++) {
                tag = -tag.ones();
                tag(inod) = 1;
                idx = find(nodeTags == 1);
                tag(idx).fill(1);
                dfsBasic(nodnod(jnod, inod),1,tag);
                tag(idx).fill(0);
                tag(inod) = 0;
                // Check if flowing boundary condition exists
                cntr = 0;
                for (int inodbc = 0; inodbc < nnodbc; inodbc++) {
                    if (tag(bcnod(inodbc)) == 1)   {
                        if (bctyp(inodbc) == 1 && bcprfl(inodbc) == 0.0)    {}
                        else {
                            cntr = 1;
                            inodbc = nnodbc;
                        }
                    }
                }
                if (cntr == 0)  {nodeTags(find(tag == 1)).fill(1);}
            }
            nodeTags(inod) = 2;
        }
    }

    for (int iseg = 0; iseg < nseg; iseg++) {
        if (nodeTags(ista(iseg)) == 1 || nodeTags(iend(iseg)) == 1) {deadEnds(iseg) = 1;}
    }
    printNum("No. of dead end segments =", accu(deadEnds));
    printNum("% of total segments =", 100*accu(deadEnds) / float(nseg));

    return deadEnds;

}