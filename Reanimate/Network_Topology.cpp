#include "Network.hpp"
#include <list>

using namespace std;
using namespace reanimate;

void Network::dfsBasic(int v, int tag, ivec &track, bool nodeCondition, int type) {

    int to{};
    int n = nodtyp(v);
    if (nodeCondition)  {
        if (n <= type)  {
            track(v) = tag;
            for (int i = 0; i < n; i++)    {
                to = nodnod(i, v);
                if (track(to) == -1)   {dfsBasic(to, tag,track, nodeCondition, type);}
            }
        }
    }
    else {
        track(v) = tag;
        for (int i = 0; i < n; i++)    {
            to = nodnod(i, v);
            if (track(to) == -1)   {dfsBasic(to, tag,track);}
        }
    }

}


void Network::doubleDfs(int v, int tag, double val, ivec &track, vec &param, string cond)  {

    // GT - greater than
    // ST - smaller than

    int to{};
    int n = nodtyp(v);
    double nod = param(v);
    if (cond == "GT")    {
        if (val <= nod)  {
            track(v) = tag;
            for (int i = 0; i < n; i++)    {
                to = nodnod(i, v);
                if (track(to) == -1)   {doubleDfs(to, tag, val, track, param, cond);}
            }
        }
    }
    else {
        if (nod <= val)  {
            track(v) = tag;
            for (int i = 0; i < n; i++)    {
                to = nodnod(i, v);
                if (track(to) == -1)   {doubleDfs(to, tag, val, track, param, cond);}
            }
        }
    }

}

void Network::dfsBranch(int v, int tag, ivec &track, int maxBranch)  {

    int to{};
    int n = nodtyp(v);
    if (n > 2)  {branch += 1;}
    if (branch <= maxBranch)  {
        track(v) = tag;
        for (int i = 0; i < n; i++)    {
            to = nodnod(i, v);
            if (abs(diam(nodseg(i, v))) > 9.)  {
                if (track(to) == -1)   {dfsBranch(to, tag, track, maxBranch);}
            }
        }
    }

}

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

// Breadth-first search which keeps track of nodal depth
ivec Network::breadthFirstSearch(int nod)  {

    // Mark all the vertices as not visited
    visited = -ones<ivec>(nnod);

    int depth{};
    ivec order = zeros<ivec>(nnod);
    list<int> bfsQueue;
    list<int> depthQueue;
    bfsQueue.push_back(nod);
    depthQueue.push_back(depth);
    visited(nod) = 1;
    while (!bfsQueue.empty()) {
        nod = bfsQueue.front();
        depth = depthQueue.front();
        order(nod) = depth;
        bfsQueue.pop_front();
        depthQueue.pop_front();

        int childNod{};
        for (int inod = 0; inod < nodtyp(nod); inod++)  {
            childNod = nodnod(inod, nod);
            if (visited(childNod) == -1)   {
                visited(childNod) = 1;
                bfsQueue.push_back(childNod);
                depthQueue.push_back(depth+1);
            }
        }
    }

    return order;

}


// Dijkstra's shortest path algorithm
ivec Network::findShortestPath(int startnode, int endNode, vec &edgeWeight) {

    int cntr{1},nextNode{};
    double infinity{1.e20},minDistance{};
    visited = -ones<ivec>(nnod);
    vec distance = zeros<vec>(nnod);
    vec pred = zeros<vec>(nnod);
    mat cost = infinity*ones<mat>(nnod, nnod);
    for (int iseg = 0; iseg < nseg; iseg++) {
        cost(ista(iseg), iend(iseg)) = edgeWeight(iseg);
        cost(iend(iseg), ista(iseg)) = edgeWeight(iseg);
    }

    for (int inod = 0; inod < nnod; inod++) {
        distance(inod) = cost(startnode, inod);
        pred(inod) = startnode;
    }
    distance(startnode) = 0.;
    visited(startnode) = 1;

    while (cntr < nnod-1) {
        minDistance = infinity;
        for (int inod = 0; inod < nnod; inod++) {
            if (distance(inod) < minDistance && visited(inod) == -1) {
                minDistance = distance(inod);
                nextNode = inod;
            }
        }

        visited(nextNode) = 1;
        for (int inod = 0; inod < nnod; inod++) {
            if (visited(inod) == -1) {
                if (minDistance + cost(nextNode, inod) < distance(inod)) {
                    distance(inod) = minDistance + cost(nextNode, inod);
                    pred(inod) = nextNode;
                }
            }
        }
        cntr += 1;
    }

    int knod{};
    for (int inod = 0; inod < nnod; inod++) {
        if (inod != startnode) {
            cout << "\nDistance of node " << inod << " = " << distance(inod);
            cout << "\nPath=" << inod;
            knod = inod;
            do {
                knod = pred(knod);
                cout << "<-" << knod;
            } while (knod != startnode);
        }
    }

    int nod{endNode};
    cntr = 0;
    if (nod != startnode) {
        cout << "\nDistance of node " << nod << " = " << distance(nod);
        cout << "\nPath=" << nodname(nod);
        knod = nod;
        cntr += 1;
        do {
            knod = pred(knod);
            cout << "<-" << nodname(knod);
            cntr += 1;
        } while (knod != startnode);
    }
    ivec path = zeros<ivec>(cntr);
    cntr = 0;
    if (nod != startnode) {
        path(cntr) = nod;
        cntr += 1;
        knod = nod;
        do {
            knod = pred(knod);
            path(cntr) = knod;
            cntr += 1;
        } while (knod != startnode);
    }

    return path;

}