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

    printText("Initiating bread first search");

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
ivec Network::findShortestPath(int startnode, int endNode, vec &edgeWeight, bool printPaths) {

    printText("Finding shortest path");

    int cntr{1},nextNode{},nod1{},nod2{};
    double infinity{1.e20},minDistance{};
    ivec mapped = -ones<ivec>(nnod);
    ivec pred = zeros<ivec>(nnod);
    vec distance = infinity*ones<vec>(nnod);

    for (int inod = 0; inod < nnod; inod++) {pred(inod) = startnode;}
    for (int iseg = 0; iseg < nseg; iseg++) {
        nod1 = ista(iseg);
        nod2 = iend(iseg);
        if (nod1 == startnode)    {
            distance(nod2) = edgeWeight(iseg);
        }
        else if (nod2 == startnode)   {
            distance(nod1) = edgeWeight(iseg);
        }
    }
    distance(startnode) = 0.;
    mapped(startnode) = 1;

    while (cntr < nnod-1) {
        minDistance = infinity;
        for (int inod = 0; inod < nnod; inod++) {
            if (distance(inod) < minDistance && mapped(inod) == -1) {
                minDistance = distance(inod);
                nextNode = inod;
            }
        }

        mapped(nextNode) = 1;
        for (int iseg = 0; iseg < nseg; iseg++) {
            nod1 = ista(iseg);
            nod2 = iend(iseg);
            if (nextNode == nod1 && mapped(nod2) == -1)   {
                if (minDistance + edgeWeight(iseg) < distance(nod2))    {
                    distance(nod2) = minDistance + edgeWeight(iseg);
                    pred(nod2) = nextNode;
                }
            }
            else if (nextNode == nod2 && mapped(nod1) == -1)   {
                if (minDistance + edgeWeight(iseg) < distance(nod1))    {
                    distance(nod1) = minDistance + edgeWeight(iseg);
                    pred(nod1) = nextNode;
                }
            }
        }
        cntr += 1;
    }

    if (printPaths == true) {printShortestPaths(startnode, endNode, distance, pred);}

    printText("Shortest path found. Populating route",2,0);
    cntr = 0;
    int knod{},nod{endNode};
    ivec path = zeros<ivec>(nnod);
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
    if (cntr < nnod)    {path.shed_rows(cntr,nnod-1);}

    return path;

}

void Network::printShortestPaths(int startNode, int endNode, vec distance, ivec pred)   {

    int knod{};
    for (int inod = 0; inod < nnod; inod++) {
        if (inod != startNode) {
            cout << "\nDistance of node " << inod << " = " << distance(inod);
            cout << "\nPath=" << inod;
            knod = inod;
            do {
                knod = pred(knod);
                cout << "<-" << knod;
            } while (knod != startNode);
        }
    }

    int cntr{};
    int nod{endNode};
    if (nod != startNode) {
        cout << "\nDistance of node " << nod << " = " << distance(nod);
        cout << "\nPath=" << nodname(nod);
        knod = nod;
        cntr += 1;
        do {
            knod = pred(knod);
            cout << "<-" << nodname(knod);
            cntr += 1;
        } while (knod != startNode);
    }

}
