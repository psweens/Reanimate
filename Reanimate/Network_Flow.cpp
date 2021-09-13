#include "Network.hpp"

using namespace reanimate;

void Network::putrank(Network &sGraph)  {

    //printText("Calculating flow ranking",2, 0);

    sGraph.q.zeros();
    uvec idx, node1, node2;
    for (int iseg = 0; iseg < sGraph.getNseg(); iseg++) {
        idx = find(sGraph.segname(iseg) == edgeLabels);
        node1 = find(sGraph.segnodname(0, iseg) == nodname);
        node2 = find(sGraph.segnodname(1, iseg) == nodname);

        if (nodpress(node1(0)) > nodpress(node2(0)))  {sGraph.q(iseg) = abs(q(idx(0)));}
        else if (nodpress(node1(0)) < nodpress(node2(0)))  {sGraph.q(iseg) = -abs(q(idx(0)));}
        //else {sGraph.q(iseg) = 0.0;}
    }

    sGraph.nodtyp.zeros();
    sGraph.nodout.zeros();
    // Rearranging 'nodseg' and 'nodnod' in terms of flow -> output nodes precede input nodes
    // Inflowing nodes first -> 'nodout' indicates the no. of outflowing segments connected to a given node
    int nod1{}, nod2{};
    for (int iseg = 0; iseg < sGraph.getNseg(); iseg++) {
        if (sGraph.q(iseg) >= 0.){
            nod1 = (int) sGraph.ista(iseg);
            nod2 = (int) sGraph.iend(iseg);
        }
        else {
            nod1 = (int) sGraph.iend(iseg);
            nod2 = (int) sGraph.ista(iseg);
        }
        sGraph.nodtyp(nod1) += 1;
        sGraph.nodseg(sGraph.nodtyp(nod1) - 1, nod1) = iseg;
        sGraph.nodnod(sGraph.nodtyp(nod1) - 1, nod1) = nod2;
        sGraph.nodout(nod1) += 1;
    }

    // Store outflowing nodes second
    for (int iseg = 0; iseg < sGraph.getNseg(); iseg++) {
        if (sGraph.q(iseg) >= 0.){
            nod1 = (int) sGraph.ista(iseg);
            nod2 = (int) sGraph.iend(iseg);
        }
        else {
            nod1 = (int) sGraph.iend(iseg);
            nod2 = (int) sGraph.ista(iseg);
        }
        sGraph.nodtyp(nod2) += 1;
        sGraph.nodseg(sGraph.nodtyp(nod2) - 1, nod2) = iseg;
        sGraph.nodnod(sGraph.nodtyp(nod2) - 1, nod2) = nod1;
    }


    // Assign low ranks to inflow boundary nodes
    sGraph.nnodfl = 0; // Counts inflowing nodes
    sGraph.nk.zeros();
    for (int inod = 0; inod < sGraph.getNnod(); inod++){
        if (sGraph.nodtyp(inod) == 1 && sGraph.nodout(inod) == 1){
            sGraph.nk(inod) = 1;
            sGraph.nodrank(sGraph.nnodfl) = inod;
            sGraph.nnodfl += 1;
        }
    }

    // Assign increasing ranks to downstream connected nodes (removed 'goto' step due to instability)
    // 'While' loop terminates when all downstream nodes have been flagged
    int flag = 1;
    int seg;
    while (flag == 1)   {
        flag = 0;
        for (int inod = 0; inod < sGraph.getNnod(); inod++){
            // Identify nodes not flagged
            if (sGraph.nk(inod) == 0 && sGraph.nodtyp(inod) > 0){
                // Cycle through outflowing segments/nodes
                for (int j = (int) sGraph.nodout(inod); j < sGraph.nodtyp(inod); j++){
                    seg = (int) sGraph.nodseg(j, inod);
                    if (inod == sGraph.iend(seg) && (sGraph.nk(sGraph.ista(seg)) == 0 || sGraph.q(seg) < 0.)) {goto skipnode;}
                    else if (inod == sGraph.ista(seg) && (sGraph.nk(sGraph.iend(seg)) == 0 || sGraph.q(seg) >= 0.)) {goto skipnode;}
                }
                sGraph.nk(inod) = 1;
                sGraph.nodrank(sGraph.nnodfl) = inod;
                sGraph.nnodfl += 1;
                flag = 1;
                skipnode:;
            }
        }
    }

    int cnt_error = 0;
    for(int inod = 0; inod < sGraph.getNnod(); inod++)	{
        if(sGraph.nk(inod) == 0)  {cnt_error += 1;}
    }
    if (cnt_error > 0)  {printText(to_string(cnt_error)+" Unprocessed nodes in putrank",4);}

}


void Network::computeBoundaryFlow()  {

    int nod{};
    BCpress = zeros<vec>(nnodbc);
    BCflow = zeros<vec>(nnodbc);
    for (int inodbc = 0; inodbc < nnodbc; inodbc++) {
        nod = bcnod(inodbc);
        BCpress(inodbc) = nodpress(nod);
        BCflow(inodbc) = abs(q(nodseg(0,nod)));
        if (nodpress(nod) < nodpress(nodnod(0,nod)))    {BCflow(inodbc) *= -1;}
    }

}