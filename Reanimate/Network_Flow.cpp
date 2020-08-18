#include "Network.hpp"

using namespace reanimate;

void Network::putrank(Network &graph)  {

    printText("Calculating flow ranking",2, 0);

    int cntr{};
    flag = zeros<ivec>(nseg);
    for (int iseg = 0; iseg < graph.getNseg(); iseg++) {
        uvec idx = find(graph.segname(iseg) == edgeLabels);
        uvec nod1 = find(graph.segnodname(0,iseg) == nodname);
        uvec nod2 = find(graph.segnodname(1,iseg) == nodname);

        if (nodpress(nod1(0)) > nodpress(nod2(0)))  {graph.q(iseg) = abs(q(idx(0)));}
        else  {graph.q(iseg) = -abs(q(idx(0)));}
    }

    graph.nodtyp.zeros();
    graph.nodout.zeros();
    // Rearranging 'nodseg' and 'nodnod' in terms of flow -> output nodes precede input nodes
    // Inflowing nodes first -> 'nodout' indicates the no. of outflowing segments connected to a given node
    int nod1{}, nod2{};
    for (int iseg = 0; iseg < graph.getNseg(); iseg++) {
        if (graph.q(iseg) >= 0.){
            nod1 = (int) graph.ista(iseg);
            nod2 = (int) graph.iend(iseg);
        }
        else {
            nod1 = (int) graph.iend(iseg);
            nod2 = (int) graph.ista(iseg);
        }
        graph.nodtyp(nod1) += 1;
        graph.nodseg(graph.nodtyp(nod1)-1,nod1) = iseg;
        graph.nodnod(graph.nodtyp(nod1)-1,nod1) = nod2;
        graph.nodout(nod1) += 1;
    }

    // Store outflowing nodes second
    for (int iseg = 0; iseg < graph.getNseg(); iseg++) {
        if (graph.q(iseg) >= 0.){
            nod1 = (int) graph.ista(iseg);
            nod2 = (int) graph.iend(iseg);
        }
        else {
            nod1 = (int) graph.iend(iseg);
            nod2 = (int) graph.ista(iseg);
        }
        graph.nodtyp(nod2) += 1;
        graph.nodseg(graph.nodtyp(nod2)-1,nod2) = iseg;
        graph.nodnod(graph.nodtyp(nod2)-1,nod2) = nod1;
    }


    // Assign low ranks to inflow boundary nodes
    graph.nnodfl = 0; // Counts inflowing nodes
    graph.nk.zeros();
    for (int inod = 0; inod < graph.getNnod(); inod++){
        if (graph.nodtyp(inod) == 1 && graph.nodout(inod) == 1){
            graph.nk(inod) = 1;
            graph.nodrank(graph.nnodfl) = inod;
            graph.nnodfl += 1;
        }
    }

    // Assign increasing ranks to downstream connected nodes (removed 'goto' step due to instability)
    int flag = 1;
    int ignore = 0;
    // 'While' loop terminates when all downstream nodes have been flagged
    int seg;
    while (flag == 1)   {
        flag = 0;
        for (int inod = 0; inod < graph.getNnod(); inod++){
            // Identify nodes not flagged
            if (graph.nk(inod) == 0 && graph.nodtyp(inod) > 0){
                // Cycle through outflowing segments/nodes
                for (int j = (int) graph.nodout(inod); j < graph.nodtyp(inod); j++){
                    seg = (int) graph.nodseg(j,inod);
                    if (inod == graph.iend(seg) && (graph.nk(graph.ista(seg)) == 0 || graph.q(seg) <= 0.)) {
                        goto skipnode;
                    }
                    else if (inod == graph.ista(seg) && (graph.nk(graph.iend(seg)) == 0 || graph.q(seg) >= 0.)) {
                        goto skipnode;
                    }
                }
                //if (ignore == 0)    {
                graph.nk(inod) = 1;
                graph.nodrank(graph.nnodfl) = inod;
                graph.nnodfl += 1;
                flag = 1;
                skipnode:;
                //}
                //else {
                ignore = 0;
                // }
            }
        }
    }

    int cnt_error = 0;
    for(int inod = 0; inod < graph.getNnod(); inod++)	{
        if(graph.nk(inod) == 0)  {
            printText( "Unprocessed node "+to_string(inod)+" in putrank",4);
            cnt_error += 1;
        }
    }

}