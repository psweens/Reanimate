#include "Network.hpp"

using namespace reanimate;

void Network::putrank()  {

    nodtyp.zeros();
    nodout.zeros();
    // Rearranging 'nodseg' and 'nodnod' in terms of flow -> output nodes precede input nodes
    // Inflowing nodes first -> 'nodout' indicates the no. of outflowing segments connected to a given node
    int nod1, nod2;
    for (int iseg = 0; iseg < nseg; iseg++) {
        if (q(iseg) >= 0.){
            nod1 = (int) ista(iseg);
            nod2 = (int) iend(iseg);
        }
        else {
            nod1 = (int) iend(iseg);
            nod2 = (int) ista(iseg);
        }
        nodtyp(nod1) += 1;
        nodseg(nodtyp(nod1)-1,nod1) = iseg;
        nodnod(nodtyp(nod1)-1,nod1) = nod2;
        nodout(nod1) += 1;
    }

    // Store outflowing nodes second
    for (int iseg = 0; iseg < nseg; iseg++) {
        if (q(iseg) >= 0.){
            nod1 = (int) ista(iseg);
            nod2 = (int) iend(iseg);
        }
        else {
            nod1 = (int) iend(iseg);
            nod2 = (int) ista(iseg);
        }
        nodtyp(nod2) += 1;
        nodseg(nodtyp(nod2)-1,nod2) = iseg;
        nodnod(nodtyp(nod2)-1,nod2) = nod1;
    }

    // Assign low ranks to inflow boundary nodes
    nnodfl = 0; // Counts inflowing nodes
    nk.zeros();
    for (int inod = 0; inod < nnod; inod++){
        if (nodtyp(inod) == 1 && nodout(inod) == 1){
            nk(inod) = 1;
            nodrank(nnodfl) = inod;
            nnodfl += 1;
        }
    }

    // Assign increasing ranks to downstream connected nodes (removed 'goto' step due to instability)
    int flag = 1;
    int ignore = 0;
    // 'While' loop terminates when all downstream nodes have been flagged
    int seg;
    while (flag == 1)   {
        flag = 0;
        for (int inod = 0; inod < nnod; inod++){
            // Identify nodes not flagged
            if (nk(inod) == 0 && nodtyp(inod) > 0){
                // Cycle through outflowing segments/nodes
                for (int j = (int) nodout(inod); j < nodtyp(inod); j++){
                    seg = (int) nodseg(j,inod);
                    if (inod == iend(seg) && (nk(ista(seg)) == 0 || q(seg) <= 0.)) {
                        goto skipnode;
                    }
                    else if (inod == ista(seg) && (nk(iend(seg)) == 0 || q(seg) >= 0.)) {
                        goto skipnode;
                    }
                }
                //if (ignore == 0)    {
                nk(inod) = 1;
                nodrank(nnodfl) = inod;
                nnodfl += 1;
                flag = 1;
                skipnode:;
                //}
                //else {
                ignore = 0;
                // }
            }
        }
    }

    /*int cnt_error = 0;
    for(int inod = 0; inod < nnod; inod++)	{
        if(nk(inod) == 0)  {
            printf("*** Error: unprocessed node %i in putrank\n",inod);
            for (int j = (int) nodout(inod); j < nodtyp(inod); j++) {
                seg = (int) nodseg(j,inod);
                cout<<segpress(seg)<<"\t"<<q(seg)<<"\t"<<qq(seg)<<endl;
            }
            cnt_error += 1;
        }
    }*/
    //cout<<cnt_error<<endl;
    //cout<<nodtyp.max()<<"\t"<<nodsegm<<endl;

}