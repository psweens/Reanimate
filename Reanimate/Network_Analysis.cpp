#include "Network.hpp"

using namespace reanimate;

void Network::edgeNetwork() {

    printText( "Analysing network edges");

    int order = 1;
    int edgIdx = 1;//0;
    ivec flag = zeros<ivec>(nseg);
    edgeLabels = zeros<ivec>(nseg);
    elseg = zeros<vec>(nseg);
    ediam = zeros<vec>(nseg);

    uvec bnod = find(nodtyp == 1 || nodtyp > 2); // Find index of boundary or bifurcating nodes
    for (int inod = 0; inod < bnod.n_elem; inod++)  {   // Cycle through nodtyp != 2

        for (int jnod = 0; jnod < nodtyp(bnod(inod)); jnod++) { // Cycle through segments at each node

            int nod1 = bnod(inod); // Start node
            int nod2 = nodnod(jnod,nod1); // End Node
            int seg = nodseg(jnod,nod1); // Current segment
            if (flag(seg) == 0)  { // If segment has not been processed -> proceed

                // Initialise first segment
                flag(seg) = 1;
                edgeLabels(seg) = edgIdx;
                nod1 = nod2; // End node becomes start node of next segment

                while (nodtyp(nod1) == 2)   {

                    uvec next = find(flag == 0 && (ista == nod1 || iend == nod1)); // Find next node
                    if (next.n_elem > 1)    {printText( "More than one edge segment detected",4);}

                    // Tag & update indices
                    seg = next(0);
                    if (nod1 == ista(seg))  {nod2 = iend(seg);}
                    else {nod2 = ista(seg);}
                    flag(seg) = 1;
                    edgeLabels(seg) = edgIdx;
                    nod1 = nod2;

                }
                edgIdx += 1;
            }

        }

    }
    for (int i = 0; i < edgeLabels.n_elem; i++) {
        if (edgeLabels(i) == 0) {
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
