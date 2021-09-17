#include "spatGraph.hpp"

using namespace reanimate;

void spatGraph::classifyNetwork(imat &baseVessels, ivec &geometry)  {

    printText("Classifying vasculature based on network topology");
    printNum("Dmin =", Dmin);
    printNum("Dr =", Dr);

    Pa = zeros<ivec>(nseg);
    Rs = zeros<ivec>(nseg);
    feedNod = zeros<ivec>(nseg);
    dFeedNod = zeros<ivec>(nseg);
    drainNod = zeros<ivec>(nseg);
    flagTree = zeros<ivec>(nseg);
    daughter = zeros<ivec>(nseg);
    Rn = zeros<ivec>(nnod);

    int cntr = 0;
    int oldcntr = 0;

    flagTree.fill(-1);

    // Ordering in terms of feeding/draining vessel diameters (descending)
    for(int m = 0; m < (int) baseVessels.n_rows; m++)    {
        for(int j = 0; j < (int) baseVessels.n_rows; j++)    {
            if(baseVessels(m, 2) > baseVessels(j, 2))   {
                int tmp1 = (int) baseVessels(m, 0);
                baseVessels(m, 0) = baseVessels(j, 0);
                baseVessels(j, 0) = tmp1;
                int tmp2 = (int) baseVessels(m, 1);
                baseVessels(m, 1) = baseVessels(j, 1);
                baseVessels(j, 1) = tmp2;
                double tmp3 = baseVessels(m, 2);
                baseVessels(m, 2) = baseVessels(j, 2);
                baseVessels(j, 2) = tmp3;
                double tmp4 = baseVessels(m, 3);
                baseVessels(m, 3) = baseVessels(j, 3);
                baseVessels(j, 3) = tmp4;
                double tmp5 = baseVessels(m, 4);
                baseVessels(m, 4) = baseVessels(j, 4);
                baseVessels(j, 4) = tmp5;
            }
        }
    }


    // Beginning with the largest-diameter starting segment, Algorithm 1 is implemented where every parent segment has a diameter >= Dmin.
    // Daughter segments with D < Dmin are added to Rs but their end nodes aren't labelled.
    // These steps are repeated for all the remaining starting segments.
    for (int iTree = 0; iTree < (int) baseVessels.n_rows; iTree++)  {

        int init_seg = (int) baseVessels(iTree, 0);
        Rs(init_seg) = 1;   // Flag as a reached segment
        Rn(ista(init_seg)) = 1; // Flag parent start node
        Rn(iend(init_seg)) = 1; // Flag parent end node
        flagTree(init_seg) = iTree; // Flag vessel as part of the jth tree

        if (nodtyp(ista(init_seg)) == 1)    {
            feedNod(init_seg) = 1;
        }
        else if (nodtyp(iend(init_seg)) == 1)   {
            feedNod(init_seg) = 2;
        }


        for (int jseg = 0; jseg < nseg; jseg++) {

            daughter.zeros();
            dFeedNod.zeros();

            // Search for max. diameter segment which hasn't been a parent
            int parent = 0;
            float temp = 0.0;
            for (int iseg = 0; iseg < nseg; iseg++)    {
                if ((Rs(iseg) == 1) && (Pa(iseg) == 0) && (diam(iseg) > temp) && (diam(iseg) >= Dmin) && (flagTree(iseg) == iTree))  {
                    temp = diam(iseg);
                    parent = iseg;
                }
            }
            Pa(parent) = 1;

            // Log drain node of parent
            if (feedNod(parent) == 1)   {
                drainNod(parent) = iend(parent);
            }
            else if (feedNod(parent) == 2)  {
                drainNod(parent) = ista(parent);
            }


            // Find daughter(s)
            int dCounter = 0;
            for (int iseg = 0; iseg < nseg; iseg++)    {
                if (parent != iseg)    {
                    if (ista(iseg) == drainNod(parent))   {
                        if (Rs(iseg) == 0)   {
                            feedNod(iseg) = 1;
                            Rs(iseg) = 1;
                        }
                        dFeedNod(iseg) = 1;
                        daughter(iseg) = 1;
                        flagTree(iseg) = iTree;
                        dCounter += 1;
                    }
                    else if (iend(iseg) == drainNod(parent))  {
                        if (Rs(iseg) == 0)  {
                            feedNod(iseg) = 2;
                            Rs(iseg) = 1;
                        }
                        dFeedNod(iseg) = 2;
                        daughter(iseg) = 1;
                        flagTree(iseg) = iTree;
                        dCounter += 1;
                    }
                }
            }


            if (dCounter == 0 && geometry(parent) == 0) {
                geometry(parent) = 1;
            }
            else if (dCounter != 0) {

                uvec temp = find(daughter == 1);
                vec d_order = conv_to< vec >::from(temp);
                vec d_diam = diam(find(daughter == 1));
                mat d_matrix = join_rows(d_order,d_diam);
                for(int m = 0; m < dCounter; m++)    {
                    for(int j = 0; j < dCounter; j++)    {
                        if(d_matrix(m,1) > d_matrix(j,1))   {
                            double t1 = d_matrix(m,0);
                            d_matrix(m,0) = d_matrix(j,0);
                            d_matrix(j,0) = t1;
                            double t2 = d_matrix(m,1);
                            d_matrix(m,1) = d_matrix(j,1);
                            d_matrix(j,1) = t2;
                        }
                    }
                }


                for (int j = 0; j < dCounter; j++) {
                    int bigD = d_matrix(j,0);
                    int inod = 0;
                    if (dFeedNod(bigD) == 1) {
                        inod = (int) iend(bigD);
                    }
                    else if (dFeedNod(bigD) == 2)    {
                        inod = (int) ista(bigD);
                    }
                    if (Rn(inod) == 1)  {
                        if (geometry(parent) != 1 && geometry(parent) != 3) {
                            geometry(bigD) = 2;
                            if (diam(parent) < (Dr*diam(bigD))) {
                                geometry(parent) = 2;
                            }
                        }
                    }
                    else    {
                        if (diam(bigD) > Dmin) {
                            Rn(inod) = 1;
                        }
                    }
                }

            }


            if (geometry(parent) == 2)    {
                for (int iseg = 0; iseg < nseg; iseg++) {
                    if (daughter(iseg) == 1 && geometry(iseg) == 0)    {
                        geometry(iseg) = 2;
                    }
                }
            }
            else if (geometry(parent) == 0) {
                geometry(parent) = 1;
            }

            int check = 0;
            for (int iseg = 0; iseg < nseg; iseg++)    {
                if ((Rs(iseg) == 1) && (Pa(iseg) == 0) && (diam(iseg) >= Dmin) && (flagTree(iseg) == iTree))  {
                    check = 1;
                }
            }
            if (check == 0) {
                goto here;
            }

        }
        here:;

        // Altering vector based on classification
        if (baseVessels(iTree, 1) == 2)  {
            for (int iseg = 0; iseg < nseg; iseg++)    {
                if (geometry(iseg) == 1 && flagTree(iseg) == iTree)    {
                    geometry(iseg) = 3;
                }
            }
        }


    }

    loop:;
    for (int jTree = 0; jTree < (int) baseVessels.n_rows; jTree++)   {
        internalClassificationLoop(0, (int) baseVessels(jTree, 1), geometry, 1, jTree, flagTree);
    }

    // Check to see if there any remaining candidate parents segments
    cntr = 0;
    for (int iseg = 0; iseg < nseg; iseg++) {
        if ((Rs(iseg) == 1) && (Pa(iseg) == 0)) {
            cntr += 1;
        }
    }

    if (oldcntr == cntr)    {
        goto jump;
    }

    for (int iseg = 0; iseg < nseg; iseg++)    {
        if ((Rs(iseg) == 1) && (Pa(iseg) == 0))  {
            goto loop;
        }
    }
    jump:;

}

void spatGraph::internalClassificationLoop(const int &init_seg, const int &classify, ivec &geometry, const int &algo_2, const int &tree, ivec &fTree)    {

    int loop;
    if (algo_2 == 1)    {
        loop = 1;
    }
    else    {

        Pa.zeros();
        Rs.zeros();
        feedNod.zeros();    // Vessel feed node
        drainNod.zeros();   // Vessel drain node
        dFeedNod.zeros();   // Daughter feed node
        daughter.zeros();   // List of daughter vessels
        Rn.zeros();

        loop = nseg;
        Rs(init_seg) = 1;   // Flag as a reached segment
        Rn(ista(init_seg)) = 1; // Flag parent start node
        Rn(iend(init_seg)) = 1; // Flag parent end node

        if (nodtyp(ista(init_seg)) == 1)    {
            feedNod(init_seg) = 1;
        }
        else if (nodtyp(iend(init_seg)) == 1)   {
            feedNod(init_seg) = 2;
        }
        else    {
            printText("Input not a boundary node",4);
            printText("Exciting classifier");
            goto exit;
        }
    }

    for (int jseg = 0; jseg < loop; jseg++) {

        daughter.zeros();
        dFeedNod.zeros();

        // Search for max. diameter segment which hasn't been a parent
        int parent = 0;
        double tmp = 0.0;
        for (int iseg = 0; iseg < nseg; iseg++)    {
            if ((Rs(iseg) == 1) && (Pa(iseg) == 0) && (diam(iseg) > tmp) && (fTree(iseg) == tree))  {
                tmp = diam(iseg);
                parent = iseg;
            }
        }
        Pa(parent) = 1;

        // Log drain node of parent
        if (feedNod(parent) == 1)   {
            drainNod(parent) = iend(parent);
        }
        else if (feedNod(parent) == 2)  {
            drainNod(parent) = ista(parent);
        }


        // Find daughter(s)
        int dCounter = 0;
        for (int iseg = 0; iseg < nseg; iseg++)    {
            if (parent != iseg && Pa(iseg) == 0)    {
                if (ista(iseg) == drainNod(parent))   {
                    if (Rs(iseg) == 0)   {
                        feedNod(iseg) = 1;
                        Rs(iseg) = 1;
                    }
                    dFeedNod(iseg) = 1;
                    daughter(iseg) = 1;
                    fTree(iseg) = tree;
                    dCounter += 1;
                }
                else if (iend(iseg) == drainNod(parent))  {
                    if (Rs(iseg) == 0)  {
                        feedNod(iseg) = 2;
                        Rs(iseg) = 1;
                    }
                    dFeedNod(iseg) = 2;
                    daughter(iseg) = 1;
                    fTree(iseg) = tree;
                    dCounter += 1;
                }
            }
        }


        if (dCounter == 0 && geometry(parent) == 0) {
            if (classify == 1)  {
                geometry(parent) = 1;
            }
            else if (classify == 2) {
                geometry(parent) = 3;
            }
        }
        else if (dCounter != 0) {

            uvec tmp2 = find(daughter == 1);
            vec d_order = conv_to< vec >::from(tmp2);
            vec d_diam = diam(find(daughter == 1));
            mat d_matrix = join_rows(d_order,d_diam);
            for(int m = 0; m < dCounter; m++)    {
                for(int j = 0; j < dCounter; j++)    {
                    if(d_matrix(m,1) > d_matrix(j,1))   {
                        double t1 = d_matrix(m,0);
                        d_matrix(m,0) = d_matrix(j,0);
                        d_matrix(j,0) = t1;
                        double t2 = d_matrix(m,1);
                        d_matrix(m,1) = d_matrix(j,1);
                        d_matrix(j,1) = t2;
                    }
                }
            }


            for (int j = 0; j < dCounter; j++) {
                int bigD = d_matrix(j,0);
                int inod = 0;
                if (dFeedNod(bigD) == 1) {
                    inod = (int) iend(bigD);
                }
                else if (dFeedNod(bigD) == 2)    {
                    inod = (int) ista(bigD);
                }
                if (Rn(inod) == 1)  {
                    if (geometry(parent) != 2) {
                        geometry(bigD) = 2;
                        if (diam(parent) < (Dr*diam(bigD))) {
                            geometry(parent) = 2;
                            /*int inod2 = 0;
                            int findEdge = 0;
                            if (feedNod(parent) == 1)   {
                                inod2 = (int) ista(parent);
                            }
                            else if (feedNod(parent) == 2)  {
                                inod2 = (int) iend(parent);
                            }
                            while (findEdge == 0)   {
                                for (int iseg = 0; iseg < nseg; iseg++) {
                                    if (inod2 == (int) ista(iseg))  {
                                        if (nodtyp(inod2) < 3) {
                                            int cntr = 0;
                                            for (int ives = 0; ives < NoOfTrees; ives++)    {
                                                if (InOutlets(ives,3) == inod2) {
                                                    cntr += 1;
                                                }
                                            }
                                            if (cntr == 0)  {
                                                vesselGeometry(iseg) = 2;
                                            }
                                        }
                                        inod2 = (int) iend(iseg);
                                        if (nodtyp(inod2) > 2)  {
                                            findEdge = 1;
                                        }
                                    }
                                    else if (inod2 == (int) iend(iseg))  {
                                        if (nodtyp(inod2) < 3) {
                                            int cntr = 0;
                                            for (int ives = 0; ives < NoOfTrees; ives++)    {
                                                if (InOutlets(ives,3) == inod2) {
                                                    cntr += 1;
                                                }
                                            }
                                            if (cntr == 0)  {
                                                vesselGeometry(iseg) = 2;
                                            }
                                        }
                                        inod2 = (int) ista(iseg);
                                        if (nodtyp(inod2) > 2)  {
                                            findEdge = 1;
                                        }
                                    }
                                }
                            }
                             */
                        }
                    }
                }
                else    {
                    Rn(inod) = 1;
                }
            }

        }


        if (geometry(parent) == 2)    {
            for (int iseg = 0; iseg < nseg; iseg++) {
                if (daughter(iseg) == 1 && geometry(iseg) == 0)    {
                    geometry(iseg) = 2;
                }
            }
        }
        else if (geometry(parent) == 0) {
            if (classify == 1)  {
                geometry(parent) = 1;
            }
            else if (classify == 2) {
                geometry(parent) = 3;
            }
        }


        if (algo_2 == 1)    {
            if (Rn(drainNod(parent)) == 0)  {
                Rn(drainNod(parent)) = 1;
            }
        }
    }



    exit:;
}


void spatGraph::mapClassification(Network &net, bool fullgraph) {

    // Vessel classification
    net.flagTree = zeros<ivec>(net.getNseg());
    if (fullgraph)  {
        uvec idx;
        for (int iseg = 0; iseg < nseg; iseg++) {
            idx = find(segname(iseg) == net.edgeLabels);
            net.vesstyp(idx).fill(vesselGeometry(iseg));
            net.flagTree(idx).fill(flagTree(iseg));
        }
    }
    else {
        net.vesstyp.ones() * 2;
        flagTree.zeros();
        ivec tag = -ones<ivec>(net.getNnod());
        for (int i = 0; i < nTrees; i++)    {
            int nod = InOutlets(i,4);
            uvec idx = find(net.nodname == nod);
            net.dfsBasic(idx(0), 1, tag);
            idx = find(tag == 1);
            for (int iseg = 0; iseg < net.getNseg(); iseg++)    {
                for (int inod = 0; inod < (int) idx.n_elem; inod++) {
                    if (net.ista(iseg) == idx(inod) || net.iend(iseg) == idx(inod)) {
                        net.flagTree(iseg) = nod;
                        if (InOutlets(i,1) == 1)    {net.vesstyp(iseg) = 1;}
                        else {net.vesstyp(iseg) = 3;}
                    }
                }
            }
            tag.fill(-1);
        }
    }

    // Boundary vesselGeometry
    net.analyseBoundaryType();

}