//
//  Script to track propagating front across vascular network
//  Blood haematocrit is assigned accordingly
//
//  Created by Paul Sweeney on 18/02/2020.
//

#include <discrete_flow.h>
#include "global_variables.h"
#include "global_params.h"
#include "output_data.h"

// INPUT
// vvar - variable to track across vascular network - e.g. blood pressure, concentration etc. (size: nnod)
// blabel - branching vessel labels (size: nseg)

// memeffects = 0 -> empirical law (Pries & Secomb 2015)
// memeffects = 1 -> semi-empirical law (Bernabeu et al. (bioRxiv, 2019))

void vhaem(const int &memeffects, ivec &zflow, const ivec &blabel)    {

    diam *= 1e3;

    if (accu(nodpress) == 0.0)  {
        //printf("*** Warning: nodal pressures have not been computed ***");
    }

    // Preamble
    // If one segment in branching vessel has zero flow -> all do (can be slight variance)
    for (int iseg = 0; iseg < nseg; iseg++) {
        if (q(iseg) == 0.0) {
            zflow(find(blabel == blabel(iseg))).fill(1);
        }
    }
    // Compute updated nodtyp list which excludes stagnant branching vessels and store labels
    int inod1, inod2;



    // Store start/end nodes of branching vessels - based on nodal pressure
    // Col. index: 0: high pressure, 1: low pressure
    int rowcnt = 0;
    double press1,press2;
    uvec findseg, vertice;
    imat edgnod = zeros<imat>(nseg,2);
    edgnod.fill(-1);
    for (int iseg = 0; iseg < nseg; iseg++) {
        if (edgnod(iseg,0) == -1 || edgnod(iseg,1) == -1)   {
            rowcnt = 0;
            vertice = zeros<uvec>(0);
            findseg = find(blabel == blabel(iseg)); // Finds edge segments
            // Check if edge only consists of one segment
            // Order of start/end nodes for stagnant flow not important
            if (findseg.n_elem == 1) {
                if (iseg == 8251)   {
                    //cout<<"oh no!"<<endl;
                }
                if (nodpress(ista(iseg)) > nodpress(iend(iseg)))    {
                    edgnod(iseg,0) = ista(iseg);
                    edgnod(iseg,1) = iend(iseg);
                }
                else if (nodpress(iend(iseg)) > nodpress(ista(iseg)))   {
                    edgnod(iseg,1) = ista(iseg);
                    edgnod(iseg,0) = iend(iseg);
                }
            }
            // > 1 segment
            else {
                for (int jseg = 0; jseg < findseg.n_elem; jseg++)   { // Cycle through edge segments
                    if (nodtyp(ista(findseg(jseg))) == 1 || nodtyp(ista(findseg(jseg))) > 2)    { // Store boundary or bifurcating nodes
                        vertice.insert_rows(rowcnt, 1);
                        vertice(rowcnt) = ista(findseg(jseg));
                        rowcnt += 1;
                    }
                    else if (nodtyp(iend(findseg(jseg))) == 1 || nodtyp(iend(findseg(jseg))) > 2)    {
                        vertice.insert_rows(rowcnt, 1);
                        vertice(rowcnt) = iend(findseg(jseg));
                        rowcnt += 1;
                    }
                }
                if (iseg == 8251)   {
                    //cout<<"check"<<endl;
                    //cout<<size(vertice)<<endl;
                }
                if (vertice.n_rows == 1)    {
                    //printf("\t*** Warning: looping branching vessel detected ***\n");
                    ////cout<<blabel(findseg(0))<<"\t"<<findseg.n_elem<<endl;
                    for (int i = 0; i < findseg.n_elem; i++)    {
                        ////cout<<"s: "<<ista(findseg(i))<<"\t"<<nodtyp(ista(findseg(i)))<<"\t"<<snodtyp(ista(findseg(i)))<<endl;
                        ////cout<<"e: "<<iend(findseg(i))<<"\t"<<nodtyp(iend(findseg(i)))<<"\t"<<snodtyp(iend(findseg(i)))<<endl;
                        //edgnod(findseg(i),0) = vertice(0);
                        //edgnod(findseg(i),1) = vertice(0);
                    }
                    for (int i = 0; i < nseg; i++)  {
                        if (ista(i) == iend(findseg(findseg.n_elem-1)) || iend(i) == iend(findseg(findseg.n_elem-1)))   {
                            ////cout<<"ss: "<<ista(i)<<"\t"<<nodtyp(ista(i))<<"\t"<<snodtyp(ista(i))<<"\t"<<blabel(i)<<endl;
                            ////cout<<"ee: "<<iend(i)<<"\t"<<nodtyp(iend(i))<<"\t"<<snodtyp(iend(i))<<"\t"<<blabel(i)<<endl;
                        }
                    }
                    ////cout<<"vertice:"<<endl;
                    ////cout<<vertice(0)<<"\t"<<nodtyp(vertice(0))<<"\t"<<snodtyp(vertice(0))<<endl;
                    //snodtyp(vertice(0)) -= 2;
                    zflow(findseg).fill(1);

                }
                else {
                    press1 = nodpress(vertice(0));
                    press2 = nodpress(vertice(1));
                    if (iseg == 8251)   {
                        //cout<<"even stranger!"<<endl;
                        //printf("%.24lf\n",press1);
                        //printf("%.24lf\n",press1);
                    }
                    if (press1 >= press2)    {
                        for (int jseg = 0; jseg < findseg.n_elem; jseg++)   {
                            edgnod(findseg(jseg), 0) = vertice(0);
                            edgnod(findseg(jseg), 1) = vertice(1);
                        }
                        //edgnod.submat(findseg,zeros<uvec>(1)).fill(vertice(0));
                        //edgnod.submat(findseg,ones<uvec>(1)).fill(vertice(1));
                    }
                    else if (press1 < press2)   {
                        for (int jseg = 0; jseg < findseg.n_elem; jseg++)   {
                            edgnod(findseg(jseg), 1) = vertice(0);
                            edgnod(findseg(jseg), 0) = vertice(1);
                        }
                        //edgnod.submat(findseg,ones<uvec>(1)).fill(vertice(0));
                        //edgnod.submat(findseg,zeros<uvec>(1)).fill(vertice(1));
                    }
                }
            }
        }
    }
    // Need to update for update snodtyp!!
    ivec snodtyp = zeros<ivec>(nnod);
    imat snodnod = zeros<imat>(nodsegm,nnod);
    imat snodseg = zeros<imat>(nodsegm,nnod);
    for (int iseg = 0; iseg < nseg; iseg++) {
        if (zflow(iseg) == 0 && nodpress(edgnod(iseg, 0)) != nodpress(edgnod(iseg, 1)))   {
            inod1 = (int) ista(iseg);
            inod2 = (int) iend(iseg);
            snodtyp(inod1) += 1;
            snodtyp(inod2) += 1;
            if (snodtyp(inod1) > nodtyp(inod1)) {
                //printf("*** Error: Too many segments connected to node %i\n", inod1);
            }
            if (snodtyp(inod2) > nodtyp(inod2)) {
                //printf("*** Error: Too many segments connected to node %i\n", inod2);
            }
            snodseg(snodtyp(inod1) - 1,inod1) = iseg;
            snodseg(snodtyp(inod2) - 1,inod2) = iseg;
            snodnod(snodtyp(inod1) - 1,inod1) = inod2;
            snodnod(snodtyp(inod2) - 1,inod2) = inod1;
        }
        if (zflow(iseg) == 0 && nodpress(edgnod(iseg, 0)) == nodpress(edgnod(iseg, 1)))   {
            zflow(iseg) = 1;
        }
    }

    // Focusing on snodtyp =  3, calculating whether the node is diverging or converging
    ivec bifurtyp = zeros<ivec>(nnod); // 1: diverging, 2: converging bifurcation
    /*for (int inod = 0; inod < nnod; inod++) {
        if (snodtyp(inod) == 3) { // If condition excluded accu(bifurtyp) = nseg
            for (int j = 0; j < snodtyp(inod); j++) {
                if (nodpress(snodseg(j, inod)) > nodpress(inod)) {
                    bifurtyp(inod) += 1;
                }
            }
        }
    }*/
    for (int inod = 0; inod < nnod; inod++) {
        if (snodtyp(inod) == 3) {
            for (int j = 0; j < 3; j++) {
                int nod1 = edgnod(snodseg(j,inod),0);
                int nod2 = edgnod(snodseg(j,inod), 1);
                if (nod1 == inod)  {
                    if (nodpress(nod2) > nodpress(inod))    {
                        bifurtyp(inod) += 1;
                    }
                }
                else if (nod2 == inod)  {
                    if (nodpress(nod1) > nodpress(inod))    {
                        bifurtyp(inod) += 1;
                    }
                }
            }
        }
    }


    // Main algorithm
    // Find inflowing and zero flow boundary branching vessels
    uvec tempvec, segs=zeros<uvec>(3);
    vec flow=zeros<vec>(3);
    ivec vflag = zeros<ivec>(nseg); // Track vessels
    ivec vfront = zeros<ivec>(nnod); // Keeps track of propagating nodal front
    for (int inodbc = 0; inodbc < nnodbc; inodbc++) {
        int nod = bcnod(inodbc);
        if (nodpress(nodnod(nodtyp(nod)-1,nod)) < nodpress(nod))   {
            tempvec = find(blabel == blabel(nodseg(nodtyp(nod)-1,nod))); // find corresponding branching vessels
            vflag(tempvec).fill(1);
            hd(tempvec).fill(bchd(inodbc));
            vfront(edgnod(nodseg(nodtyp(nod)-1,nod),1)) = 1; // Track propagating front
        }
        else if (nodpress(nodnod(nodtyp(nod)-1,nod)) == nodpress(nod))  { // Track zero flow branching vessels at boundary
            tempvec = find(blabel == blabel(nodseg(nodtyp(nod)-1,nod)));
            vflag(tempvec).fill(-1); // zero flow segment
            hd(tempvec).fill(0.0); // Assuming zero Hd
        }
    }

    // Fill zero flow vessels with zero Hd
    hd(find(zflow == 1)).fill(0.0);
    int esc_fac = 0;
    while (accu(vflag(find(vflag >= 0))) < nseg - accu(zflow) || esc_fac <= 5)  {

        for (int inod = 0; inod < nnod; inod++) {
            if (vfront(inod) == 1)  {
                //cout<<"test0"<<endl;
                // Check node type
                if ((snodtyp(inod) == 1 || snodtyp(inod) == 0)  && nodtyp(inod) == 1)   { // Outflowing boundary
                    vfront(inod) = 0;
                }
                else {
                    // Find all connected edges (branching vessels) exc. stagnant
                    int pcntr = 0;
                    int dcntr = 0;
                    uvec ndaught = zeros<uvec>(0);
                    uvec nparent = zeros<uvec>(0);
                    //cout<<"test1"<<endl;
                    for (int nodt = 0; nodt < snodtyp(inod); nodt++) { // Cycle through segments associated to node
                        //cout<<zflow(nodseg(nodt,inod))<<"\t"<<edgnod(nodseg(nodt, inod), 1)<<"\t"<<edgnod(nodseg(nodt, inod), 0)<<"\t"<<inod<<endl;
                        //cout<<"\t"<<nodtyp(inod)<<"\t"<<snodtyp(inod)<<endl;
                        if (edgnod(nodseg(nodt,inod),0) != -1)  {
                            //cout<<"\t"<<nodtyp(edgnod(nodseg(nodt, inod), 0))<<"\t"<<snodtyp(edgnod(nodseg(nodt, inod), 0))<<endl;
                        }
                        if (zflow(snodseg(nodt, inod)) == 0)    {
                            // Log parent and daughter branching segments
                            if (edgnod(snodseg(nodt, inod), 1) == inod)    { // Parent found
                                nparent.insert_rows(pcntr, 1);
                                nparent(pcntr) = snodseg(nodt, inod);
                                if (nodpress(edgnod(snodseg(nodt, inod), 0)) < nodpress(edgnod(snodseg(nodt, inod), 1)))  {
                                    //printf("*** Error: daughter vessel incorrectly assigned as parent ***\n");
                                    //printf("\t\t Edge start/end pressures (mmHg): %.3f / %.3f\n",nodpress(edgnod(nodseg(nodt, inod), 0)),nodpress(edgnod(nodseg(nodt, inod), 1)));
                                }
                                pcntr += 1;
                            }
                            else if (edgnod(snodseg(nodt, inod), 0) == inod) { // Daughter found
                                ndaught.insert_rows(dcntr, 1);
                                ndaught(dcntr) = snodseg(nodt, inod); // Store segment indices
                                if (nodpress(edgnod(snodseg(nodt, inod), 0)) < nodpress(edgnod(snodseg(nodt, inod), 1)))  {
                                    //printf("*** Error: parent vessel incorrectly assigned as daughter ***\n");
                                    //printf("\t\t Edge start/end pressures (mmHg): %.3f / %.3f\n",nodpress(edgnod(nodseg(nodt, inod), 0)),nodpress(edgnod(nodseg(nodt, inod), 1)));
                                }
                                dcntr += 1;
                            }
                       }
                    }

                    // Propagate downstream if only one connected edge (exc. stagnant)
                    if (snodtyp(inod) == 2) {
                        //cout<<"test2"<<endl;
                        //cout<<size(ndaught)<<"\t"<<size(nparent)<<endl;
                        //cout<<"\t\t"<<edgnod(nparent(0),0)<<"\t"<<edgnod(nparent(0),1)<<"\t"<<nodtyp(edgnod(nparent(0),1))<<"\t"<<snodtyp(edgnod(nparent(0),1))<<endl;
                        if (ndaught.n_elem > 0) { // THIS IS KINDA CHEATING! - excluding zero flow can create loops
                            //cout<<"\t\t"<<edgnod(ndaught(0),0)<<"\t"<<edgnod(ndaught(0),1)<<"\t"<<nodtyp(edgnod(ndaught(0),0))<<"\t"<<snodtyp(edgnod(ndaught(0),0))<<endl;
                            hd(find(blabel == blabel(ndaught(0)))).fill(hd(nparent(0))); // RBCs propagate along downstream edge
                        //cout<<"test2.1"<<endl;
                        //cout<<ndaught(0)<<endl;
                        //cout<<find(blabel == blabel(ndaught(0)))<<endl;
                            vflag(find(blabel == blabel(ndaught(0)))).fill(1); // Flag new edge
                            vflag(find(blabel == blabel(nparent(0)))).fill(1);
                        //cout<<"test2.2"<<endl;
                            vfront(edgnod(ndaught(0), 1)) = 1; // Propagate new front
                        //cout<<"test2.3"<<endl;
                        }
                        vfront(inod) = 0; // Switch off old front
                    }
                    else {
                        //cout<<"test3"<<endl;
                        if (snodtyp(inod) == 3) { //
                            if (bifurtyp(inod) == 1)    { // Diverging - apply phase separation law
                                //cout<<"test4"<<endl;
                                //cout<<size(ndaught)<<"\t"<<snodtyp(inod)<<"\t"<<nodtyp(inod)<<endl;
                                //cout<<size(nparent)<<endl;
                                //printf("%.16lf\n",nodpress(inod));
                                for (int j = 0; j < snodtyp(inod); j++) {
                                    //printf("\t %.16lf\n",nodpress(nodnod(j,inod)));
                                    //printf("\t\t %.16lf\n",segpress(nodseg(j,inod)));
                                    //printf("\t\t\t %.16lf\n",q(nodseg(j,inod)));
                                    //printf("\t\t\t\t %lli\n",zflow(nodseg(j,inod)));
                                    //cout<<edgnod(snodseg(j, inod), 0)<<"\t"<<edgnod(snodseg(j, inod), 1)<<endl;
                                    //cout<<snodseg(j, inod)<<endl;
                                }
                                segs(0) = ndaught(0);
                                segs(1) = ndaught(1);
                                segs(2) = nparent(0);
                                if (vflag(nparent(0)) == 0)   {
                                    cout<<"noo!"<<endl;
                                }
                                //cout<<"test4.1"<<endl;
                                flow = abs(q(segs));
                                if (memeffects == 0)    {
                                    womem(segs, flow, diam, hd);
                                }
                                else {
                                    wmem(segs, flow, v_lseg, diam, hd);
                                }
                                //cout<<"test4.2"<<endl;
                                if (hd(segs(0)) > 1.0 || hd(segs(1)) > 1.0)   {
                                    ////cout<<hd(segs(0))<<"\t"<<hd(segs(1))<<"\t"<<hd(segs(2))<<endl;
                                }
                                for (int i = 0; i < ndaught.n_elem; i++)   {
                                    hd(find(blabel == blabel(ndaught(i)))).fill(hd(ndaught(i)));
                                    vflag(find(blabel == blabel(ndaught(i)))).fill(1);
                                    vfront(edgnod(ndaught(i), 1)) = 1;
                                }
                                vflag(find(blabel == blabel(nparent(0)))).fill(1);
                                vfront(inod) = 0;
                                //cout<<"test4.3"<<endl;
                            }
                            else if (bifurtyp(inod) == 2)   { // Converging - use conservation
                                //cout<<"test5"<<endl;
                                //cout<<size(ndaught)<<"\t"<<snodtyp(inod)<<"\t"<<nodtyp(inod)<<endl;
                                segs(0) = ndaught(0);
                                segs(1) = nparent(0);
                                segs(2) = nparent(1);
                                flow = abs(q(segs));
                                if (accu(vflag(nparent)) == nparent.n_elem) {
                                    // Apply conservation
                                    hd(segs(0)) = (flow(1)*hd(segs(1)) + flow(2)*hd(segs(2))) / flow(0);
                                    hd(find(blabel == blabel(ndaught(0)))).fill(hd(ndaught(0)));
                                    vflag(find(blabel == blabel(ndaught(0)))).fill(1);
                                    vfront(edgnod(ndaught(0), 1)) = 1;
                                    vfront(inod) = 0;
                                    for (int i = 0; i < nparent.n_elem; i++)    {
                                        vflag(find(blabel == blabel(nparent(i)))).fill(1);
                                    }
                                } // Otherwise wait
                            }
                        }
                        else {  // No phase separation - split Hd based on flow
                            if (accu(vflag(nparent)) == nparent.n_elem) { // Check parent vessels have been processed
                                //cout<<"test6"<<endl;
                                double flowsum = 0.;
                                double hq = 0.;
                                for (int i = 0; i < nparent.n_elem; i++)   {
                                    flowsum += abs(q(nparent(i)));
                                    hq += hd(nparent(i))*abs(q(nparent(i)));
                                }
                                for (int i = 0; i < ndaught.n_elem; i++)  {
                                    hd(ndaught(i)) = hq/flowsum;
                                    hd(find(blabel == blabel(ndaught(i)))).fill(hd(ndaught(i)));
                                }
                                for (int jnod = 0; jnod < ndaught.n_elem; jnod++)   {
                                    vflag(find(blabel == blabel(ndaught(jnod)))).fill(1);
                                    vfront(edgnod(ndaught(jnod),1)) = 1;
                                }
                                vfront(inod) = 0;
                            }
                        }
                    }
                }
            }
        }


        cout<<accu(vflag(find(vflag >= 0)))<<"\t"<<nseg - accu(zflow)<<endl;
        cout<<esc_fac<<endl;
        //cout<<"\t\t"<<nseg<<"\t"<<accu(zflow)<<endl;
        if (accu(vflag(find(vflag >= 0))) != nseg - accu(zflow)) {
            esc_fac += 1;
        }
    }
    cout<<accu(vflag(find(vflag >= 0)))<<"\t"<<nseg - accu(zflow)<<endl;

    diam *= 1e-3;
}