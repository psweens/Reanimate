//
//  resetNet.cpp
//  Network Tools 2.0
//
//  Resets network data if segments / nodes need to be removed
//
//  Created by Paul Sweeney on 25/11/2015.
//  Copyright © 2015 Paul Sweeney - University College London. All rights reserved.
//

#include <armadillo>
#include <iostream>

#include "global_variables.h"
#include "global_params.h"
#include "output_data.h"

using namespace arma;
using namespace std;

void resetNet(ivec &index) {
    
    cout<<"Updating Network..."<<endl;
    for (int iseg = 0; iseg < nseg; iseg++) {
        if (index(iseg) == 1)   {
            segname.shed_row(iseg);
            vesstyp.shed_row(iseg);
            segnodname.shed_col(iseg);
            diam.shed_row(iseg);
            lseg.shed_row(iseg);
            q.shed_row(iseg);
            hd.shed_row(iseg);
            index.shed_row(iseg);
            nseg -= 1;
            iseg = 0;
        }
    }
    
    
    
    ista = zeros<uvec>(nseg);
    iend = zeros<uvec>(nseg);
    for (int iseg = 0; iseg < nseg; iseg++)	{
        //Search for nodes corresponding to this segment
        for (int i = 0; i < 2; i++) {
            for (int inod = 0; inod < nnod; inod++) {
                if (nodname(inod) == segnodname(i,iseg))    {
                    if (i == 0) {
                        ista(iseg) = inod;
                        goto foundit;
                    }
                    else if (i == 1)    {
                        iend(iseg) = inod;
                        goto foundit;
                    }
                }
            }
            printf("*** Error: No matching node found for segname %lli\n", segname(iseg));
        foundit:;
        }
    }
    
    
    //Setup nodtyp, nodseg and nodnod
    nodtyp.zeros();
    for (int iseg = 0; iseg < nseg; iseg++) {
        int inod1 = (int) ista(iseg);
        int inod2 = (int) iend(iseg);
        nodtyp(inod1) += 1;
        nodtyp(inod2) += 1;
        if(nodtyp(inod1) > nodsegm) {
            printf("*** Error: Too many segments connected to node %i\n", inod1);
        }
        if(nodtyp(inod2) > nodsegm) {
            printf("*** Error: Too many segments connected to node %i\n", inod2);
        }
    }
    
    for (int inod = 0; inod < nnod; inod++) {
        if (nodtyp(inod) == 0)  {
            nodname.shed_row(inod);
            nodtyp.shed_row(inod);
            cnode.shed_col(inod);
            inod = 0;
            nnod -= 1;
        }
    }
    
    uvec copyBCnodname = bcnodname;
    uvec copyBCtyp = bctyp;
    vec copyBCprfl = bcprfl;
    
    nnodbc = 0;
    for (int inod = 0; inod < nnod; inod++) {
        if (nodtyp(inod) == 1)  {
            nnodbc += 1;
        }
    }
    
    bcnodname = zeros<uvec>(nnodbc);
    bctyp = zeros<uvec>(nnodbc);
    bcprfl = zeros<vec>(nnodbc);
    bchd = zeros<vec>(nnodbc);
    bchd.fill(0.4);
    
    int jnodbc = 0;
    for (int inod = 0; inod < nnod; inod++) {
        if (nodtyp(inod) == 1)  {
            bcnodname(jnodbc) = nodname(inod);
            for (int knodbc = 0; knodbc < copyBCnodname.n_elem; knodbc++)   {
                if (bcnodname(jnodbc) == copyBCnodname(knodbc)) {
                    bctyp(jnodbc) = copyBCtyp(knodbc);
                    bcprfl(jnodbc) = copyBCprfl(knodbc);
                    goto here;
                }
                else {
                    bctyp(jnodbc) = 3;
                }
            }
        here:;
            jnodbc += 1;
        }
    }
    
    
    
    printf("\t No. of segments = %i\n",nseg);
    printf("\t No. of nodes = %i\n",nnod);
    printf("\t No. of boundary nodes = %i\n",nnodbc);
    
    outputNet("ZeroFlowOutput.txt", nseg, nnod, nnodbc, segname, vesstyp, segnodname, diam, q, hd, nodname, cnode, bcnodname, bctyp, bcprfl, bchd);
}