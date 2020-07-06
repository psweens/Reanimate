//
//  tort_l.cpp
//  Vascular-Flow
//
//  Created by Paul Sweeney on 20/10/2016.
//  Copyright © 2016 Paul Sweeney - University College London. All rights reserved.
//

#include "global_variables.h"
#include "output_data.h"

ivec tort_l()   {

    int order = 1;
    int edge = 0;
    vec ledge = zeros<vec>(nseg);
    vec v_lseg = zeros<vec>(0);
    vec v_diam = zeros<vec>(0);
    ivec flag = zeros<ivec>(nseg);
    ivec master = zeros<ivec>(nseg);
    master.fill(-1);
    for (int inod = 0; inod < nnod; inod++) {
        if (nodtyp(inod) != 2)   {

            int cntr = 0;
            int branches = (int) nodtyp(inod);
            imat feedNod = zeros<imat>(2,branches);
            //flag.zeros();
            for (int iseg = 0; iseg < nseg; iseg++) {
                if (cntr == branches)   {
                    iseg = nseg;
                }
                else if (flag(iseg) == 0)    {
                    if (ista(iseg) == inod)  {
                        flag(iseg) = 1;
                        feedNod(0,cntr) = iend(iseg);
                        feedNod(1,cntr) = order;
                        master(iseg) = feedNod(1,cntr);
                        order += 1;
                        cntr += 1;
                    }
                    else if (iend(iseg) == inod)    {
                        flag(iseg) = 1;
                        feedNod(0,cntr) = ista(iseg);
                        feedNod(1,cntr) = order;
                        master(iseg) = feedNod(1,cntr);
                        order += 1;
                        cntr += 1;
                    }
                }
                else if (flag(iseg) != 0 && (ista(iseg) == inod || iend(iseg) == inod)) {
                    feedNod.shed_col(cntr);
                    branches -= 1;
                }
            }

            if (branches > 0 && cntr != 0)    {
                for (int i = 0; i < feedNod.n_cols; i++)  {
                    edge = 0;
                    for (int iseg = 0; iseg < nseg; iseg++) {
                        if (ista(iseg) == feedNod(0,i) && flag(iseg) == 0 && nodtyp(feedNod(0,i)) == 2)   {
                            master(iseg) = feedNod(1,i);
                            flag(iseg) = 2;
                            feedNod(0,i) = iend(iseg);
                            iseg = 0;
                        }
                        else if (iend(iseg) == feedNod(0,i) && flag(iseg) == 0 && nodtyp(feedNod(0,i)) == 2)   {
                            master(iseg) = feedNod(1,i);
                            flag(iseg) = 2;
                            feedNod(0,i) = ista(iseg);
                            iseg = 0;
                        }
                    }
                    v_lseg.insert_rows(edge,1);
                    v_lseg(edge) = accu(lseg(find(master == feedNod(1,i))));
                    v_diam.insert_rows(edge,1);
                    v_diam(edge) = mean(diam(find(master == feedNod(1,i))));
                    edge += 1;
                    for (int iseg = 0; iseg < nseg; iseg++) {
                        if (master(iseg) == feedNod(1,i))   {
                            ledge(iseg) = v_lseg(edge-1);
                        }
                    }
                }
            }

        }
    }

    cout<<"Edge lengths: "<<mean(v_lseg)<<" pm "<<stddev(v_lseg)<<endl;
    cout<<"Edge diams: "<<mean(v_diam)<<" pm "<<stddev(v_diam)<<endl;
    cout<<"Length / diam ratio: "<<mean(v_lseg / v_diam)<<" pm "<<stddev(v_lseg/v_diam)<<endl;

    int nedge = v_lseg.n_rows;
    v_lseg = ledge;


    /*int cntr = 0;
    for (int iseg = 1; iseg <= v_lseg.n_elem; iseg++) {
        int count = 0;
        for (int jseg = 0; jseg < nseg; jseg++) {
            if (iseg == (int) master(jseg))   {
                v_diam(iseg-1) += diam(jseg);
                count += 1;
            }
        }
        v_diam(iseg-1) /= count;
    }*/

    //vessel_order(master,v_lseg);

    net2amira("Interstitial/Branch2Amira.txt","Branch", nnod, nseg, cnode, ista, iend, rseg, conv_to<vec>::from(master));

    /*for (int iseg= 0; iseg < nseg; iseg++)  {
        if (master(iseg) == -1) {
            cout<<flag(iseg)<<"\t"<<q(iseg)<<"\t"<<ista(iseg)<<"\t"<<iend(iseg)<<"\t"<<endl;
            cout<<"\t\t"<<nodtyp(ista(iseg))<<"\t"<<nodtyp(iend(iseg))<<endl;
        }
    }*/

    //cout<<master<<endl;

    return master;
}