//
//  vorder_vessel.cpp
//  Vascular-Flow
//
//  Created by Paul Sweeney on 25/10/2016.
//  Copyright © 2016 Paul Sweeney - University College London. All rights reserved.
//

#include <stdio.h>
#include "global_variables.h"
#include "output_data.h"

void vessel_vorder(ivec master,vec v_lseg)  {
    
    vorder = zeros<ivec>(nseg);
    vorder.fill(-1);
    
    // Locate penetrating arterioles by stripping back SMA classification
    uvec temp_typ = vesstyp;
    ivec nod = zeros<ivec>(nnod);
    for (int iseg = 0; iseg < nseg; iseg++) {
        if (vesstyp(iseg) == 1)   {
            nod(ista(iseg)) += 1;
            nod(iend(iseg)) += 1;
        }
    }
    nod(find(nodtyp == 1)).fill(0);
    for (int iseg = 0; iseg < nseg; iseg++) {
        if (diam(iseg) < 10. && (nod(ista(iseg)) == 1 || nod(iend(iseg)) == 1))   {
            for (int jseg = 0; jseg < nseg; jseg++) {
                if (master(iseg) == master(jseg))   {
                    temp_typ(jseg) = 2;
                }
            }
        }
    }
    
    nod.zeros();
    for (int iseg = 0; iseg < nseg; iseg++) {
        if (temp_typ(iseg) == 1)   {
            nod(ista(iseg)) += 1;
            nod(iend(iseg)) += 1;
        }
    }
    nod(find(nodtyp == 1)).fill(0);
    // Remove isolated branches
    int cntr = 0;
    for (int iseg = 0; iseg < nseg; iseg++) {
        uvec temp = find(master == iseg);
        for (int i = 0; i < temp.n_elem; i++)   {
            if (nod(ista(i)) == 1 || nod(iend(i)) == 1) {
                cntr += 1;
            }
        }
        if (cntr == 2)  {
            temp_typ(temp).fill(2);
        }
    }
    // Accounting for arterioles of nodtyp > 1
    for (int inod = 0; inod < nnod; inod++) {
        if (nodtyp(inod) - nod(inod) == 1 && nodtyp(inod) != 1) {
            nod(inod) = 1;
        }
    }
    
    // Define penetrating arterioles
    vorder(find(temp_typ == 1)).fill(0);
    
    // Define vessel vorders
    for (int iseg = 0; iseg < nseg; iseg++) {
        if (nod(ista(iseg)) == 1)   {
            for (int jseg = 0; jseg < nseg; jseg++) {
                if (iseg != jseg && vorder(jseg) == -1 && (ista(iseg) == ista(jseg) || ista(iseg) == iend(jseg)))   {
                    vorder(find(master == master(jseg))).fill(1);
                }
            }
        }
        else if (nod(iend(iseg)) == 1)   {
            for (int jseg = 0; jseg < nseg; jseg++) {
                if (iseg != jseg && vorder(jseg) == -1 && (iend(iseg) == ista(jseg) || iend(iseg) == iend(jseg)))   {
                    vorder(find(master == master(jseg))).fill(1);
                }
            }
        }
    }
    
    for (int i = 2; i <= 20; i++)    {
        
        nod.zeros();
        for (int iseg = 0; iseg < nseg; iseg++) {
            if (vorder(iseg) > -1)   {
                nod(ista(iseg)) += 1;
                nod(iend(iseg)) += 1;
            }
        }
        nod(find(nodtyp == 1)).fill(0);
        
        for (int iseg = 0; iseg < nseg; iseg++) {
            if (nod(ista(iseg)) == 1)   {
                for (int jseg = 0; jseg < nseg; jseg++) {
                    if (iseg != jseg && vorder(jseg) == -1 && (ista(iseg) == ista(jseg) || ista(iseg) == iend(jseg)) && vesstyp(jseg) != 3 && diam(jseg) < 10.)   {
                        vorder(find(master == master(jseg))).fill(i);
                    }
                }
            }
            else if (nod(iend(iseg)) == 1)   {
                for (int jseg = 0; jseg < nseg; jseg++) {
                    if (iseg != jseg && vorder(jseg) == -1 && (iend(iseg) == ista(jseg) || iend(iseg) == iend(jseg)) && vesstyp(jseg) != 3 && diam(jseg) < 10.0)   {
                        vorder(find(master == master(jseg))).fill(i);
                    }
                }
            }
        }
    }

    
    //vorder(find(vesstyp == 2 && vorder == -1)).fill(4);
    for (int iseg = 0; iseg < nseg; iseg++) {
        if (vesstyp(iseg) == 1 && (cnode(2,ista(iseg)) > 375. || cnode(2,iend(iseg)) > 375.))   {
            //vorder(iseg) = 5;
        }
    }
    //vorder(find(vesstyp == 3 && vorder == -1)).fill(6);
    net2amira("vorder2Amira.txt","VessTyp", nnod, nseg, cnode, ista, iend, conv_to<vec>::from(rseg.row(0)), conv_to<vec>::from(vorder));
    vesstyp = conv_to<uvec>::from(vorder);
    
    nod.zeros();
    for (int iseg = 0; iseg < nseg; iseg++) {
        if (vesstyp(iseg) == 3) {
            nod(ista(iseg)) += 1;
            nod(iend(iseg)) += 1;
        }
    }
    for (int inodbc = 0; inodbc < nnodbc; inodbc++) {
        nod(bcnod(inodbc)) = 0;
    }
    double cnt_this = 0;
    double vorder_num = 0;
    for (int iseg = 0; iseg < nseg; iseg++) {
        if (vesstyp(iseg) == 2 && vorder(iseg) != -1 && (nod(ista(iseg)) == 1 || nod(iend(iseg)) == 1)) {
            cnt_this += 1;
            vorder_num += vorder(iseg);
        }
    }
    //cout<<vorder<<endl;

    printf("Average No. of Branches = %.1f\n",vorder_num/cnt_this);


}
