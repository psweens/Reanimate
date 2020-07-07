#include "Network.hpp"

using namespace reanimate;

void Network::edgeNetwork() {

    cout<<"Calculating edge / vertex network ... "<<endl;

    int order = 1;
    int edgIdx = 0;
    vec ledge = zeros<vec>(nseg);
    ivec flag = zeros<ivec>(nseg);
    edgeLabels = zeros<ivec>(nseg);
    edgeLabels.fill(-1);
    elseg = zeros<vec>(nseg);
    ediam = zeros<vec>(nseg);

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
                        edgeLabels(iseg) = feedNod(1, cntr);
                        order += 1;
                        cntr += 1;
                    }
                    else if (iend(iseg) == inod)    {
                        flag(iseg) = 1;
                        feedNod(0,cntr) = ista(iseg);
                        feedNod(1,cntr) = order;
                        edgeLabels(iseg) = feedNod(1, cntr);
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
                for (int i = 0; i < (int) feedNod.n_cols; i++)  {
                    edgIdx = 0;
                    for (int iseg = 0; iseg < nseg; iseg++) {
                        if (ista(iseg) == feedNod(0,i) && flag(iseg) == 0 && nodtyp(feedNod(0,i)) == 2)   {
                            edgeLabels(iseg) = feedNod(1, i);
                            flag(iseg) = 2;
                            feedNod(0,i) = iend(iseg);
                            iseg = 0;
                        }
                        else if (iend(iseg) == feedNod(0,i) && flag(iseg) == 0 && nodtyp(feedNod(0,i)) == 2)   {
                            edgeLabels(iseg) = feedNod(1, i);
                            flag(iseg) = 2;
                            feedNod(0,i) = ista(iseg);
                            iseg = 0;
                        }
                    }
                    edgIdx += 1;
                }
            }

        }
    }
    ivec edgeIdx = unique(edgeLabels);
    for (int iseg = 0; iseg < (int) edgeIdx.n_elem; iseg++)  {
        double d = mean(diam(find(edgeLabels == edgeIdx(iseg))));
        double l = mean(lseg(find(edgeLabels == edgeIdx(iseg))));
        for (int jseg = 0; jseg < nseg; jseg++) {
            if (edgeIdx(iseg) == edgeLabels(jseg))    {
                ediam(jseg) = d;
                elseg(jseg) = l;
            }
        }
    }

    cout << "\tEdge lengths: " << mean(elseg) << " pm " << stddev(elseg) << endl;
    cout << "\tEdge diams: " << mean(ediam) << " pm " << stddev(ediam) << endl;
    cout << "\tLength / diam ratio: " << mean(elseg / ediam) << " pm " << stddev(elseg / ediam) << endl;

    nedge = elseg.n_rows;

}
