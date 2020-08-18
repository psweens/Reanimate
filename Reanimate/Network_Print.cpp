#include "Network.hpp"

using namespace reanimate;

void Network::printNetwork(const string& filename, bool resetDomain)    {

    FILE *ofp;

    double minx{},miny{},minz{};
    if (resetDomain)    {

        // Update tissue dimensions
        alx = max(cnode.row(0));
        aly = max(cnode.row(1));
        alz = max(cnode.row(2));
        minx = min(cnode.row(0));
        miny = min(cnode.row(1));
        minz = min(cnode.row(2));

        cnode.row(0) -= minx;
        cnode.row(1) -= miny;
        cnode.row(2) -= minz;
        alx -= minx;
        aly -= miny;
        alz -= minz;

    }


    ofp = fopen((buildPath + filename).c_str(),"w");
    fprintf(ofp,"%s",networkName.c_str());
    fprintf(ofp,"%f %f %f  Box dimensions in microns \n",alx-minx,aly-miny,alz-minz);
    fprintf(ofp,"%i %i %i  No. of tissue points in x,y,z directions \n",mxx,myy,mzz);
    fprintf(ofp,"%f    Outer bound distance \n",lb);
    fprintf(ofp,"%f    Max. segment length \n",maxl);
    fprintf(ofp,"%i    Max. segments per node \n",nodsegm);
    fprintf(ofp,"%i    Total number of segments\n",nseg);
    fprintf(ofp,"Segname Type Start End Diam Flow[qL/min] Hd\n");
    for(int iseg = 0; iseg < nseg; iseg++)     {
        fprintf(ofp,"%lli %lli %lli %lli %f %f %f\n",segname(iseg),vesstyp(iseg),segnodname(0,iseg),segnodname(1,iseg),diam(iseg),q(iseg),hd(iseg));
    }
    fprintf(ofp,"%i Total Number of nodes\n", nnod);
    fprintf(ofp,"Nodname x y z\n");
    for(int inod = 0; inod < nnod; inod++)  {
        fprintf(ofp, "%lli %f %f %f %f\n",nodname(inod),cnode(0,inod)-minx,cnode(1,inod)-miny,cnode(2,inod)-minz,nodpress(inod));
    }
    fprintf(ofp,"%i Total Number of boundary nodes\n", nnodbc);
    fprintf(ofp,"Bcnodname Bctyp Bcprfl BcHd\n");
    for (int inodbc = 0; inodbc < nnodbc; inodbc++)    {
        fprintf(ofp,"%lli %lli %f %f\n",bcnodname(inodbc),bctyp(inodbc),bcprfl(inodbc),bchd(inodbc));
    }

    fclose(ofp);

}

void Network::printHistogram(const string &filename, mat &data, const field<string> &headers)   {

    FILE *ofp;
    ofp = fopen(filename.c_str(),"w");
    fprintf(ofp,"Histogram Data\n");

    data = conv_to<mat>::from(data);
    for (int i = 0; i < (int) data.n_cols; i++)   {

        double dmin = min(data.col(i));
        double dmax = max(data.col(i));
        double dmean = mean(data.col(i));
        double dSD = stddev(data.col(i));

        int range = round(dmax)-round(dmin)+1;
        vec drange = linspace<vec>(round(dmin),round(dmax),range);
        vec dhist = conv_to<vec>::from(hist(data.col(i),drange));

        // Output as percentages
        dhist /= (float) data.n_rows;

        fprintf(ofp,"%s data:\n",headers(i,0).c_str());
        fprintf(ofp,"Mean = %g , S.D. = %g , min = %g , max  = %g\n", dmean,dSD,dmin,dmax);
        fprintf(ofp,"value  %% cumul. %%\n");
        for (int j = 0; j < (int) drange.n_elem; j++)  {
            fprintf(ofp,"%f \t %f\n", drange(j), dhist(j));
        }

    }

    fclose(ofp);

}