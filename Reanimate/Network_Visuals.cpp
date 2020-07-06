#include "Network.hpp"

using namespace reanimate;

void Network::pictureNetwork(const string &filename, vec vector, bool logdist, int nl, bool nodes, bool segs)  {

    double xs,ys;

    if (logdist)    {
        vector = log(vector);
    }

    FILE *ofp;

    double xmin = min(cnode.row(0));
    double xmax = max(cnode.row(0));
    double ymin = min(cnode.row(1));
    double ymax = max(cnode.row(1));

    double picfac = 500./fmax(xmax,ymax);

    ofp = fopen(filename.c_str(), "w");
    fprintf(ofp, "%%!PS\n");
    fprintf(ofp, "%%%%Pages: 3\n");
    fprintf(ofp, "/mx {%g mul 50 add} def\n",picfac);
    fprintf(ofp, "/my {%g mul 100 add} def\n",picfac);
    fprintf(ofp, "/m {moveto} def\n");
    fprintf(ofp, "/l {lineto} def\n");
    fprintf(ofp, "/n {newpath} def\n");
    fprintf(ofp, "/c {setrgbcolor} def\n");
    fprintf(ofp, "/cf {closepath fill} def\n");
    fprintf(ofp, "/cs {closepath stroke} def\n");
    fprintf(ofp, "/Times-Roman findfont\n");
    fprintf(ofp, "6 scalefont\n");
    fprintf(ofp, "setfont\n");

    fprintf(ofp, "newpath\n");
    fprintf(ofp, "%g mx %g my m\n",xmin-xmin,ymin-ymin);
    fprintf(ofp, "%g mx %g my l\n",xmax-xmin,ymin-ymin);
    fprintf(ofp, "%g mx %g my l\n",xmax-xmin,ymax-ymin);
    fprintf(ofp, "%g mx %g my l\n",xmin-xmin,ymax-ymin);
    fprintf(ofp, "closepath\n");
    fprintf(ofp, "stroke\n");

    fprintf(ofp, "%%%%Page: %i \n",1);

    // Plot vessels (in red)
    for (int iseg = 0; iseg < nseg; iseg++){
        fprintf(ofp,"1 0 0 setrgbcolor\n"); //red
        fprintf(ofp,"%g setlinewidth\n",picfac*diam(iseg));
        xs = cnode(0,ista(iseg)) - xmin;
        ys = cnode(1,ista(iseg)) - ymin;
        fprintf(ofp, "%g mx %g my m ", xs,ys);
        xs = cnode(0,iend(iseg)) - xmin;
        ys = cnode(1,iend(iseg)) - ymin;
        fprintf(ofp, "%g mx %g my l ", xs,ys);
        fprintf(ofp, "stroke\n");
    }
    // Label boundary nodes in black
    fprintf(ofp,"0 0 0 setrgbcolor\n");//black
    for (int inod = 0; inod < nnod; inod++){
        if (nodtyp(inod) == 1) {
            for (int iseg = 0; iseg < nseg; iseg++)    {
                if (inod==ista(iseg) || inod==iend(iseg))    {
                    xs = cnode(0,inod) - xmin;
                    ys = cnode(1,inod) - ymin;
                    fprintf(ofp, "%g mx %g my m ", xs + 0.5/picfac,ys);
                    fprintf(ofp, "(%i) show\n",(int) nodname(inod));
                }
            }
        }
    }
    fprintf(ofp, "showpage\n");


    fprintf(ofp, "%%%%Page: %i \n",2);
    //plot vessels using Matlab 'jet' scheme
    for(int i = 0; i < (int) vector.n_elem; i++)  {
        double plot = (vector(i) - vector.min())/(vector.max()-vector.min());
        double blue = min(max(1.5-4*abs(plot-0.25),0.),1.);
        double green = min(max(1.5-4*abs(plot-0.5),0.),1.);
        double red = min(max(1.5-4*abs(plot-0.75),0.),1.);
        if(plot < 0. || plot > 1.){
            red = 1.0;
            green = 1.0;
            blue = 1.0;
        }
        fprintf(ofp,"%6.3f %6.3f %6.3f c\n",red,green,blue);

        //fprintf(ofp,"1 setlinewidth\n"); // vessels have uniform diameter
        fprintf(ofp,"%g setlinewidth\n",picfac*diam(i)); // vessels are drawn as proportional to their actual diameters
        xs = cnode(0,ista(i)) - xmin;
        ys = cnode(1,ista(i)) - ymin;
        fprintf(ofp, "%g mx %g my m ", xs,ys);
        xs = cnode(0,iend(i)) - xmin;
        ys = cnode(1,iend(i)) - ymin;
        fprintf(ofp, "%g mx %g my l ", xs,ys);
        fprintf(ofp, "stroke\n");
    }

    //label nodes in black
    if (nodes)  {
        fprintf(ofp,"0 0 0 setrgbcolor\n");//black
        for(int inod = 0; inod < nnod; inod++) {
            xs = cnode(0,inod)-xmin;
            ys = cnode(1,inod)-ymin;
            fprintf(ofp, "%g mx %g my m ", xs + 0.5/picfac,ys);
            fprintf(ofp, "(%lld) show\n",nodname(inod));
        }
    }

    //label segments in blue
    if (segs)   {
        fprintf(ofp,"0 0 1 setrgbcolor\n");//blue
        for(int iseg = 0; iseg < nseg; iseg++) {
            xs = (cnode(0,ista(iseg)) + cnode(0,iend(iseg))-2*xmin)/2.;
            ys = (cnode(1,ista(iseg)) + cnode(1,iend(iseg))-2*ymin)/2.;
            fprintf(ofp, "%g mx %g my m ", xs + 0.5*picfac,ys);
            fprintf(ofp, "(%.1f) show\n",vector(iseg));
        }
    }

    // Create a colour bar
    // nl  - the number of divisions in the colour bar
    vec redvect = zeros<vec>(nl+1);
    vec greenvect = zeros<vec>(nl+1);
    vec bluevect = zeros<vec>(nl+1);
    vec vals = zeros<vec>(nl+1);
    for(int k = 0; k <= nl; k++){
        float xz = float(k)/float(nl);
        bluevect(k) = min(max(1.5-4*abs(xz-0.25), 0.), 1.);
        greenvect(k)= min(max(1.5-4*abs(xz-0.5), 0.), 1.);
        redvect(k)  = min(max(1.5-4*abs(xz-0.75), 0.), 1.);
        vals(k)= vector.min() + k*(vector.max()-vector.min())/nl;
    }

    fprintf(ofp,"%d setlinewidth\n",1);
    double cbx = 525.; // Origin of colour bar
    double cby = 100.;   // Origin of colour bar
    double cbbox = 21.; // Size of boxes
    for(int k = 0; k <= nl; k++){
        fprintf(ofp, "%f %f %f setrgbcolor\n",redvect(k),greenvect(k),bluevect(k));
        fprintf(ofp, "n %g %g m %g %g l %g %g l %g %g l cf\n",
                cbx,cby+k*cbbox,cbx+cbbox,cby+k*cbbox,cbx+cbbox,cby+(k+1)*cbbox,cbx,cby+(k+1)*cbbox);
        fprintf(ofp, "%g %f m 0 0 0 setrgbcolor (%g) show\n",cbx+cbbox*1.1,cby+cbbox*(k-0.1),vals(k));
    }
    fprintf(ofp, "n %g %g m %g %g l %g %g l %g %g l cs\n",
            cbx,cby,cbx+cbbox,cby,cbx+cbbox,cby+cbbox*(nl+1),cbx,cby+cbbox*(nl+1));

    fprintf(ofp, "showpage\n");


    fclose(ofp);
}