#include "Network.hpp"

using namespace reanimate;

void Network::pictureNetwork(const string &filename, vec vector, bool logdist, int nl, bool nodes, bool segs)  {

    double xs,ys;

    if (logdist)    {vector = log(vector);}

    FILE *ofp;

    double xmin = min(cnode.row(0));
    double xmax = max(cnode.row(0));
    double ymin = min(cnode.row(1));
    double ymax = max(cnode.row(1));

    double picfac = 500./fmax(xmax,ymax);

    ofp = fopen((buildPath + filename).c_str(), "w");
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
    for (int inodbc = 0; inodbc < nnodbc; inodbc++) {
        if (cnode(2,bcnod(inodbc)) > 0.95*alz)    {
            xs = cnode(0,bcnod(inodbc)) - xmin;
            ys = cnode(1,bcnod(inodbc)) - ymin;
            fprintf(ofp, "%g mx %g my m ", xs + 0.5/picfac,ys);
            fprintf(ofp, "(%i) show\n",(int) nodname(bcnod(inodbc)));
        }
    }
    fprintf(ofp, "showpage\n");


    fprintf(ofp, "%%%%Page: %i \n",2);
    //plot vessels using Matlab 'jet' scheme
    double plot{},blue{},green{},red{},vmin{},vmax{},nod1{},nod2{};
    vmin = vector.min();
    vmax = vector.max();
    for(int i = 0; i < (int) vector.n_elem; i++)  {
        if (vmax != vmin)   {plot = (vector(i) - vmin)/(vmax-vmin);}
        else {plot = vector(i) - vmin;}
        blue = min(max(1.5-4*abs(plot-0.25),0.),1.);
        green = min(max(1.5-4*abs(plot-0.5),0.),1.);
        red = min(max(1.5-4*abs(plot-0.75),0.),1.);
        if(plot < 0. || plot > 1.){
            red = 1.0;
            green = 1.0;
            blue = 1.0;
        }
        fprintf(ofp,"%6.3f %6.3f %6.3f c\n",red,green,blue);

        nod1 = ista(i);
        nod2 = iend(i);
        //fprintf(ofp,"1 setlinewidth\n"); // vessels have uniform diameter
        fprintf(ofp,"%g setlinewidth\n",picfac*diam(i)); // vessels are drawn as proportional to their actual diameters
        xs = cnode(0,nod1) - xmin;
        ys = cnode(1,nod1) - ymin;
        fprintf(ofp, "%g mx %g my m ", xs,ys);
        xs = cnode(0,nod2) - xmin;
        ys = cnode(1,nod2) - ymin;
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


void Network::plotContour(const string filename, Network &graph, double maxval, double minval, bool vector, bool overlay, const int xgrid, const int ygrid, const int NL) {

    string file = buildPath + filename;

    int isl1{},isl2{};
    double xs{},ys{},xmax{},ymax{},xmin{},ymin{};
    double red{},green{},blue{},pressureplot{};

    FILE *ofp;

    xmin = min(graph.cnode.row(0));
    xmax = max(graph.cnode.row(0));
    ymin = min(graph.cnode.row(1));
    ymax = max(graph.cnode.row(1));

    double picfac = 500./fmax(xmax,ymax);

    vec xsl0 = zeros<vec>(3);
    vec xsl1 = zeros<vec>(3);
    vec xsl2 = zeros<vec>(3);
    xsl0(0) = xmin;
    xsl0(1) = ymin;
    //xsl0(2) = 0.5 * max(graph.cnode.row(2));
    xsl1(0) = xmax;
    xsl1(1) = ymin;
    //xsl1(2) = 0.5 * max(graph.cnode.row(2));
    xsl2(0) = xmin;
    xsl2(1) = ymax;
    //xsl2(2) = 0.5 * max(graph.cnode.row(2));

    int nlmax = 1;
    if (NL > nlmax) {nlmax = NL;}

    vec cl = zeros<vec>(nlmax+1);
    vec x = zeros<vec>(3);
    mat zv = zeros<mat>(xgrid+1,ygrid+1);


    //Calculate P on a planar slice through the region, specified by three corners and number of points along each edge
    int cntr = 0;
    double p{};
    if (vector) {printText("Generating speed contour plot");}
    else {printText("Generating pressure contour plot");}
    for(isl1 = 1; isl1 <= xgrid; isl1++) {
        for(isl2 = 1; isl2 <= ygrid; isl2++){
            for(int i = 0; i < 3; i++) {x(i) = xsl0(i) + (isl1-1)*(xsl1(i)-xsl0(i))/(xgrid-1) + (isl2-1)*(xsl2(i)-xsl0(i))/(ygrid-1);}
            if (vector) {p = (evalTissVel(x));}
            else {p = evalTissPress(x);}
            zv(isl1-1,isl2-1) = p;
        }
    }

    if (maxval == 0. && minval == 0.)   {
        maxval = zv.max();
        minval = zv.min();
    }
    else {
        zv(find(zv > maxval)).fill(maxval);
        zv(find(zv < minval)).fill(minval);
    }
    double pint = (maxval - minval) / NL;
    printNum("Max. contour value =", maxval);
    printNum("Min. contour value =", minval);

    ofp = fopen(file.c_str(), "w");
    fprintf(ofp, "%%!PS-Adobe-2.0\n");
    fprintf(ofp, "/Times-Roman findfont\n");
    fprintf(ofp, "12 scalefont\n");
    fprintf(ofp, "setfont\n");
    for(int i = 0; i <= NL; i++) cl(i) = minval + (i)*pint;

    shadeContour(ofp, xgrid, ygrid, picfac, const_cast<int &>(NL), pint, xmin, xmax, ymin, ymax, cl, zv);

    // Network projection
    if (overlay)    {
        for(int iseg = 0; iseg < graph.getNseg(); iseg++){
            pressureplot = (graph.segpress(iseg) - minval)/(maxval - minval);
            blue = fmin(fmax(1.5-4*fabs(pressureplot-0.25),0.),1.);
            green = fmin(fmax(1.5-4*fabs(pressureplot-0.5),0.),1.);
            red = fmin(fmax(1.5-4*fabs(pressureplot-0.75),0.),1.);
            if(pressureplot < 0. || pressureplot > 1.){
                red = 1.0;
                green = 1.0;
                blue = 1.0;
            }

            fprintf(ofp,"%6.3f %6.3f %6.3f setrgbcolor\n",red,green,blue);
            //		fprintf(ofp,"1 setlinewidth\n"); // vessels have uniform diameter
            fprintf(ofp,"%g setlinewidth\n",picfac*graph.diam(iseg)); // vessels are drawn as proportional to their actual diameters
            xs = graph.cnode(0,graph.ista(iseg)) - xmin;
            ys = graph.cnode(1,graph.ista(iseg)) - ymin;
            fprintf(ofp, "%g mx %g my m ", xs,ys);
            xs = graph.cnode(0,graph.iend(iseg)) - xmin;
            ys = graph.cnode(1,graph.iend(iseg)) - ymin;
            fprintf(ofp, "%g mx %g my l ", xs,ys);
            fprintf(ofp, "stroke\n");
        }
    }


    // Create a colour bar
    double cbbox = 21.;    //size of boxes
    double cbx = 550;  //origin of color bar
    double cby = 50;   //origin of color bar
    vec bluevect = zeros<vec>(NL/3+1);
    vec greenvect = zeros<vec>(NL/3+1);
    vec redvect = zeros<vec>(NL/3+1);
    vec pvals = zeros<vec>(NL/3+1);
    for(int k = 0; k <= int(NL/3.); k++){
        double xz = 3*float(k)/float(NL);
        bluevect(k) = min(max(1.5-4*abs(xz-0.25), 0.), 1.);
        greenvect(k)= min(max(1.5-4*abs(xz-0.5), 0.), 1.);
        redvect(k)  = min(max(1.5-4*abs(xz-0.75), 0.), 1.);
        pvals(k)= minval + 3*k*(maxval - minval)/NL;
    }
    fprintf(ofp,"%d setlinewidth\n",1);
    for(int k = 0; k <= int(NL/3.); k++){
        fprintf(ofp, "%f %f %f setrgbcolor\n",redvect(k),greenvect(k),bluevect(k));
        fprintf(ofp, "n %g %g m %g %g l %g %g l %g %g l cf\n",cbx,cby+k*cbbox,cbx+cbbox,cby+k*cbbox,cbx+cbbox,cby+(k+1)*cbbox,cbx,cby+(k+1.)*cbbox);
        fprintf(ofp, "%g %f m 0 0 0 setrgbcolor (%.3g) show\n",cbx+cbbox*1.1,cby+cbbox*(k-0.1),pvals(k));
    }
    fprintf(ofp, "n %g %g m %g %g l %g %g l %g %g l cs\n",cbx,cby,cbx+cbbox,cby,cbx+cbbox,cby+cbbox*(NL/3.+1.),cbx,cby+cbbox*(NL/3.+1.));

    fprintf(ofp, "showpage\n");
    fclose(ofp);

}


void Network::shadeContour(FILE *ofp, const int &m, const int &n, double &scalefac, int &nl, const double pint, const double &xmin, const double &xmax, const double &ymin, const double &ymax, const vec &cl, const mat &zv) {

    fprintf(ofp, "%%!PS\n");
    int iwsp, in1, in2, in3, in4, inh;
    const int ncmax = 100000;
    double xz = 0.0, xv12 = 0.0, yv12 = 0.0, xv23 = 0.0, yv23 = 0.0, xv34 = 0.0, yv34 = 0.0, xv41 = 0.0, yv41 = 0.0;
    double cx, cy, cbx, cby, cbbox, eps;

    cx = 50;    //25; //origin of contour plot
    cy = 100;   //50;//origin of contour plot
    cbbox = 21.;    //size of boxes
    cbx = 550;  //origin of color bar
    cby = 50;   //origin of color bar

    // vectors were indexed from 0 to nl --> nl+1 entries so we now index from 0 to nl+1
    vec xv = zeros<vec>(m + 1);
    vec yv = zeros<vec>(n + 1);
    vec xvv = zeros<vec>(5);
    vec yvv = zeros<vec>(5);
    vec red = zeros<vec>(nl + 1);
    vec green = zeros<vec>(nl + 1);
    vec blue = zeros<vec>(nl + 1);
    mat wsp = zeros<mat>(ncmax + 1, 5);
    imat corners = zeros<imat>(5, 2);
    vec dd = zeros<vec>(5);
    ivec ii = zeros<ivec>(5);
    ivec jj = zeros<ivec>(5);
    corners(1, 0) = 0;
    corners(1, 1) = 0;
    corners(2, 0) = 1;
    corners(2, 1) = 0;
    corners(3, 0) = 1;
    corners(3, 1) = 1;
    corners(4, 0) = 0;
    corners(4, 1) = 1;


    for (int i = 1; i <= m; i++) xv(i) = xmin + (i - 1) * (xmax - xmin) / (m - 1) - xmin;
    for (int j = 1; j <= n; j++) yv(j) = ymin + (j - 1) * (ymax - ymin) / (n - 1) - ymin;

    fprintf(ofp, "/mx {%g sub %g mul %g add} def\n", xmin - xmin, scalefac, cx);
    fprintf(ofp, "/my {%g sub %g mul %g add} def\n", ymin - ymin, scalefac, cy);
    fprintf(ofp, "/m {moveto} def\n");
    fprintf(ofp, "/l {lineto} def\n");
    fprintf(ofp, "/n {newpath} def\n");
    fprintf(ofp, "/s {stroke} def\n");
    fprintf(ofp, "/cf {closepath fill} def\n");
    fprintf(ofp, "/cs {closepath stroke} def\n");
    fprintf(ofp, "0.5 setlinewidth\n");
    fprintf(ofp, "/Times-Roman findfont\n");
    fprintf(ofp, "12 scalefont\n");
    fprintf(ofp, "setfont\n");

    //eps defines extra width added to colored rectangles.
    //otherwise, small lines of different color appear between rectangles
    eps = fmax(xmax - xmin, ymax - ymin) * 1.e-3;

    //Set up colors using Matlab 'jet' scheme
    for (int k = 0; k <= nl; k++) {
        xz = float(k) / float(nl);
        blue(k) = fmin(fmax(1.5 - 4 * fabs(xz - 0.25), 0.), 1.);
        green(k) = fmin(fmax(1.5 - 4 * fabs(xz - 0.5), 0.), 1.);
        red(k) = fmin(fmax(1.5 - 4 * fabs(xz - 0.75), 0.), 1.);
    }
    //Color whole region with lowest color
    fprintf(ofp, "%f %f %f setrgbcolor\n", red(0), green(0), blue(0));
    fprintf(ofp, "n %g mx %g my m %g mx %g my l %g mx %g my l %g mx %g my l cf\n",
            xmin - xmin, ymin - ymin, xmax - xmin, ymin - ymin, xmax - xmin, ymax - ymin, xmin - xmin, ymax - ymin);
    //Analyze each rectangle separately. Overwrite lower colors
    iwsp = 0;
    for (int k = 0; k <= nl; k++) {
        fprintf(ofp, "%f %f %f setrgbcolor\n", red(k), green(k), blue(k));
        for (int i = 1; i < m; i++)
            for (int j = 1; j < n; j++) {
                in1 = 0;
                for (int in = 0; in <= 4; in++) {
                    ii(in) = i + corners(in, 0);
                    jj(in) = j + corners(in, 1);
                    dd(in) = zv(ii(in)-1, jj(in)-1) - cl(k);
                    if (dd(in) >= 0.) in1 = in;    //find a corner this color or higher
                }
                inh = 0;
                if (k < nl)
                    for (int in = 0; in <= 4; in++)
                        if (zv(ii(in)-1, jj(in)-1) - cl(k + 1) < 0.)
                            inh = 1;//check that not all corners are above the next contour
                if (in1 > 0 && inh > 0) {
                    in2 = in1 % 4 + 1;
                    in3 = in2 % 4 + 1;
                    in4 = in3 % 4 + 1;
                    for (int in = 1; in <= 4; in++) {
                        xvv(in) = xv(ii(in));
                        yvv(in) = yv(jj(in));
                        if (ii(in) == i + 1) xvv(in) += eps; else xvv(in) -= eps;
                        if (jj(in) == j + 1) yvv(in) += eps; else yvv(in) -= eps;
                    }
                    if (dd(in1) != dd(in2)) {
                        xv12 = (dd(in1) * xv(ii(in2)) - dd(in2) * xv(ii(in1))) / (dd(in1) - dd(in2));
                        yv12 = (dd(in1) * yv(jj(in2)) - dd(in2) * yv(jj(in1))) / (dd(in1) - dd(in2));
                    }
                    if (dd(in2) != dd(in3)) {
                        xv23 = (dd(in2) * xv(ii(in3)) - dd(in3) * xv(ii(in2))) / (dd(in2) - dd(in3));
                        yv23 = (dd(in2) * yv(jj(in3)) - dd(in3) * yv(jj(in2))) / (dd(in2) - dd(in3));
                    }
                    if (dd(in3) != dd(in4)) {
                        xv34 = (dd(in3) * xv(ii(in4)) - dd(in4) * xv(ii(in3))) / (dd(in3) - dd(in4));
                        yv34 = (dd(in3) * yv(jj(in4)) - dd(in4) * yv(jj(in3))) / (dd(in3) - dd(in4));
                    }
                    if (dd(in4) != dd(in1)) {
                        xv41 = (dd(in4) * xv(ii(in1)) - dd(in1) * xv(ii(in4))) / (dd(in4) - dd(in1));
                        yv41 = (dd(in4) * yv(jj(in1)) - dd(in1) * yv(jj(in4))) / (dd(in4) - dd(in1));
                    }
                    fprintf(ofp, "n %g mx %g my m ", xvv(in1), yvv(in1));
                    if (dd(in2) > 0) {                                                    //corners 1,2 are this color
                        fprintf(ofp, "%g mx %g my l ", xvv(in2), yvv(in2));
                        if (dd(in3) > 0) {                                                //corners 1,2,3 are this color
                            fprintf(ofp, "%g mx %g my l ", xvv(in3), yvv(in3));
                            if (dd(in4) > 0)
                                fprintf(ofp, "%g mx %g my l ", xvv(in4),
                                        yvv(in4));        //corners 1,2,3,4 are this color
                            else {                                                        //corners 1,2,3,not 4 are this color
                                fprintf(ofp, "%g mx %g my l ", xv34, yv34);
                                fprintf(ofp, "%g mx %g my l ", xv41, yv41);
                                iwsp += 1;
                                wsp(iwsp, 1) = xv34;
                                wsp(iwsp, 2) = yv34;
                                wsp(iwsp, 3) = xv41;
                                wsp(iwsp, 4) = yv41;
                            }
                        } else {                                                            //corners 1,2,not 3 are this color
                            fprintf(ofp, "%g mx %g my l ", xv23, yv23);
                            iwsp += 1;
                            wsp(iwsp, 1) = xv23;
                            wsp(iwsp, 2) = yv23;
                            if (dd(in4) >
                                0) {                                            //corners 1,2,not 3,4 are this color
                                fprintf(ofp, "%g mx %g my l ", xv34, yv34);
                                wsp(iwsp, 3) = xv34;
                                wsp(iwsp, 4) = yv34;
                                fprintf(ofp, "%g mx %g my l ", xvv(in4), yvv(in4));
                            } else {
                                fprintf(ofp, "%g mx %g my l ", xv41,
                                        yv41);                //corners 1,2,not 3,not 4 are this color
                                wsp(iwsp, 3) = xv41;
                                wsp(iwsp, 4) = yv41;
                            }
                        }
                    } else {                                                                //corners 1,not 2 are this color
                        fprintf(ofp, "%g mx %g my l ", xv12, yv12);
                        iwsp += 1;
                        wsp(iwsp, 1) = xv12;
                        wsp(iwsp, 2) = yv12;
                        if (dd(in3) >
                            0) {                                                //corners 1,not 2,3 are this color
                            fprintf(ofp, "%g mx %g my l ", xv23, yv23);
                            wsp(iwsp, 3) = xv23;
                            wsp(iwsp, 4) = yv23;
                            fprintf(ofp, "%g mx %g my l ", xvv(in3), yvv(in3));
                            if (dd(in4) > 0)
                                fprintf(ofp, "%g mx %g my l ", xvv(in4),
                                        yvv(in4));        //corners 1,not 2,3,4 are this color
                            else {                                                        //corners 1,not 2,3,not 4 are this color
                                fprintf(ofp, "%g mx %g my l ", xv34, yv34);
                                fprintf(ofp, "%g mx %g my l ", xv41, yv41);
                                iwsp += 1;
                                wsp(iwsp, 1) = xv34;
                                wsp(iwsp, 2) = yv34;
                                wsp(iwsp, 3) = xv41;
                                wsp(iwsp, 4) = yv41;
                            }
                        } else {                                                            //corners 1,not 2,not 3 are this color
                            if (dd(in4) >
                                0) {                                            //corners 1,not 2,not 3,4 are this color
                                fprintf(ofp, "%g mx %g my l ", xv34, yv34);
                                wsp(iwsp, 3) = xv34;
                                wsp(iwsp, 4) = yv34;
                                fprintf(ofp, "%g mx %g my l ", xvv(in4), yvv(in4));
                            } else {
                                fprintf(ofp, "%g mx %g my l ", xv41,
                                        yv41);                //corners 1,not 2,not 3,not 4 are this color
                                wsp(iwsp, 3) = xv41;
                                wsp(iwsp, 4) = yv41;
                            }
                        }
                    }
                    if (iwsp > ncmax - 4) printf("*** Error: ncmax too small in contr\n");
                    fprintf(ofp, "cf\n");
                }
            }
    }
    //Now outline contours
    /*fprintf(ofp, "0 0 0 setrgbcolor\n");//black
    for(int in=1; in<=iwsp; in++) fprintf(ofp, "n %g mx %g my m %g mx %g my l s\n",
                                      wsp(in,1),wsp(in,2),wsp(in,3),wsp(in,4));*/
    //Draw a box
    //fprintf(ofp, "0 0 0 setrgbcolor\n");//black
    //fprintf(ofp, "n %g mx %g my m %g mx %g my l %g mx %g my l %g mx %g my l cs\n",
    //xmin-xmin,ymin-ymin,xmax-xmin,ymin-ymin,xmax-xmin,ymax-ymin,xmin-xmin,ymax-ymin);



    /*fprintf(ofp,"%d setlinewidth\n",1);
    for(int k=0; k<=nl; k++){
        fprintf(ofp, "%f %f %f setrgbcolor\n",red(k),green(k),blue(k));
        fprintf(ofp, "n %g %g m %g %g l %g %g l %g %g l cf\n",
                cbx,cby+k*cbbox,cbx+cbbox,cby+k*cbbox,cbx+cbbox,cby+(k+1)*cbbox,cbx,cby+(k+1)*cbbox);
        //if(k>0) {
            fprintf(ofp, "%g %f m 0 0 0 setrgbcolor (%.3g) show\n",cbx+cbbox*1.1,cby+cbbox*(k-0.1),cl(k));
            cout<<cl(k)<<endl;
        //}

    }
    fprintf(ofp, "n %g %g m %g %g l %g %g l %g %g l cs\n",
            cbx,cby,cbx+cbbox,cby,cbx+cbbox,cby+cbbox*(nl+1),cbx,cby+cbbox*(nl+1));*/

}
