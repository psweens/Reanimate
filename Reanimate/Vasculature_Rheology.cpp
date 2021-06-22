//
// Created by sweene01 on 21/06/2020.
//

#include <cmath>

#include "Vasculature.hpp"
#include "spatGraph.hpp"

using namespace reanimate;

void Vasculature::rheolParams()  {

    printText("Assigning rheological parameters");

    // Bifurcation law
    bifpar(0) = 0.964;
    bifpar(1) = 6.98;
    bifpar(2) = -13.29;

    // C - viscosity variable
    cpar(0) = 0.80;
    cpar(1) = -0.075;
    cpar(2) = -11.0;
    cpar(3) = 12.0;


    // In vitro visocity params
    viscpar(0) = 220;
    viscpar(1) = -1.3;
    viscpar(2) = 3.2;
    viscpar(3) = -2.44;
    viscpar(4) = -0.06;
    viscpar(5) = 0.645;

    // Constant visc
    constvisc = 3.;

    // Plasma viscosity (cP)
    vplas = 1.0466;

    // Constant hematocrit
    consthd = 0.4;

    // Mean cell vol.
    mcv = 55.;

}

double Vasculature::viscor(const double &d, const double &hd) {

    double dOff = 2.4;
    double dCrit = 10.5;
    double d50 = 100;
    double eAmp = 1.1;
    double eWidth = 0.03;
    double ePeak = 0.6;
    double eHD = 1.18;
    double wMax = 2.6;

    double wAs = 0.;
    if (dOff < d)   {
        wAs = wMax*(d-dOff)/(d+d50-2*dOff);
    }

    double wPeak = 0.;
    if (d > dOff && d <= dCrit) {
        wPeak = eAmp*(d-dOff)/(dCrit-dOff);
    }
    else if (dCrit < d) {
        wPeak = eAmp*exp(-eWidth*(d-dCrit));
    }

    double wPH = wAs + wPeak*ePeak;
    double wEFF = wAs + wPeak*(1 + hd*eHD);
    double dPH = d - 2*wPH;

    float hdref = 0.45F;
    double C = (cpar(0) + exp(cpar(1)*dPH)) * (-1. + 1./(1. + powf(10.F,cpar(2)) * powf(dPH,cpar(3)))) + 1./(1. + powf(10.F,cpar(2)) * powf(dPH,cpar(3))); // Curvature of relationship between relative apparent viscosity and hematrocrit with tube diameter.
    double eta45 = viscpar(0) * exp(viscpar(1)*dPH) + viscpar(2) + viscpar(3) * exp(viscpar(4) * pow(dPH,viscpar(5))); // mu_0.45
    double hdfac = (powf(1.F - hd,C) - 1.)/(powf(1.F - hdref,C) - 1.); // Discharge haematocrit fraction in mu_vitro equation
    double etaVitro = 1. + (eta45 - 1.) * hdfac;

    double etaVivo = etaVitro*pow(d/(d-2*wEFF),4)*vplas;

    return etaVivo;
}

void Vasculature::dishem(bool &memoryeffects, Network &graph)    {

    //printText("Calculating haematocrit distribution",2, 0);
    //cout<<"dishem"<<endl;
    //timecheck();

    int segfltot = graph.getNseg();
    ivec segs = zeros<ivec>(graph.nodsegm);
    vec flow = zeros<vec>(graph.nodsegm);

    // Assign Hd to boundary nodes
    int isegk{}; // Counts inflowing nodes
    for (int inodbc = 0; inodbc < graph.getNnodbc(); inodbc++) {
        if (graph.nodout(graph.bcnod(inodbc)) >= 1)  {
            graph.hd(graph.nodseg(0, graph.bcnod(inodbc))) = graph.bchd(inodbc);
            isegk += 1;
        }
    }
    //timecheck();

    // Cycle through interior inflowing nodes
    int nodt{},nout{},nrank{};
    double flowsum{},hq{};
    for (int in = 0; in < graph.nnodfl; in++) {

        // If node is an interior node:
        // Store absolute flow value of segment in increasing rank order
        nrank = (int) graph.nodrank(in);
        nodt = (int) graph.nodtyp(nrank);
        nout = (int) graph.nodout(nrank);
        if (nodt >= 2)  {
            for (int i = 0; i < nodt; i++)  {
                segs(i) = graph.nodseg(i,nrank);
                flow(i) = abs(graph.q(segs(i)));
            }
        }

        // 2-segment nodes: nout = 2 could happen in iterations
        if (nodt == 2)  {
            if (nout == 1)  {graph.hd(segs(0)) = graph.hd(segs(1));}
            else if (nout == 2)  {
                graph.hd(segs(0)) = graph.bchd(0);
                graph.hd(segs(1)) = graph.bchd(0);
            }
        }

        // 3-segment nodes: convergent, divergent
        if (nodt == 3)  {
            if (!memoryeffects)  {
                woMemory(nout, nodt, segfltot, segs, flow, graph);
            }
            else    {
                // With memory effects - daughter with: high flow (favourable, f), low flow (unfavourable, u)
                wMemory(nout, nodt, segfltot, segs, flow, graph);
            }
        }

        // No phase separation for nodes with more than three segments
        if (nodt > 3)  {
            if (nout == nodt)   {
                for (int i = 0; i < nodt; i++) {graph.hd(segs(i)) = graph.bchd(0);}
            }
            else    {
                flowsum = 0.;
                hq = 0.;
                for (int i = nout; i < nodt; i++)   {
                    flowsum += abs(graph.q(segs(i)));
                    hq += graph.hd(segs(i)) * abs(graph.q(segs(i)));
                }
                if (nout >= 1)  {
                    for (int i = 0; i < nout; i++)  {
                        graph.hd(segs(i)) = hq / flowsum;
                    }
                }
            }
        }
        if (nodt >= 2)  {
            isegk += nout;
        }
    }
    if(isegk != segfltot) {
        printText("Hd computed in " + to_string(isegk) + " of " + to_string(segfltot) + " segments processed", 4);
    }
    //timecheck();

    double tmp{};
    for (int iseg = 0; iseg < graph.getNseg(); iseg++) {
        tmp = graph.hd(iseg);
        for (int jseg = 0; jseg < graph.edgetyp(iseg); jseg++)    {
            hd(graph.edgeseg(jseg, iseg)) = tmp;
        }
    }
    //timecheck();

}

void Vasculature::woMemory(int &nout, int &nodt, int &segfltot, ivec &segs, vec &flow, Network &graph)    {

    if (nout == 1)  { // Convergent - use conservation
        graph.hd(segs(0)) = (flow(1)*graph.hd(segs(1)) + flow(2)*graph.hd(segs(2))) / flow(0);
    }
    if (nout == 2) { // Divergent - apply empirical law
        double hdd = (1. - graph.hd(segs(2)))/graph.diam(segs(2));
        double diaquot = pow(graph.diam(segs(0))/graph.diam(segs(1)),2);
        double a = bifpar(2)*(diaquot - 1.)/(diaquot + 1.)*hdd;
        double b = 1. + bifpar(1)*hdd;
        double x0 =  bifpar(0)*hdd;
        double qikdash = (flow(0)/flow(2) - x0)/(1. - 2.*x0);
        if (qikdash <= 0.){
            graph.hd(segs(0)) = 0.;
            graph.hd(segs(1)) = graph.hd(segs(2))*flow(2)/flow(1);
        }
        else if (qikdash >= 1.){
            graph.hd(segs(1)) = 0.;
            graph.hd(segs(0)) = graph.hd(segs(2))*flow(2)/flow(0);
        }
        else {
            double rbcrat = 1./(1. + exp(-a-b*log(qikdash/(1.-qikdash))));
            graph.hd(segs(0)) = rbcrat*graph.hd(segs(2))*flow(2)/flow(0);
            graph.hd(segs(1)) = (1.-rbcrat)*graph.hd(segs(2))*flow(2)/flow(1);
        }
    }
    else if (nout == 3) {
        for (int i = 0; i < 3; i++) {
            graph.hd(segs(i)) = graph.bchd(0);
        }
    }

}

void Vasculature::wMemory(int &nout, int &nodt, int &segfltot, ivec &segs, vec &flow, Network &graph)    {


    if (nout == 1)  { // Convergent - use conservation
        graph.hd(segs(0)) = (flow(1)*graph.hd(segs(1)) + flow(2)*graph.hd(segs(2))) / flow(0);
    }
    if (nout == 2) { // Divergent - apply empirical law
        double dp = graph.diam(segs(2));
        double l = graph.lseg(segs(2)); // Distance to previous bifurcation - calculated using edge/vertex graph
        double rfn{};
        if (l < dp) {rfn = 0.;}
        else {rfn = recovfn(l, dp);}
        double x0 = bifpar(0) * (1. - graph.hd(segs(2))) / dp;
        double x0f = x0 * (1. - rfn);
        double x0u = x0 * (1. + rfn);
        double qikdash = (flow(0)/flow(2) - x0f) / (1. - x0u - x0f);

        if (qikdash <= 0.){
            graph.hd(segs(0)) = 0.;
            graph.hd(segs(1)) = graph.hd(segs(2))*flow(2)/flow(1);
        }
        else if (qikdash >= 1.){
            graph.hd(segs(1)) = 0.;
            graph.hd(segs(0)) = graph.hd(segs(2))*flow(2)/flow(0);
        }
        else {
            double A = -6.96 * log(graph.diam(segs(0)) / graph.diam(segs(1))) / dp;
            double B = 1. + 6.98 * (1. - graph.hd(segs(2))) / dp;
            double Ashift = 0.5;
            double Af = A + Ashift * rfn;
            double rbcrat = 1./(1. + exp(-Af-B*log(qikdash/(1.-qikdash))));
            graph.hd(segs(0)) = rbcrat*graph.hd(segs(2))*flow(2)/flow(0);
            graph.hd(segs(1)) = (1.-rbcrat)*graph.hd(segs(2))*flow(2)/flow(1);
        }
    }
    else if (nout == 3) {
        for (int i = 0; i < 3; i++) {
            graph.hd(segs(i)) = bchd(0);
        }
    }

}

// Cell-free layer 90% recovered -> recovfn = 0.1
double Vasculature::recovfn(double &len, double &dp) {

    double omeg = 4.; // Cell-free layer is ~90% recovered at a distance 10*dp -> omeg = ~4.
    return exp(-len / (omeg * dp));

}

