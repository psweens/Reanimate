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
    consthd = 0.45;

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
    //vplas *= 1e-9; // Viscosity of plasma in kg/mu s

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

    printText("Calculating haematocrit distribution",2, 0);

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


    // Cycle through interior inflowing nodes
    for (int in = 0; in < graph.nnodfl; in++) {

        // If node is an interior node:
        // Store associated segment indices in increasing rank order
        // Store absolute flow value of segment in increasing rank order
        int nodt = (int) graph.nodtyp(graph.nodrank(in));
        int nout = (int) graph.nodout(graph.nodrank(in));
        if (nodt >= 2)  {
            for (int i = 0; i < nodt; i++)  {
                segs(i) = graph.nodseg(i, graph.nodrank(in));
                flow(i) = abs(graph.q(segs(i)));
            }
        }

        // 2-segment nodes: nout = 2 could happen in iterations
        if (nodt == 2)  {
            if (nout == 1)  {
                graph.hd(segs(0)) = graph.hd(segs(1));
            }
            else if (nout == 2)  {
                graph.hd(segs(0)) = graph.bchd(0);
                graph.hd(segs(1)) = graph.bchd(0);
            }
        }


        // 3-segment nodes: convergent, divergent
        if (!memoryeffects)  {
            woMemory(nout, nodt, segfltot, segs, flow, graph);
        }
        else    {
            // With memory effects - daughter with: high flow (favourable, f), low flow (unfavourable, u)
            wMemory(nout, nodt, segfltot, segs, flow, graph);
        }

        // No phase separation for nodes with more than three segments
        if (nodt > 3)  {
            if (nout == nodt)   {
                for (int i = 0; i < nodt; i++) {
                    graph.hd(segs(i)) = graph.bchd(0);
                }
            }
            else    {
                double flowsum = 0.;
                double hq = 0.;
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

    for (int iseg = 0; iseg < graph.getNseg(); iseg++) {
        uvec idx = find(graph.segname(iseg) == edgeLabels);
        hd(idx).fill(graph.hd(iseg));
    }

}

void Vasculature::woMemory(int &nout, int &nodt, int &segfltot, ivec &segs, vec &flow, Network &eNetwork)    {

    if (nodt == 3)  {
        if (nout == 1)  { // Convergent - use conservation
            eNetwork.hd(segs(0)) = (flow(1)*eNetwork.hd(segs(1)) + flow(2)*eNetwork.hd(segs(2))) / flow(0);
        }
        if (nout == 2) { // Divergent - apply empirical law
            double hdd = (1. - eNetwork.hd(segs(2)))/eNetwork.diam(segs(2));
            double diaquot = sqrt(eNetwork.diam(segs(0))/eNetwork.diam(segs(1)));
            double a = bifpar(2)*(diaquot - 1.)/(diaquot + 1.)*hdd;
            double b = 1. + bifpar(1)*hdd;
            double x0 =  bifpar(0)*hdd;
            double qikdash = (flow(0)/flow(2) - x0)/(1. - 2.*x0);
            if (qikdash <= 0.){
                eNetwork.hd(segs(0)) = 0.;
                eNetwork.hd(segs(1)) = eNetwork.hd(segs(2))*flow(2)/flow(1);
            }
            else if (qikdash >= 1.){
                eNetwork.hd(segs(1)) = 0.;
                eNetwork.hd(segs(0)) = eNetwork.hd(segs(2))*flow(2)/flow(0);
            }
            else {
                double rbcrat = 1./(1. + exp(-a-b*log(qikdash/(1.-qikdash))));
                eNetwork.hd(segs(0)) = rbcrat*eNetwork.hd(segs(2))*flow(2)/flow(0);
                eNetwork.hd(segs(1)) = (1.-rbcrat)*eNetwork.hd(segs(2))*flow(2)/flow(1);
            }
        }
        else if (nout == 3) {
            for (int i = 0; i < 3; i++) {
                eNetwork.hd(segs(i)) = eNetwork.bchd(0);
            }
        }
    }

}

void Vasculature::wMemory(int &nout, int &nodt, int &segfltot, ivec &segs, vec &flow, Network &graph)    {

    if (nodt == 3)  {

        if (nout == 1)  { // Convergent - use conservation
            graph.hd(segs(0)) = (flow(1)*graph.hd(segs(1)) + flow(2)*graph.hd(segs(2))) / flow(0);
        }
        if (nout == 2) { // Divergent - apply empirical law
            double dp = graph.diam(segs(2));
            double length = graph.lseg(segs(2)); // Distance to previous bifurcation - calculated using edge/vertex graph
            double rfn = recovfn(length, dp);
            double x0 = bifpar(0) * (1. - graph.hd(segs(2))) / dp;
            double x0f = x0 * (1. - rfn);
            double x0u = x0 * (1. + rfn);
            double qikdash = (flow(0)/flow(2) - x0f) / (1. - x0u - x0f);
            double qikdash2 = (flow(0) - x0f * flow(2)) / (flow(2) - flow(0) - x0u * flow(2));
            double qikdash3 = (flow(0) - x0f * flow(2)) / (flow(1) - x0u * flow(2));
            if (qikdash <= 0.){
                graph.hd(segs(0)) = 0.;
                graph.hd(segs(1)) = graph.hd(segs(2))*flow(2)/flow(1);
            }
            else if (qikdash >= 1.){
                graph.hd(segs(1)) = 0.;
                graph.hd(segs(0)) = graph.hd(segs(2))*flow(2)/flow(0);
            }
            else {
                double A;
                A = -6.96 * log(graph.diam(segs(0)) / graph.diam(segs(1))) / dp;
                double B = 1. + 6.98 * (1. - graph.hd(segs(2))) / dp;
                double Ashift = 0.5;
                double Af = A + Ashift * rfn;
                double rbfac = exp(Af) * pow(qikdash2, B);
                double rbfac2 = exp(Af) * pow(qikdash3, B);
                graph.hd(segs(0)) = (rbfac / (1. + rbfac)) * (flow(2) / flow(0)) * graph.hd(segs(2));
                graph.hd(segs(1)) = (graph.hd(segs(0)) / rbfac2) * (flow(0) / flow(1));
                if (isnan(graph.hd(segs(1))) && qikdash3 <= 0.)  {
                    graph.hd(segs(1)) = 0.;
                }
                if (graph.hd(segs(0)) >= 1.)  {
                    graph.hd(segs(0)) = consthd;
                }
                if (graph.hd(segs(1)) >= 1.)  {
                    graph.hd(segs(1)) = consthd;
                }
                /*cout<<"3?: "<<graph.hd(segs(1))<<"\t"<<endl;
                cout<<"rbfac2: "<<rbfac2<<endl;
                cout<<"flows: "<<flow(0)<<"\t"<<flow(1)<<"\t"<<flow(2)<<endl;
                cout<<"x0: "<<x0f<<"\t"<<x0u<<endl;
                cout<<"B: "<<B<<endl;
                cout<<"Base: "<<(flow(0) - x0f * flow(2)) / (flow(1) - x0u * flow(2))<<endl;
                cout<<"num.denom: "<<flow(0) - x0f * flow(2)<<"\t"<<flow(1) - x0u * flow(2)<<endl;
                cout<<"HDs: "<<graph.hd(segs(0))<<"\t"<<graph.hd(segs(2))<<endl;*/
            }
        }
        else if (nout == 3) {
            for (int i = 0; i < 3; i++) {
                graph.hd(segs(i)) = bchd(0);
            }
        }
    }

}

// Cell-free layer 90% recovered -> recovfn = 0.1
double Vasculature::recovfn(double len, double dp) {

    double omeg = 4.; // Cell-free layer is ~90% recovered at a distance 10*dp -> omeg = ~4.
    return exp(-len / (omeg * dp));

}

