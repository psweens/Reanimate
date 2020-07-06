//
// Created by sweene01 on 21/06/2020.
//

#include "Vasculature.hpp"

using namespace reanimate;

void Vasculature::rheolParams()  {

    cout<<"Assigning rheology parameters..."<<endl;

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

void Vasculature::dishem(bool &memoryeffects)    {

    diam *= 1e3;

    int segfltot = nseg;
    ivec segs = zeros<ivec>(nodsegm);
    vec flow = zeros<vec>(nodsegm);

    // Assign Hd to boundary nodes
    int isegk{}; // Counts inflowing nodes
    for (int inodbc = 0; inodbc < nnodbc; inodbc++) {
        if (nodout(bcnod(inodbc)) >= 1)  {
            hd(nodseg(0,bcnod(inodbc))) = bchd(inodbc);
            isegk += 1;
        }
    }


    // Cycle through interior inflowing nodes
    for (int in = 0; in < nnodfl; in++) {

        // If node is an interior node:
        // Store associated segment indices in increasing rank order
        // Store absolute flow value of segment in increasing rank order
        int nodt = (int) nodtyp(nodrank(in));
        int nout = (int) nodout(nodrank(in));
        if (nodt >= 2)  {
            for (int i = 0; i < nodt; i++)  {
                segs(i) = nodseg(i,nodrank(in));
                flow(i) = abs(q(segs(i)));
            }
        }

        // 2-segment nodes: nout = 2 could happen in iterations
        if (nodt == 2)  {
            if (nout == 1)  {
                hd(segs(0)) = hd(segs(1));
            }
            else if (nout == 2)  {
                hd(segs(0)) = bchd(0);
                hd(segs(1)) = bchd(0);
            }
        }

        // 3-segment nodes: convergent, divergent
        if (!memoryeffects)  {
            woMemory(nout, nodt, segfltot, segs, flow);
        }
        else{
            // With memory effects - daughter with: high flow (favourable, f), low flow (unfavourable, u)
            wMemory(nout, nodt, segfltot, segs, flow);
        }

        // No phase separation for nodes with more than three segments
        if (nodt > 3)  {
            if (nout == nodt)   {
                for (int i = 0; i < nodt; i++) {
                    hd(segs(i)) = bchd(0);
                }
            }
            else    {
                double flowsum = 0.;
                double hq = 0.;
                for (int i = nout; i < nodt; i++)   {
                    flowsum += abs(q(segs(i)));
                    hq += hd(segs(i))*abs(q(segs(i)));
                }
                if (nout >= 1)  {
                    for (int i = 0; i < nout; i++)  {
                        hd(segs(i)) = hq/flowsum;

                    }
                }
            }
        }
        if (nodt >= 2)  {
            isegk += nout;
        }
    }
    if(isegk != segfltot) {
        printf("*** ERROR: Hd computed in %i of %i segments processed ***\n", isegk, segfltot);
    }

    diam *= 1e-3;
}

void Vasculature::woMemory(int &nout, int &nodt, int &segfltot, ivec &segs, vec &flow)    {

    if (nodt == 3)  {
        if (nout == 1)  { // Convergent - use conservation
            hd(segs(0)) = (flow(1)*hd(segs(1)) + flow(2)*hd(segs(2))) / flow(0);
        }
        if (nout == 2) { // Divergent - apply empirical law
            double hdd = (1. - hd(segs(2)))/diam(segs(2));
            double diaquot = sqrt(diam(segs(0))/diam(segs(1)));
            double a = bifpar(2)*(diaquot - 1.)/(diaquot + 1.)*hdd;
            double b = 1. + bifpar(1)*hdd;
            double x0 =  bifpar(0)*hdd;
            double qikdash = (flow(0)/flow(2) - x0)/(1. - 2.*x0);
            if (qikdash <= 0.){
                hd(segs(0)) = 0.;
                hd(segs(1)) = hd(segs(2))*flow(2)/flow(1);
            }
            else if (qikdash >= 1.){
                hd(segs(1)) = 0.;
                hd(segs(0)) = hd(segs(2))*flow(2)/flow(0);
            }
            else {
                double rbcrat = 1./(1. + exp(-a-b*log(qikdash/(1.-qikdash))));
                hd(segs(0)) = rbcrat*hd(segs(2))*flow(2)/flow(0);
                hd(segs(1)) = (1.-rbcrat)*hd(segs(2))*flow(2)/flow(1);
            }
        }
        else if (nout == 3) {
            for (int i = 0; i < 3; i++) {
                hd(segs(i)) = bchd(0);
            }
        }
    }

}

void Vasculature::wMemory(int &nout, int &nodt, int &segfltot, ivec &segs, vec &flow)    {

    if (nodt == 3)  {
        if (nout == 1)  { // Convergent - use conservation
            hd(segs(0)) = (flow(1)*hd(segs(1)) + flow(2)*hd(segs(2))) / flow(0);
        }
        if (nout == 2) { // Divergent - apply empirical law
            double dp = diam(segs(2));
            double length = elseg(segs(2)); // Distance to previous bifurcation - calculated using tort_l fn
            double rfn = recovfn(length, dp);
            double x0 = (1. - hd(segs(2))) / dp;
            double x0f = x0 * (1. - rfn);
            double x0u = x0 * (1. + rfn);
            double qikdash = (flow(0)/flow(2) - x0f) / (1. - x0u - x0f);
            if (qikdash <= 0.){
                hd(segs(0)) = 0.;
                hd(segs(1)) = hd(segs(2))*flow(2)/flow(1);
            }
            else if (qikdash >= 1.){
                hd(segs(1)) = 0.;
                hd(segs(0)) = hd(segs(2))*flow(2)/flow(0);
            }
            else {
                double A = -6.96 * log(diam(segs(0)) / diam(segs(1))) / dp; // Typo in preprint - may not be correct
                double B = 1. + 6.98 * (1. - hd(segs(2))) / dp;
                double Ashift = 0.5;
                double Af = A + Ashift * rfn;
                double rbfac = exp(Af) * pow((flow(0) - x0f * flow(2)) / (flow(2) - flow(0) - x0u * flow(2)), B);
                double rbfac2 = exp(Af) * pow((flow(0) - x0f * flow(2)) / (flow(1) - x0u * flow(2)), B);
                //cout<<hd(segs(0))<<"\t"<<hd(segs(1))<<"\t"<<hd(segs(2))<<endl;
                hd(segs(0)) = (rbfac / (1. + rbfac)) * (flow(2) / flow(0)) * hd(segs(2));
                if (hd(segs(0)) >= 1.)  {
                    hd(segs(0)) = consthd;
                }
                hd(segs(1)) = (hd(segs(0)) / rbfac2) * (flow(0) / flow(1));
                if (hd(segs(1)) >= 1.)  {
                    hd(segs(1)) = consthd;
                }
            }
        }
        else if (nout == 3) {
            for (int i = 0; i < 3; i++) {
                hd(segs(i)) = bchd(0);
            }
        }
    }

}

// Cell-free layer 90% recovered -> recovfn = 0.1
double Vasculature::recovfn(double len, double dp) {

    double omeg = 4.; // Cell-free layer is ~90% recovered at a distance 10*dp -> omeg = ~4.
    return exp(-len / (omeg * dp));

}
