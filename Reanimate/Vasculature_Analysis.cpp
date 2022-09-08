//
// Created by sweene01 on 08/03/2021.
//

#include "Vasculature.hpp"
#include <math.h>

using namespace reanimate;

void Vasculature::analyseVascularFlow()  {

    printText("Analysing vascular flow", 6);
    rheolParams();

    // Blood velocity
    if (cuboidVessels)  {vel = qq * gamma / (pow(1.e-3,2) * width % height);} // Average velocity
    else {vel = qq * gamma / (M_PI*pow(0.5*diam*1e-3,2));} // mm/s

    // Reynolds number
    vec Re = zeros<vec>(nseg);
    double visc{};
    for (int iseg = 0; iseg < nseg; iseg++) {
        if (cuboidVessels) {visc = constvisc * xi;}
        else {visc = viscor((diam(iseg)),hd(iseg))*xi;}
        Re(iseg) = bloodDensity * vel(iseg) * lseg(iseg)*1e-3 / visc;
    }

    // Vascular segment transit time
    vec vstt = zeros<vec>(nseg);
    uvec idx = find(qq > 0.);
    if (cuboidVessels)  {vstt = 60*(pow(1.e-3,2) * width(idx) % height(idx) % lseg(idx)*1e-3) / (4.*gamma*qq((idx)));}
    else {vstt = 60*(M_PI*pow(1e-3*diam((idx)),2) % lseg((idx))*1e-3) / (4.*gamma*qq((idx)));}

    // Get vessel class indices
    idx = find(vesstyp == 1);
    uvec jdx = find(vesstyp == 2);
    uvec kdx = find(vesstyp == 3);
    
    printText("Blood pressure = "+to_string(mean(segpress))+" ± "+to_string(stddev(segpress))+" mmHg", 1, 0);
    if (max(vesstyp) == 3)  {
        printText("Arteriole pressure = "+to_string(mean(segpress(idx)))+" ± "+to_string(stddev(segpress(idx)))+" mmHg", 1, 0);
        printText("Capillary pressure = "+to_string(mean(segpress(jdx)))+" ± "+to_string(stddev(segpress(jdx)))+" mmHg", 1, 0);
        printText("Venule pressure = "+to_string(mean(segpress(kdx)))+" ± "+to_string(stddev(segpress(kdx)))+" mmHg",1, 0);
    }
    printText("Blood flow = "+to_string(mean(qq))+" ± "+to_string(stddev(qq))+" nl/min",1, 0);
    if (max(vesstyp) == 3)  {
        printText("Arteriole flow = "+to_string(mean(qq(idx)))+" ± "+to_string(stddev(qq(idx)))+" nl/min",1, 0);
        printText("Capillary flow = "+to_string(mean(qq(jdx)))+" ± "+to_string(stddev(qq(jdx)))+" nl/min",1, 0);
        printText("Venule flow = "+to_string(mean(qq(kdx)))+" ± "+to_string(stddev(qq(kdx)))+" nl/min",1, 0);
    }
    printText("Blood velocity = "+to_string(mean(vel))+" ± "+to_string(stddev(vel))+" mm/s",1, 0);
    if (max(vesstyp) == 3)  {
        printText("Arteriole velocity = "+to_string(mean(vel(idx)))+" ± "+to_string(stddev(vel(idx)))+" mm/s",1, 0);
        printText("Capillary velocity = "+to_string(mean(vel(jdx)))+" ± "+to_string(stddev(vel(jdx)))+" mm/s",1, 0);
        printText("Venule velocity = "+to_string(mean(vel(kdx)))+" ± "+to_string(stddev(vel(kdx)))+" mm/s",1, 0);
    }
    printText("Vessel Wall Shear Stress = "+to_string(mean(tau))+" ± "+to_string(stddev(tau))+" dyn/cm2",1, 0);
    if (max(vesstyp) == 3)  {
        printText("Arteriole WSS = "+to_string(mean(tau(idx)))+" ± "+to_string(stddev(tau(idx)))+" dyn/cm2",1, 0);
        printText("Capillary WSS = "+to_string(mean(tau(jdx)))+" ± "+to_string(stddev(tau(jdx)))+" dyn/cm2",1, 0);
        printText("Venule WSS = "+to_string(mean(tau(kdx)))+" ± "+to_string(stddev(tau(kdx)))+" dyn/cm2",1, 0);
    }
    printText("Reynolds number = "+to_string(mean(Re))+" ± "+to_string(stddev(Re)),1, 0);
    printText("Vascular segment transit time = "+to_string(mean(vstt))+" ± "+to_string(stddev(vstt))+" s",1, 0);

}

void Network::cuboidWSS()    {

    vec int_tauY = zeros<vec>(nseg);
    vec int_tauZ = zeros<vec>(nseg);
    double a{},b{},sumY{},sumZ{},coef,bn{},grad{};
    width *= 1.e-3;
    height *= 1.e-3;
    double visc = constvisc * xi;
    for (int iseg = 0; iseg < nseg; iseg++) {
        sumY = 0.;
        sumZ = 0.;
        a = width(iseg);
        b = height(iseg);
        grad = nodpress(ista(iseg)) - nodpress(iend(iseg));
        int_tauY(iseg) = (grad * b) / (2. * visc);
        coef = (4 * grad * pow(b,2)) / (visc * pow(M_PI,3));
        for (int n = 0; n < 5; n++) {
            bn = (2.*n - 1.) * M_PI / a;
            sumY += (2. / pow(2.*n - 1,3)) * tanh(0.5*bn*a);
            sumZ += (1. / pow(2*n - 1,3)) * (cosh(bn*a) / sinh(bn*a)) * (1 - cos(bn*a));
        }
        int_tauY(iseg) -= coef * sumZ;
        int_tauZ(iseg) = coef * sumZ;
    }
    width *= 1.e3;
    height *= 1.e3;

    // Average wall shear stress
    tau = (1. / (width + height)) % (int_tauY + int_tauZ) * (gamma/beta);

}
