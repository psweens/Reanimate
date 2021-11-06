//
// Created by sweene01 on 08/03/2021.
//

#include "Vasculature.hpp"

using namespace reanimate;

void Vasculature::analyseVascularFlow()  {

    // PROVIDE ANALYSIS FOR EDGES TOO!

    printText("Analysing vascular flow", 6);
    rheolParams();

    // Blood velocity
    vel = qq * gamma / (M_PI*pow(0.5*diam*1.e-3,2)); // mm/s

    // Reynolds number
    vec Re = zeros<vec>(nseg);
    double visc{};
    for (int iseg = 0; iseg < nseg; iseg++) {
        visc = viscor((diam(iseg)),hd(iseg))*xi;
        Re(iseg) = bloodDensity * vel(iseg) * lseg(iseg)*1.e-3 / visc;
    }

    // Vascular segment transit time
    vec vstt = zeros<vec>(nseg);
    vstt = 60*(M_PI*pow(1e-3*diam(find(qq > 0.0)),2) % lseg(find(qq > 0.0))*1e-3) / (4.*gamma*qq(find(qq > 0.0)));

    printText("Blood pressure = "+to_string(mean(segpress))+" ± "+to_string(stddev(segpress))+" mmHg",1, 0);
    if (max(vesstyp) == 3)  {
        printText("Arteriole flow = "+to_string(mean(segpress(find(vesstyp == 1))))+" ± "+to_string(stddev(segpress(find(vesstyp == 1))))+" mmHg",1, 0);
        printText("Capillary flow = "+to_string(mean(segpress(find(vesstyp == 2))))+" ± "+to_string(stddev(segpress(find(vesstyp == 2))))+" mmHg",1, 0);
        printText("Venule flow = "+to_string(mean(segpress(find(vesstyp == 3))))+" ± "+to_string(stddev(segpress(find(vesstyp == 3))))+" mmHg",1, 0);
    }
    printText("Blood flow = "+to_string(mean(qq))+" ± "+to_string(stddev(qq))+" nl/min",1, 0);
    if (max(vesstyp) == 3)  {
        printText("Arteriole flow = "+to_string(mean(qq(find(vesstyp == 1))))+" ± "+to_string(stddev(qq(find(vesstyp == 1))))+" nl/min",1, 0);
        printText("Capillary flow = "+to_string(mean(qq(find(vesstyp == 2))))+" ± "+to_string(stddev(qq(find(vesstyp == 2))))+" nl/min",1, 0);
        printText("Venule flow = "+to_string(mean(qq(find(vesstyp == 3))))+" ± "+to_string(stddev(qq(find(vesstyp == 3))))+" nl/min",1, 0);
    }
    printText("Blood velocity = "+to_string(mean(vel))+" ± "+to_string(stddev(vel))+" mm/s",1, 0);
    if (max(vesstyp) == 3)  {
        printText("Arteriole velocity = "+to_string(mean(vel(find(vesstyp == 1))))+" ± "+to_string(stddev(vel(find(vesstyp == 1))))+" mm/s",1, 0);
        printText("Capillary velocity = "+to_string(mean(vel(find(vesstyp == 2))))+" ± "+to_string(stddev(vel(find(vesstyp == 2))))+" mm/s",1, 0);
        printText("Venule velocity = "+to_string(mean(vel(find(vesstyp == 3))))+" ± "+to_string(stddev(vel(find(vesstyp == 3))))+" mm/s",1, 0);
    }
    printText("Vessel Wall Shear Stress = "+to_string(mean(tau))+" ± "+to_string(stddev(tau))+" dyn/cm2",1, 0);
    if (max(vesstyp) == 3)  {
        printText("Arteriole WSS = "+to_string(mean(tau(find(vesstyp == 1))))+" ± "+to_string(stddev(tau(find(vesstyp == 1))))+" dyn/cm2",1, 0);
        printText("Capillary WSS = "+to_string(mean(tau(find(vesstyp == 2))))+" ± "+to_string(stddev(tau(find(vesstyp == 2))))+" dyn/cm2",1, 0);
        printText("Venule WSS = "+to_string(mean(tau(find(vesstyp == 3))))+" ± "+to_string(stddev(tau(find(vesstyp == 3))))+" dyn/cm2",1, 0);
    }
    printText("Reynolds number = "+to_string(mean(Re))+" ± "+to_string(stddev(Re)),1, 0);
    printText("Vascular segment transit time = "+to_string(mean(vstt))+" ± "+to_string(stddev(vstt))+" s",1, 0);

}
