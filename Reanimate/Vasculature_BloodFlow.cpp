//
//  Network_BloodFlow.cpp
//  TEEPEE
//
//  Created by Paul Sweeney on 20/06/2020.
//  Copyright © 2020 Paul Sweeney. All rights reserved.
//

#include "Vasculature.hpp"

using namespace reanimate;

void Vasculature::bloodFlow(bool varviscosity, bool phaseseparation, bool memoryeffects){

    Vasculature::varviscosity = varviscosity;
    Vasculature::phaseseparation = phaseseparation;
    Vasculature::memoryeffects = memoryeffects;

    cout<<"\n\nSimulating blood flow ..."<<endl;

    // If Hd splitting w/ mem. -> compute edge lengths (microns)
    if (memoryeffects && phaseseparation && nedge == 0)    {edgeNetwork();}

    // Millimetre scaling
    diam *= 1e-3;
    lseg *= 1e-3;
    rseg *= 1e-3;
    
    // Detect unknown boundary conditions
    if (any(bctyp == 3))    {unknownBCs = true;}

    setup_flowArrays();
    rheolParams();
    if (unknownBCs)    {
        cout<<"Unknown boundary conditions detected -> Running flow estimation ..."<<endl;
        setup_estimationArrays();
        iterateFlowDir();
    }
    else {splitHD(&Network::fullSolver);}


    // Micron scaling
    diam *= 1e3;
    lseg *= 1e3;
    rseg *= 1e3;

}


template <typename CallAnother>
void Vasculature::splitHD(CallAnother solver) {

    // Nonlinear iterations.  If convergence is slow, hematocrit is increasingly underrelaxed.
    // This eventually forces convergence even when system is unstable due to rheological
    // effects (see work of R.L. Carr et al.).
    uword errsegq{},errseghd{};
    int nitmax = 1e2;
    double relax = 1., maxqerr{}, maxhderr{};
    vec qchange,hdchange;

    if (!phaseseparation)   {
        nitmax = 1;
    }
    for (int iter = 1; iter <= nitmax; iter++)  {

        printf("%i/%i: ",iter, nitmax);

        if (iter % 5 == 0)  {
            relax *= 0.8;
        }

        qold = q;
        hdold = hd;
        double visc{};
        double tdiam{};
        for (int iseg = 0; iseg < nseg; iseg++) {
            tdiam = diam(iseg);
            if (varviscosity) {
                visc = viscor((tdiam)*1e3,hd(iseg))*xi;
                if (isnan(visc))    {
                    printf("*** Viscosity error at segment %lli, h'crit = %f ***\n",segname(iseg),hd(iseg));
                }
            }
            else {
                visc = constvisc*xi;
            }
            conductance(iseg) = M_PI*pow(tdiam,4)/(128*visc*lseg(iseg));
            c(iseg) = 4*visc/(M_PI*pow(tdiam*0.5,3));
        }


        // Flow solver
        (this->*solver)();


        // Calculate segment hematocrits
        if (phaseseparation){
            if (unknownBCs) {
                //zeroflow(find(q == 0.0)).fill(1);
                //q(find(q == 0.0)).fill(1e-12);
            } // Odd behaviour otherwise (between estimation/phase separation)
            putrank();
            dishem(memoryeffects);
        }


        // Assuming zero flow --> no Hd .. or at least zero --> RBCs not passing through lumen
        //hd(find(abs(q) <= 1e-9)).fill(0.0);
        //hd(find(hd > 1.0)).fill(consthd);


        // Compare hd and q with previous values
        errsegq = 0;
        errseghd = 0;
        qchange = q - qold;
        hdchange = hd - hdold;
        hd = hdold + relax*hdchange;

        maxqerr = abs(qchange).max(errsegq);
        maxhderr = abs(hdchange).max(errseghd);
        if (phaseseparation) {
            printf("Flow error = %f, h'crit error = %f\n", iter, maxqerr, maxhderr);
        }
        if(maxqerr < qtol && maxhderr < hdtol)  {
            goto converged;
        }

    }
    errsegq = (int) errsegq;
    errseghd = (int) errseghd;

    if (phaseseparation)   {
        printf("*** WARNING: Non-linear iteration not converged ***\n");
        printf("Flow error = %f at segment %lli, Hd error = %f at segment %lli\n",maxqerr,segname(errsegq),maxhderr,segname(errseghd));
    }
    converged:;


    q(find(zeroflow == 1)).fill(0.0);
    // Calculate absolute flow, shear stress and segment pressure
    qq = abs(q);
    if (!unknownBCs)    {
        for (int iseg = 0; iseg < nseg; iseg++) {
            segpress(iseg) = (nodpress(ista(iseg)) + nodpress(iend(iseg)))/2.;
        }

        for (int inodbc = 0; inodbc < nnodbc; inodbc++) {
            BCpress(inodbc) = nodpress(bcnod(inodbc));
        }
        tau = (c % qq)*(gamma/beta);
    }
    else    {
        tau = (c % qq);
    }



    int nzeroflow = accu(zeroflow);
    if (nzeroflow > 0) {printf("*** WARNING: Zero flow detected in %i segments ***\n",nzeroflow);}
}

void Vasculature::iterateFlowDir()   {

    bool stableflow = false;
    while (!stableflow)    {

        splitHD(&Network::estimationSolver);

        // Calculate flow sign and assign tau0 accordingly
        flowsign = sign(q);
        double K = mean(tau % tau) / pow(mean(tau),2);
        tau0 = targStress * K * flowsign;


        // Count flow reversal
        int nflowreversal = (int) accu(abs(flowsign(find(flowsign != oldFlowsign))));
        cout <<"Total reversed flow = "<<nflowreversal<<",\t"<<"ktau/kp = "<<ktau/kp<<",\t"<<"Tissue perfusion = "<<tissperfusion<<" ml/min/100g\n"<<endl;


        // Condition for final flow solution
        if (ktau / kp > 1 && nflowreversal == 0 && nodpress.min() > 0.) {

            stableflow = true;
            tissperfusion = 0.;
            inflow = 0.;
            ktau  = oldktau;
            tau = abs(oldTau) / beta;
            q = oldq / gamma;
            qq = abs(qq);
            nodpress = oldNodpress / alpha;
            hd = oldHd;

            cout<<"Final ktau/kp = "<<ktau/kp<<"\tMean wall shear stress = "<<mean(tau)<<" dyn/cm2"<<endl;

            for (int iseg = 0; iseg < nseg; iseg++) {
                segpress(iseg) = (nodpress(ista(iseg)) + nodpress(iend(iseg)))/2.;
            }

            for (int inodbc = 0; inodbc < nnodbc; inodbc++) {
                BCpress(inodbc) = nodpress(bcnod(inodbc));
            }

        }
        else {

            // Increase ktau
            oldktau = ktau;
            ktau *= 2;

            oldFlowsign = flowsign;
            oldTau = tau;
            oldHd = hd;
            oldq = q;
            oldNodpress = nodpress;

        }

    }

}
