#include "Vasculature.hpp"
#include "spatGraph.hpp"

using namespace reanimate;

void Vasculature::bloodFlow(bool varviscosity, bool phaseseparation, bool memoryeffects){

    Vasculature::varviscosity = varviscosity;
    Vasculature::phaseseparation = phaseseparation;
    Vasculature::memoryeffects = memoryeffects;

    printText("Blood Flow Module", 3);

    // Detect unknown boundary conditions
    if (any(bctyp == 3))    {unknownBCs = true;}
    setup_flowArrays();

    // If Hd splitting w/ mem. -> computer edge / vertex network to reduce runtime
    if (phaseseparation)    {findDeadends();}

    // Millimetre scaling
    diam *= 1e-3;
    lseg *= 1e-3;
    rseg *= 1e-3;

    rheolParams();
    printText("Simulating blood flow");
    if (unknownBCs)    {
        printText("Unknown boundary conditions detected -> Running flow estimation",2, 0);
        setup_estimationArrays();
        iterateFlowDir();
    }
    else {
        printText("Full boundary conditions detected -> Running full solver",2, 0);
        splitHD(&Network::fullSolver);
    }


    // Micron scaling
    diam *= 1e3;
    lseg *= 1e3;
    rseg *= 1e3;

}


template <typename Call>
void Vasculature::splitHD(Call solver) {

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

    ivec noflowOld = zeros<ivec>(nseg);
    if (phaseseparation)    {printText("Applying phase separations laws",2, 0);}
    for (int iter = 1; iter <= nitmax; iter++)  {

        if(phaseseparation) {printText(to_string(iter)+"/"+to_string(nitmax),1);}

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
                    printText( "Viscosity at segment "+to_string(segname(iseg))+", h'crit = "+to_string(hd(iseg)),4);
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
        /*cout<<q.n_elem<<endl;
        cout<<nseg<<endl;
        uvec idx;
        for (int iseg = 0; iseg < nseg; iseg++) {
            cout<<"sGs:"<<subGraphs(iseg)<<endl;
            cout<<max(subGraphs)<<"\t"<<min(subGraphs)<<endl;
            cout<<deadends(iseg)<<endl;
            cout<<subGraphs(iseg)<<endl;
            if (deadends(iseg) == -3)   {
                idx = find(subGraphs = subGraphs(iseg));
                cout<<idx.n_elem<<endl;
                if (idx.n_elem > 0) {
                    cout << q(idx) << endl;
                    cout << "next!" << endl;
                }
                else {cout << "whoops!" << endl;}
            }
        }
        cout<<"hi"<<endl;*/

        // Calculate segment hematocrits
        if (phaseseparation){
            // Odd behaviour otherwise (between estimation/phase separation)
            //hd(find(noflow == 1)).fill(0.0);
            bool flowSwitch = checkNoflow(noflowOld);
            spatGraph hdGraph;
            if (iter == 1 || flowSwitch)  {
                hdGraph.generate(*this, true);
            }
            putrank(hdGraph);
            dishem(memoryeffects, hdGraph);
            if (any(hd > 1.0))  {
                printText( "Unphysiological haematocrit detected",5);
                hd(find(hd > 1.0)).fill(consthd);
                hd(find(deadends == 1)).fill(0.);
            }
        }


        // Compare hd and q with previous values
        errsegq = 0;
        errseghd = 0;
        qchange = q - qold;
        hdchange = hd - hdold;
        hd = hdold + relax*hdchange;

        maxqerr = abs(qchange).max(errsegq);
        maxhderr = abs(hdchange).max(errseghd);
        if (phaseseparation) {
            printText( "Flow error = "+to_string(maxqerr)+" at segment "+to_string(segname(errsegq))+", h'crit error = "+to_string(maxhderr)+" at segment "+to_string(segname(errseghd)),1,0);
        }
        if(maxqerr < qtol && maxhderr < hdtol)  {
            goto converged;
        }

    }

    if (phaseseparation)   {printText("Non-linear iteration not converged",5);}
    converged:;
    if (phaseseparation) {printText("Iterative h'crit solver converged");}


    // Calculate absolute flow, shear stress and segment pressure
    q(find(noflow == 1)).fill(0.0);
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



    int nnoflow = accu(noflow);
    if (nnoflow > 0) {
        printText("No flow detected in "+to_string(nnoflow)+" segments",5);
        if (phaseseparation) {hd(find(noflow == 1)).fill(0.0);}
    }
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
        printText( "Total reversed flow = "+to_string(nflowreversal)+", ktau/kp = "+to_string(ktau/kp)+", tissue perfusionr = "+to_string(tissperfusion)+" ml/min/100g",1,0);

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

            printText("Final ktau/kp = "+to_string(ktau/kp)+", mean wall shear stress = "+to_string(mean(tau))+" dyn/cm2",1);

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

bool Vasculature::checkNoflow(ivec &noflowOld) {

    noflow(find(q == 0.0)).fill(1);
    noflow(find(deadends == 1)).fill(1);

    bool flowSwitch = false;
    ivec tag = zeros<ivec>(nseg); // Preventing 'noflow' assignment to edge already accessed
    for (int iseg = 0; iseg < getNseg(); iseg++)    {
        if (noflow(iseg) == 1 && tag(iseg) == 0) {
            uvec idx = find(edgeLabels == edgeLabels(iseg));
            noflow(idx).fill(1);
            tag(idx).fill(1);
        }
    }

    for (int iseg = 0; iseg < nseg; iseg++) {
        if (noflow(iseg) != noflowOld(iseg))    {
            flowSwitch = true;
            iseg = nseg;
            printText("Flow switching detected in "+to_string(accu(noflow))+" segments ("+to_string(accu(noflowOld))+" in prior iteration)",5,0);
        }
    }
    noflowOld = noflow;

    return flowSwitch;

}

void Vasculature::findDeadends()    {

    graph.generate(*this);
    graph.findBridgeheads();
    for (int iseg = 0; iseg < nseg; iseg++)    {
        for (int jseg = 0; jseg < graph.getNseg(); jseg++)  {
            if (edgeLabels(iseg) == graph.segname(jseg))    {
                //deadends(iseg) = graph.subGraphs(jseg);
                subGraphs(iseg) = graph.subGraphs(jseg);
                if (subGraphs(iseg) < 0)    {
                    deadends(iseg) = 1;
                }
            }
        }
    }
    hd(find(deadends == 1)).fill(0.);

}
