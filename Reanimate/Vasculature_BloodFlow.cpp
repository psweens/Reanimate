#include "Vasculature.hpp"
#include "spatGraph.hpp"

using namespace reanimate;

void Vasculature::bloodFlow(bool varViscosity, bool phaseSeparation, bool memoryEffects, bool updateBCHD, bool skipAnalysis){

    varviscosity = varViscosity;
    phaseseparation = phaseSeparation;
    memoryeffects = memoryEffects;
    updateBoundaryHD = updateBCHD;

    printText("Blood Flow Module", 3);
    setup_flowArrays(false);

    // Detect unknown boundary conditions
    if (any(bctyp == 3))    {unknownBCs = true;}

    printText("Cloning vascular network");
    Vasculature networkCopy = *this;

    if (loadDeadEnds)   {
        printText("Reading dead ends");
        deadEnds.load(string(buildPath+"Network_DeadEnds.txt"), raw_ascii);
    }
    else {
        deadEnds = networkCopy.findDeadends();
        deadEnds.save(string(buildPath+"Network_DeadEnds.txt"), raw_ascii);
    }
    if (accu(deadEnds) > 0) {
        printText("Removing dead ends",2,0);
        networkCopy.subNetwork(deadEnds, false, false);
        mat extraD = zeros<mat>(networkCopy.getNseg(), 1);
    }


    spatGraph hdGraph;
    if (phaseSeparation) {
        hdGraph.generate(networkCopy, true, graphOverride); // Diameter / length dimensions are in microns (taken from edge data)
        npoint = networkCopy.npoint;
        if (any(hdGraph.nodtyp == 2) && memoryeffects)   {
            printText( "Type 2 vertex detected. Amending ...",1,0);
            hdGraph.linkEdges();
        }
    }

    // Millimetre scaling for copied network
    networkCopy.diam *= 1e-3;
    networkCopy.lseg *= 1e-3;
    networkCopy.rseg *= 1e-3;

    // Assign rheological parameters
    networkCopy.rheolParams();

    printText("Simulating blood flow");
    networkCopy.setup_flowArrays();
    networkCopy.nitmax = 100;
    if (unknownBCs)    {
        printText("Unknown boundary conditions detected -> Running flow estimation",2);
        networkCopy.setup_estimationArrays();
        networkCopy.iterateFlowDir(hdGraph);
    }
    else {
        printText("Full boundary conditions detected -> Running full solver",2);
        networkCopy.splitHD(&Network::fullSolver, hdGraph);
    }

    // Map solved flow to original network
    mapFlow(networkCopy);
    computeBoundaryFlow();

    // Analyse flow
    analyseVascularFlow();
    printNetwork("SolvedBloodFlow.txt");

}


template <typename Call>
void Vasculature::splitHD(Call solver, spatGraph &hdGraph) {

    // Nonlinear iterations.  If convergence is slow, hematocrit is increasingly underrelaxed.
    // This eventually forces convergence even when system is unstable due to rheological
    // effects (see work of R.L. Carr et al.).
    double relax = 1., maxqerr{}, maxhderr{};

    if (!phaseseparation)   {nitmax = 1e2;}

    bool converged = false;
    ivec noflowOld = zeros<ivec>(nseg);
    if (phaseseparation)    {printText("Applying phase separations law",2);}
    for (int iter = 1; iter <= nitmax; iter++)  {

        // Update relaxation
        if (iter % 5 == 0)  {relax *= 0.8;}
        track = iter;


        // Calculate conductance
        if(phaseseparation) {
            printText(to_string(iter)+"/"+to_string(nitmax)+": ",1, -1);
            computeConductance();
        }
        else {computeConductance();}

        // Flow solver
        (this->*solver)();

        if (any(q == 0.0))  {printText( "No flow detected",5);
//            vec temp = zeros<vec>(nseg);
//            temp(find(q == 0.0)).fill(1.);
//            mat extraD = zeros<mat>(nseg, 1);
//            extraD.col(0) = temp;
//            printAmira("amiraNoFlow.am", extraD);
//            cout<<accu(temp)<<endl;
            /*for (int iseg = 0; iseg < nseg; iseg++) {
                if (q(iseg) == 0.0) {
                    if (nodtyp(ista(iseg)) == 2)    {
                        cout<<"sta: "<<nodpress(ista(iseg))<<"\t"<<nodpress(iend(iseg))<<endl;
                        for (int jseg = 0; jseg < nseg; jseg++) {
                            if (ista(iseg) == ista(jseg) || ista(iseg) == iend(jseg)) {
                                if (iseg != jseg) { cout << "\t" << q(jseg) << endl; }
                            }
                        }
                    }
                    else if (nodtyp(iend(iseg)) == 2)    {
                        cout<<"end: "<<nodpress(ista(iseg))<<"\t"<<nodpress(iend(iseg))<<endl;
                        for (int jseg = 0; jseg < nseg; jseg++) {
                            if (iend(iseg) == ista(jseg) || iend(iseg) == iend(jseg)) {
                                if (iseg != jseg) { cout << "\t" << q(jseg) << endl; }
                            }
                        }
                    }
                }
            }*/

        }
        if (unknownBCs) {
            // To allow tau0 to update flow directions rather than magnitude
            flowsign = sign(q);
            tau0 = abs(tau0) % flowsign;
        }


        // Calculate segment hematocrits
        if (phaseseparation){
            putrank(hdGraph);
            dishem(memoryeffects, hdGraph);

            if (any(hd > 1.0))  {
                uword segerr{};
                double maxhderr= hd.max(segerr);
                printText( "Unphysiological haematocrit detected, h'crit = "+to_string(maxhderr)+" at segment "+to_string(segerr),5);
            }

            // Compare hd and q with previous values
            vec errs = computeFlowError(relax);
            maxqerr = errs(0);
            maxhderr = errs(1);
            //hd(find(hd > 1.0)).fill(consthd);
            for (int iseg = 0; iseg < nseg; iseg++) {
                if (hd(iseg) > 1.0 || isnan(hd(iseg)) || isinf(hd(iseg))) {hd(iseg) = consthd;}
            }

            qold = q;
            hdold = hd;

        }
        if((maxqerr < qtol && maxhderr < hdtol))  {
            iter = nitmax;
            converged = true;
        }

    }
    if (phaseseparation && !converged)   {printText("Non-linear iteration not converged",5);}

}

void Vasculature::iterateFlowDir(spatGraph &hdGraph)   {

    printText("Starting flow estimation",2);

    int iteration = 1;
    bool stableflow = false;
    bool relaxFlow = false;
    if (phaseseparation)    {
        relaxFlow = true;
        phaseseparation = false;
    }
    while (!stableflow)    {

        if (ktau / kp > 1 && nodpress.min() > 0. && relaxFlow) {phaseseparation = true;}

        printText("Iteration "+to_string(iteration), 6);
        splitHD(&Network::estimationSolver, hdGraph);

        // Calculate flow sign and assign tau0 accordingly
        flowsign = sign(q);
        double newK = mean(tau % tau) / pow(mean(tau),2);
        tau0 = targStress * newK * flowsign;

        // Rough estimation of tissue perfusion - IMPROVE - below too
        computeBoundaryFlow();
        tissperfusion = accu(find(BCflow > 0.)) / (alx * aly * alz * tissDensity) * 1e6; // ml/min/100g

        // Count flow reversal
        int nflowreversal = (int) accu(abs(flowsign(find(flowsign != oldFlowsign))));
        printText( "Total reversed flow = "+to_string(nflowreversal)+", ktau/kp = "+to_string(ktau/kp)+", tissue perfusion = "+to_string(tissperfusion)+" ml/min/100g",1,0);
        cout<<"Min pressure: "<<nodpress.min() / alpha<<endl;
//        if (ktau / kp > 1 && nflowreversal < 50 && nodpress.min() > 0.0) {nitmax = 200;}
//        else {nitmax = 100;}

        // Condition for final flow solution
        if (ktau / kp > 1 && nflowreversal == 0 && nodpress.min() > 0.0) {

            stableflow = true;
            inflow = 0.;
            ktau  = oldktau;
            tau = abs(oldTau) / beta;
            q = oldq / gamma;
            qq = abs(qq);
            nodpress = oldNodpress / alpha;
            hd = oldHd;

            printText("Final ktau/kp = "+to_string(ktau/kp)+", mean wall shear stress = "+to_string(mean(tau))+" dyn/cm2",1);

            computeSegpress();
            computeBoundaryFlow();
            tissperfusion = accu(find(BCflow > 0.)) / (alx * aly * alz * tissDensity) * 1e6; // ml/min/100g

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

            if (updateBoundaryHD) {relaxBoundaryHD();}

        }
        iteration += 1;
    }

}

void Vasculature::computeConductance()   {

    double visc{};
    double tdiam{};
    conductance = zeros<vec>(nseg);
    c = zeros<vec>(nseg);
    for (int iseg = 0; iseg < nseg; iseg++) {
        tdiam = diam(iseg);
        if (varviscosity) {
            visc = viscor(tdiam*1e3,hd(iseg))*xi;
            if (isnan(visc))    {
                printText( "Viscosity at segment "+to_string(segname(iseg))+", h'crit = "+to_string(hd(iseg)),4);
            }
        }
        else {visc = constvisc*xi;}
        conductance(iseg) = M_PI*pow(tdiam,4)/(128.*visc*lseg(iseg));
        c(iseg) = 4*visc/(M_PI*pow(tdiam*0.5,3));
    }

}

void Vasculature::computeSegpress() {

    segpress = zeros<vec>(nseg);
    for (int iseg = 0; iseg < nseg; iseg++) {segpress(iseg) = (nodpress(ista(iseg)) + nodpress(iend(iseg)))/2.;}

}

vec Vasculature::computeFlowError(double &relax)    {

    uword errsegq{},errseghd{};
    vec qchange = (q - qold) / gamma;
    vec hdchange = hd - hdold;
    hd = hdold + relax*hdchange;
    if (any(hd < 0))    {
        printText("Negative h'crit",5);
        hd(find(hd < 0)).fill(0.0);
    }
    double maxqerr = abs(qchange).max(errsegq);
    double maxhderr = abs(hdchange).max(errseghd);
    vec errs = zeros<vec>(2);
    errs(0) = maxqerr;
    errs(1) = maxhderr;

    printText( "Flow error = "+to_string(maxqerr)+" at segment "+to_string(segname(errsegq))+", h'crit error = "+to_string(maxhderr)+" at segment "+to_string(segname(errseghd)),1,0);

    return errs;

}


void Vasculature::relaxBoundaryHD() {

    for (int inodbc = 0; inodbc < nnodbc; inodbc++) {
        if (bctyp(inodbc) == 3) {
            for (int iseg = 0; iseg < nseg; iseg++) {
                if (bcnod(inodbc) == ista(iseg))    {
                    if (nodpress(ista(iseg)) < nodpress(iend(iseg)))    {
                        bchd(inodbc) = hd(iseg);
                        iseg = nseg;
                    }
                }
                else if (bcnod(inodbc) == iend(iseg))   {
                    if (nodpress(iend(iseg)) < nodpress(ista(iseg)))    {
                        bchd(inodbc) = hd(iseg);
                        iseg = nseg;
                    }
                }
            }
        }
    }

}


void Vasculature::mapFlow(Vasculature &Network) {

    int found{};
    ivec nodFlag = -3*ones<ivec>(nnod);
    for (int inod = 0; inod < nnod; inod++) {
        found = 0;
        for (int jnod = 0; jnod < Network.getNnod(); jnod++) {
            if (nodname(inod) == Network.nodname(jnod)) {
                nodpress(inod) = Network.nodpress(jnod);
                if (Network.nodtyp(jnod) != nodtyp(inod))   {nodFlag(inod) = -2;}
                jnod = Network.getNnod();
                found = 1;
            }
        }
        if (found == 0) {nodFlag(inod) = -1;}
    }

    for (int inod = 0; inod < nnod; inod++) {
        if (nodFlag(inod) == -2)    {
            dfsBasic(inod,nodpress(inod),nodFlag);
        }
    }
    for (int inod = 0; inod < nnod; inod++) {
        if (nodFlag(inod) >= 0) {nodpress(inod) = nodFlag(inod);}
    }

    computeSegpress();

    q.zeros();
    tau.zeros();
    for (int iseg = 0; iseg < Network.getNseg(); iseg++)    {
        for (int jseg = 0; jseg < nseg; jseg++) {
            if (Network.segname(iseg) == segname(jseg)) {
                q(jseg) = Network.q(iseg);
                hd(jseg) = Network.hd(iseg);
                tau(jseg) = Network.c(iseg) * abs(Network.q(iseg)) * (gamma/beta);
                jseg = nseg;
            }
        }
    }

    hd(find(deadEnds == 1)).fill(0.0);
    q(find(deadEnds == 1)).fill(0.0);

    // Calculate absolute flow, shear stress and segment pressure
    q(find(noflow == 1)).fill(0.0);
    hd(find(noflow == 1)).fill(0.0);
    qq = abs(q);
    if (!unknownBCs)    {
        for (int inodbc = 0; inodbc < nnodbc; inodbc++) {BCpress(inodbc) = nodpress(bcnod(inodbc));}
        bcprfl = BCpress;
        bctyp.fill(0);
    }

    int nnoflow = accu(noflow);
    if (nnoflow > 0) {
        printText("No flow detected in "+to_string(nnoflow)+" segments",5);
        if (phaseseparation) {hd(find(noflow == 1)).fill(0.0);}
    }

}
