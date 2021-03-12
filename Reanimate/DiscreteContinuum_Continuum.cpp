//
// Created by Paul Sweeney on 27/02/2021.
//

#include "DiscreteContinuum.hpp"

using namespace reanimate;

void DiscreteContinuum::computeContinuum()  {

    printText("Setting up continuum component");

    printText("Scaling coordinates based on anisotropic micro-cell conductivity");
    if (cell.aniScaleY > 0.)   {discreteNet.cnode.row(1) /= cell.aniScaleY;}
    if (cell.aniScaleZ > 0.)   {discreteNet.cnode.row(2) /= cell.aniScaleZ;}

    printText("Maximum base pressure (mmHg) = "+to_string(Pa_Pv.max()),0);
    printText("Minimum base pressure (mmHg) = "+to_string(Pa_Pv.min()),0,0);
    printText("Mean capillary pressure (mmHg) = "+to_string(capPress),0,0);

    setup_continuumArrays();
    distribSource();

/*    optLambda = 15.87;
    optBeta = 3.95*1e-3;
    optKappa = optBeta/pow(optLambda,2);*/

    optKappa = cell.kappa; // mm3.s/kg

    // Convert mmHg to kg/mm.s2
    Pa_Pv *= alpha;
    capPress *= alpha;

    NewtRaph();

    printText("Populating tissue matrix");
    //tissueMat(optLambda);

    //capHomogenisation();

    discreteNet.cnode.row(1) *= cell.aniScaleY;
    discreteNet.cnode.row(2) *= cell.aniScaleZ;

}


void DiscreteContinuum::distribSource() {

    // Define the radial distance between sources points (working in mm)
    int cntr1 = 0;
    int cntr2 = 0;
    vec x, y;
    for (int inodbc = 0; inodbc < discreteNet.getNnodbc(); inodbc++) {
        if (sourceBCtyp(inodbc) == 2 || sourceBCtyp(inodbc) == 4)   {
            cntr2 = 0;
            for (int jnodbc = 0; jnodbc < discreteNet.getNnodbc(); jnodbc++) {
                if (sourceBCtyp(jnodbc) == 2 || sourceBCtyp(jnodbc) == 4)   {
                    x = discreteNet.cnode.col(discreteNet.bcnod(inodbc));
                    y = discreteNet.cnode.col(discreteNet.bcnod(jnodbc));
                    rnod(inodbc-cntr1,jnodbc-cntr2) = eucDistance(x, y) * 1e-3;
                }
                else    {cntr2 += 1;}
            }
            for (int iseg = 0; iseg < discreteNet.getNseg(); iseg++) {
                if (discreteNet.bcnod(inodbc) == discreteNet.ista(iseg) || discreteNet.bcnod(inodbc) == discreteNet.iend(iseg)) {
                    r0(inodbc - cntr1) = discreteNet.rseg(iseg);
                }
            }
            snode.col(inodbc - cntr1) = discreteNet.cnode.col(discreteNet.bcnod(inodbc))*1e-3;
        }
        else {cntr1 += 1;}
    }

    // Initialise source r0 based on associate vascular r0
    bool stable = false;
    int iter{},nitmax=1e3;
    ivec fixed = zeros<ivec>(nnodT);
    vec oldr0 = r0;
/*    while (!stable && iter < nitmax) {



        oldr0 = r0;

    }*/

    // Calculate r0
    double rmin = nonzeros(rnod).min();
/*    for(int i = 0; i < nnodT; i++){
        for(int j = 0; j < nnodT; j++){
            if(i != j && rnod(i,j) < rmin && rnod(i,j) > 0.){rmin = rnod(i,j);}
        }
    }*/
    //r0 = 0.5*rmin;

    //r0.fill(r0);
    //cout<<r0<<endl;

    // Pack some circles!
    //r0.fill(1e3);
    //r0 *= 1e3;
    //sphere_packing(0.95, r0, snode);

}


void DiscreteContinuum::growSource()    {



}


void DiscreteContinuum::tissueMat(const double &lambda) {

    Mtiss.zeros();
    int cntr1{},cntr2{};
    for (int inodbc = 0; inodbc < discreteNet.getNnodbc(); inodbc++) {
        cntr2 = 0;
        if (sourceBCtyp(inodbc) == 2 || sourceBCtyp(inodbc) == 4)   {
            for (int jnodbc = 0; jnodbc < discreteNet.getNnodbc(); jnodbc++) {
                if (sourceBCtyp(jnodbc) == 2 || sourceBCtyp(jnodbc) == 4)   {
                    if (sourceTree(inodbc) == sourceTree(jnodbc))   {
                        Mtiss(inodbc-cntr1,jnodbc-cntr2) = evalTiss(rnod(inodbc-cntr1,jnodbc-cntr2), r0(inodbc - cntr1), lambda);
                    }
                }
                else {cntr2 += 1;}
            }
        }
        else    {cntr1 += 1;}
    }

}


// Function to evaluate each entry for Mtiss matrix
double DiscreteContinuum::evalTiss(const double &r, const double &r0, const double &lambda)   {

    int n = 0;
    int m = 1;
    double val = lambda * r0;
    double i0 = SPHI(n, val);
    double i1 = SPHI(m, val);
    double k0 = SPHK(n, val);
    double k1 = SPHK(m, val);
    double C1 = i1 / (i1*k0+i0*k1);
    double C2 = k1 / (i1*k0+i0*k1);

    double Mij{};
    double lamR = lambda * r;
    if (r < r0) {Mij = (1-C2*SPHI(n,lamR)) * (3/(4*M_PI*pow(r0,3)));}
    else    {Mij = C1*SPHK(n,lamR) * (3/(4*M_PI*pow(r0,3)));}

    return Mij;

}


void DiscreteContinuum::capHomogenisation() {

    // Calculate inverse matrix as in eq 3.23
    mat invMat = inv(Mnet + Mtiss/optBeta);


    // Calculate the flow values at the source nodes
    qout = invMat * (Pa_Pv-capPress);
    qsum = sum(qout) / gamma;

    cout<<"\t\t\t\tTotal Inflow: "<<sum(qout)/gamma<<" nl/min"<<endl;


    // Solving pressure solutions in the continuum domain
    continuumPress = capPress + (Mtiss * qout / optBeta) / alpha;

    if (fabs(qsum - qact) > 1e-3) {
        cout<<"\t\t\t*** Warning : The total flow balance been arteriola and venula trees has not been conserved "<<qact<<",\t qsum: "<<qsum<<endl;
    }
    else    {
        cout<<"\t\t\t\tQact = "<<qact<<" nl/min"<<endl;
    }

    pout = Pa_Pv - Mnet*qout/alpha;

    cout<<accu(qout(find(qout < 0.0)))/gamma<<endl;
    cout<<accu(qout(find(qout > 0.0)))/gamma<<endl;

}


void DiscreteContinuum::NewtRaph()    {

    double logbeta{},lambda{};
    double logbetaLOW = -20.;
    double logbetaHIGH = 10.;
    double logbetaACC = 1e-3;
    double dbeta{},dQdB_SUM{};

    vec dQdB_INT = zeros<vec>(nnodT);
    vec dQdB = zeros<vec>(nnodT);
    mat N = zeros<mat>(nnodT,nnodT);
    mat NINV = zeros<mat>(nnodT,nnodT);

    dbeta = logbetaACC + 1.;
    logbeta = -1e-6;

    qsum = 0.0;
    int nitmax = 1e2;
    int iter = 1;
    double relax = 1.;
    double iBeta{};
    cout<<max(Pa_Pv/alpha)<<endl;
    while((abs(dbeta) > logbetaACC || abs(qsum-qact) > 1.e-3)  && iter <= nitmax)   {

        printText(to_string(iter)+"/"+to_string(nitmax)+": ",1, -1);

        //if (iter % 5 == 0)  {relax *= 0.8;}

        iBeta = pow(10.,logbeta);
        lambda = sqrt(iBeta/optKappa);

        tissueMat(lambda);
        NINV = inv((Mnet + Mtiss/iBeta));
        mat N = (Mnet + Mtiss/iBeta);

        // Calculating the flow values at the source nodes
        qout = NINV * (Pa_Pv-capPress);// solve(N,(Pa_Pv-capPress));//NINV * (Pa_Pv-capPress);
        qsum = accu(qout)/gamma;

        // Calculating the pressures at the source nodes
        continuumPress = capPress + (Mtiss*qout)/iBeta;
        pout = Pa_Pv - (Mnet*qout);

        double rmsPErr = sqrt(mean(pow((continuumPress-pout)/alpha,2)));
        if (rmsPErr > 1.e-8)   {printText("Mtiss/Mnet. -> RMS pressure error = " + to_string(rmsPErr),4);}

        dQdB_INT = Mtiss*qout;
        dQdB = NINV*dQdB_INT;
        dQdB_SUM = pow(1./iBeta,2) * accu(dQdB);

        dbeta = (1./(iBeta * log(10.))) * gamma * (qsum-qact) / dQdB_SUM;
        logbeta -= dbeta;

        /*if (logbeta < -10.)  {logbeta *= 1e-6;}
        else if (logbeta > 1.) {logbeta *= -0.2;}*/

        if ((logbetaLOW-logbeta)*(logbeta-logbetaHIGH) < 1e-6)   {printText("Newton-Raphson Method -> logbeta = " + to_string(logbeta),4);}
        printText( "Flow error = "+to_string(qsum-qact)+", lambda = "+to_string(lambda)+", beta =  "+to_string(pow(10,logbeta)),1,0);

        iter += 1;
    }

    optBeta = pow(10,logbeta);
    optLambda = sqrt(optBeta/optKappa);
    printText( "Final: lambda (1/mm) = "+to_string(optLambda)+", beta = "+to_string(optBeta)+", kappa (mm3.s/kg) =  "+to_string(optKappa),1,0);

    /*cout<<continuumPress/alpha<<endl;
    cout<<"hi"<<endl;
    cout<<pout/alpha<<endl;*/
}
