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

    optKappa = cell.kappa; // mm3.s/kg

    printText("Calculating capillary drainage via Newton-Raphson scheme");
    NewtRaph();
    while (terminateNM) {
        printText("Updating hydraulic conductivity");
        terminateNM = false;
        optKappa *= 0.8;
        NewtRaph();
    }
    int iter{};
    double oldmse{1.e3};
    bool optimiseHybrid = true;
    bool foundkappa = false;
    vec closestNeighbour = zeros<vec>(nnodT);
    for (int i = 0; i < nnodT; i++) {
        vec tmp = rnod.col(i);
        tmp = tmp(find(tmp > 0.));
        closestNeighbour(i) = min(tmp);
    }
    while (optimiseHybrid && iter < nitmax) {
        NewtRaph();
        vec dSourcePress = discreteNet.BCpress(sourceIdx);
        vec maxSquareError = MaxSE(dSourcePress, pout);
        if (!terminateNM)   {foundkappa = true;}
        if (terminateNM && !foundkappa)    {
            printText("Updating hydraulic conductivity");
            optKappa *= 0.8;
            terminateNM = false;
        }
        else if (maxSquareError(0) > 5.*MedianSE(dSourcePress, pout) && oldmse > maxSquareError(0)) {
            cout<<"here"<<endl;
            cout<<maxSquareError(0)<<endl;
            cout<<MedianSE(dSourcePress, pout)<<endl;
            int idx = maxSquareError(1);
            if (dSourcePress(idx) > pout(idx))  {
                if (r0(idx) * 1.2 < closestNeighbour(idx)) {r0(idx) *= 1.2;}
                else {optimiseHybrid = false;}
            }
            else {r0(idx) *= 0.8;}
            if (optimiseHybrid) {packSpheres(idx,true);}
        }
        else {optimiseHybrid = false;}
        iter += 1;
        oldmse = maxSquareError(0);
    }

    printText("Populating tissue matrix");
    tissueMat(optLambda);
    capHomogenisation(optBeta, optLambda);

    discreteNet.cnode.row(1) *= cell.aniScaleY;
    discreteNet.cnode.row(2) *= cell.aniScaleZ;

}


void DiscreteContinuum::distribSource() {

    // Define the radial distance between sources points (working in mm)
    vec x, y;
    for (int i = 0; i < nnodT; i++)    {
        x = discreteNet.cnode.col(discreteNet.bcnod(sourceIdx(i)));
        for (int j = i; j < (int) sourceIdx.n_elem; j++)    {
            y = discreteNet.cnode.col(discreteNet.bcnod(sourceIdx(j)));
            rnod(i,j) = eucDistance(x, y);
            rnod(j,i) = eucDistance(x, y);
        }
        for (int iseg = 0; iseg < discreteNet.getNseg(); iseg++) {
            if (discreteNet.bcnod(sourceIdx(i)) == discreteNet.ista(iseg) || discreteNet.bcnod(sourceIdx(i)) == discreteNet.iend(iseg)) {
                r0(i) = discreteNet.rseg(iseg);
            }
        }
    }


    // Pack source spheres to touching distance
    packSpheres();

}


void DiscreteContinuum::tissueMat(const double &lambda) {

    Mtiss.zeros();
    for (int i = 0; i < nnodT; i++)    {
        for (int j = 0; j < nnodT; j++)    {
            if (sourceTree(sourceIdx(i)) == sourceTree(sourceIdx(j)))   {
                Mtiss(i,j) = evalTiss(rnod(i,j), r0(i), lambda);
            }
        }
    }

}


// Function to evaluate each entry for Mtiss matrix
double DiscreteContinuum::evalTiss(const double &r, const double &r0, const double &lambda)   {

    double lr0 = lambda * r0;
    double i0 = SPHI(0, lr0);
    double i1 = SPHI(1, lr0);
    double k0 = SPHK(0, lr0);
    double k1 = SPHK(1, lr0);
    double C1 = i1 / (i1*k0 + i0*k1);
    double C2 = k1 / (i1*k0 + i0*k1);

    double Mij = 0.;
    double lamR = lambda * r;
    if (r < r0) {Mij = (1.-C2*SPHI(0,lamR)) * (3./(4.*M_PI*pow(r0,3)));}
    else    {Mij = C1*SPHK(0,lamR) * (3./(4.*M_PI*pow(r0,3)));}

    return Mij;

}


void DiscreteContinuum::capHomogenisation(double &iBeta, double &lambda) {

    tissueMat(lambda);
    NINV = inv((Mnet + Mtiss/iBeta));

    // Calculating the flow values at the source nodes
    qout = NINV * (Pa_Pv-capPress) * alpha;
    qsum = accu(qout) / gamma;

    // Calculating the pressures at the source nodes
    continuumPress = capPress + (Mtiss*qout)/(alpha * iBeta);
    pout = Pa_Pv - (Mnet*qout) / alpha;

}


void DiscreteContinuum::NewtRaph()    {

    double logbeta{},lambda{};
    double logbetaLOW = -20.;
    double logbetaHIGH = 10.;
    double logbetaACC = 1.e-3;
    double dbeta{},dQdB_SUM{};

    vec dQdB_INT = zeros<vec>(nnodT);
    vec dQdB = zeros<vec>(nnodT);

    dbeta = logbetaACC + 1.;
    logbeta = -1.e-8;

    qsum = 0.0;
    int nitmax = 1e2;
    int iter = 1;
    double iBeta{};
    double oldlogbeta{};
    while(abs(qsum-qact) >  1.e-3 && iter <= nitmax)   {

        printText(to_string(iter)+"/"+to_string(nitmax)+": ",1, -1);

        iBeta = pow(10.,logbeta);
        oldlogbeta = logbeta;
        lambda = sqrt(iBeta/optKappa);

        capHomogenisation(iBeta, lambda);

        double rmsPErr = sqrt(mean(pow(continuumPress-pout,2)));
        if (rmsPErr > 1.e-8)   {printText("Mtiss/Mnet. -> RMS pressure error = " + to_string(rmsPErr),4);}

        dQdB_INT = Mtiss*qout;
        dQdB = NINV*dQdB_INT;
        dQdB_SUM = log(10.) * pow(iBeta,-2) * accu(dQdB);

        dbeta = (1./iBeta) * gamma * (qsum-qact) / dQdB_SUM;
        logbeta -= dbeta;

        if ((logbetaLOW-logbeta)*(logbeta-logbetaHIGH) < 1.e-6)   {printText("Newton-Raphson Method -> logbeta = " + to_string(logbeta),4);}
        if (logbeta < -10.)  {
            iter = nitmax + 1;
            terminateNM = true;
        }
        printText( "Flow error = "+to_string(qsum-qact)+", lambda = "+to_string(lambda)+", beta =  "+to_string(pow(10.,logbeta)),1,0);

        iter += 1;

    }
    if (terminateNM) {
        printText("Floating-point error -> Newton method terminated",4);
        optBeta = pow(10.,oldlogbeta);
    }
    else {
        optBeta = pow(10.,logbeta);
        if (iter > nitmax)  {
            printText("Conservation of flow not achieved",4);
            terminateNM = true;
        }
    }
    optLambda = sqrt(optBeta/optKappa);
    printText( "Final: lambda (1/mm) = "+to_string(optLambda)+", beta (mm.s/kg) = "+to_string(optBeta)+", kappa (mm3.s/kg) =  "+to_string(optKappa),1,1);

}

double DiscreteContinuum::evalTissPress(vec &x) {

    double p{capPress},dist{},pgreen{};
    vec y;
    if (cell.aniScaleY > 0.)    {x(1) /= cell.aniScaleY;}
    if (cell.aniScaleZ > 0.)    {x(2) /= cell.aniScaleZ;}
    for (int i = 0; i < (int) sourceIdx.n_elem; i++)    {
        y = discreteNet.cnode.col(discreteNet.bcnod(sourceIdx(i)));
        if (cell.aniScaleY > 0.)    {y(1) /= cell.aniScaleY;}
        if (cell.aniScaleZ > 0.)    {y(2) /= cell.aniScaleZ;}
        dist = eucDistance(x, y);
        pgreen = evalTiss(dist, r0(i), optLambda) / (optBeta * alpha);
        if (pgreen < 0.)    {printText("Negative pressure calculated in Green's",5);}
        p += pgreen * qout(i);
    }

    return p;

}


double DiscreteContinuum::evalTissVel(vec &x) {

    double v{},dist{},vgreen{};
    vec y;
    if (cell.aniScaleY > 0.)    {x(1) /= cell.aniScaleY;}
    if (cell.aniScaleZ > 0.)    {x(2) /= cell.aniScaleZ;}
    for (int i = 0; i < nnodT; i++) {

        double lr0 = optLambda*r0(i);
        double i0 = SPHI(0, lr0);
        double i1 = SPHI(1, lr0);
        double k0 = SPHK(0, lr0);
        double k1 = SPHK(1, lr0);

        double C1 = i1/(i1*k0+i0*k1);
        double C2 = k1/(i1*k0+i0*k1);

        y = discreteNet.cnode.col(discreteNet.bcnod(sourceIdx(i)));
        if (cell.aniScaleY > 0.)    {y(1) /= cell.aniScaleY;}
        if (cell.aniScaleZ > 0.)    {y(2) /= cell.aniScaleZ;}
        dist = eucDistance(x, y);

        double ldist = optLambda*dist;
        if (dist < r0(i))   {vgreen = -C2*(3/(4*M_PI*pow(r0(i),3)*optBeta))*(ldist*cos(ldist) - sin(ldist)) / (optLambda*pow(dist,2));}
        else {vgreen = C1*optLambda*SPHK(0,ldist)*(3/(4*M_PI*pow(r0(i),3)*optBeta));}

        if (vgreen < 0.)    {printText("Negative speed calculated in Green's",5);}

        v += vgreen*abs(qout(i));
    }

    return log10(v);
}
