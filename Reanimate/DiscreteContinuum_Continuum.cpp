//
// Created by Paul Sweeney on 27/02/2021.
//

#include "DiscreteContinuum.hpp"

using namespace reanimate;

void DiscreteContinuum::computeContinuum()  {

    printText("Setting up continuum component");

    cell.aniScaleY = 0.997;
    cell.aniScaleZ = 0.995;

    printText("Scaling coordinates based on anisotropic micro-cell conductivity");
    if (cell.aniScaleY > 0.)   {discreteNet.cnode.row(1) /= cell.aniScaleY;}
    if (cell.aniScaleZ > 0.)   {discreteNet.cnode.row(2) /= cell.aniScaleZ;}

    printText("Maximum base pressure (mmHg) = "+to_string(Pa_Pv.max()),0);
    printText("Minimum base pressure (mmHg) = "+to_string(Pa_Pv.min()),0,0);
    printText("Mean capillary pressure (mmHg) = "+to_string(capPress),0,0);
    printText("Target capillary flux (nl/min) = "+to_string(qact),0,0);

    setup_continuumArrays();
    distribSource();

    optKappa = 9.5424 * 1e-4;//cell.kappa; // mm3.s/kg

    printText("Calculating capillary drainage via Newton-Raphson scheme");
    NewtRaph();

/*    int iter{},nitmax{10};
    while (terminateNM && iter < nitmax)     {
        iter += 1;
        if (qsum - qact < 0.)   {
            for (int i = 0; i < nnodT; i++) {
                if (discreteNet.BCgeo(sourceIdx(i)) == 3)    {r0(i) *= 0.98;}
            }
        }
        else if (qsum - qact > 0.)  {
            for (int i = 0; i < nnodT; i++) {
                if (discreteNet.BCgeo(sourceIdx(i)) == 1)    {r0(i) *= 0.98;}
            }
        }
        NewtRaph();
    }*/

/*    bool optimiseR0 = true;
    int iter = 0;
    while (optimiseR0 && iter < 1)  {
        for (int i = 0; i < nnodT; i++) {
            if (discreteNet.BCgeo(sourceIdx(i)) == 1 && pout(i) > discreteNet.BCpress(sourceIdx(i))) {
                for (int iseg = 0 ; iseg < discreteNet.getNseg(); iseg++)   {
                    if (discreteNet.bcnod(sourceIdx(i)) == discreteNet.ista(iseg) || discreteNet.bcnod(sourceIdx(i)) == discreteNet.iend(iseg)) {
                        //cout<<r0(i)<<"\t"<<discreteNet.rseg(iseg)<<endl;
                        //if (0.5 * r0(i) > discreteNet.rseg(iseg))   {r0(i) *= 0.9;}
                        r0(i)  = discreteNet.rseg(iseg);
                        iseg = discreteNet.getNseg();
                    }
                }
            }
            if (discreteNet.BCgeo(sourceIdx(i)) == 3 && pout(i) < discreteNet.BCpress(sourceIdx(i))) {
                for (int iseg = 0 ; iseg < discreteNet.getNseg(); iseg++)   {
                    if (discreteNet.bcnod(sourceIdx(i)) == discreteNet.ista(iseg) || discreteNet.bcnod(sourceIdx(i)) == discreteNet.iend(iseg)) {
                        //cout<<r0(i)<<"\t"<<discreteNet.rseg(iseg)<<endl;
                        //if (0.5 * r0(i) > discreteNet.rseg(iseg))   {r0(i) *= 0.9;}
                        r0(i)  = discreteNet.rseg(iseg);
                        iseg = discreteNet.getNseg();
                    }
                }
            }
        }
        NewtRaph();
        iter += 1;
        cout<<"Iteration ... "<<iter<<endl;
    }*/


 /*   int iter = 0;
    while (pout.min() < Pa_Pv.min() && iter < 5)    {
        for (int i = 0; i < nnodT; i++) {
            if (pout(i) < Pa_Pv.min())  {
                r0(i) = r0.min();
            }
        }
        NewtRaph();
        iter += 1;
    }*/

/*while (terminateNM) {
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
        else if ((maxSquareError(0) > 5.*MedianSE(dSourcePress, pout) && oldmse > maxSquareError(0)) || NewtRaphExplosion==true) {
            cout<<"here"<<endl;
            cout<<maxSquareError(0)<<endl;
            cout<<MedianSE(dSourcePress, pout)<<endl;
            int idx = maxSquareError(1);
            if (dSourcePress(idx) > pout(idx))  {
                if (r0(idx) * 1.05 < closestNeighbour(idx)) {r0(idx) *= 1.05;}
                else {optimiseHybrid = false;}
            }
            else {r0(idx) *= 0.95;}
            if (optimiseHybrid) {packSpheres(idx,true);}
        }
        else {optimiseHybrid = false;}
        iter += 1;
        oldmse = maxSquareError(0);
        if (NewtRaphExplosion == true)   {optimiseHybrid = true;}
    }*/

    tissueMat(optLambda);
    capHomogenisation(optBeta, optLambda);

    uvec idx = find(sourceGeoTyp == 1);
    printText("Mean Arteriolar Source Pressure (mmHg) = "+to_string(mean(continuumPress(idx)))+" ± "+to_string(stddev(continuumPress(idx))),0);
    printText("Min / Max Arteriolar Source Pressure (mmHg) = "+to_string(min(continuumPress(idx)))+" / "+to_string(max(continuumPress(idx))),0,0);
    idx = find(sourceGeoTyp == 3);
    printText("Mean Venular Source Pressure (mmHg) = "+to_string(mean(continuumPress(idx)))+" ± "+to_string(stddev(continuumPress(idx))),0,0);
    printText("Min / Max Venular Source Pressure (mmHg) = "+to_string(min(continuumPress(idx)))+" / "+to_string(max(continuumPress(idx))),0,0);

    discreteNet.cnode.row(1) *= cell.aniScaleY;
    discreteNet.cnode.row(2) *= cell.aniScaleZ;

}


void DiscreteContinuum::distribSource() {

    // Define the radial distance between sources points (working in mm)
    vec x, y;
    for (int i = 0; i < nnodT; i++)    {
        x = discreteNet.cnode.col(discreteNet.bcnod(sourceIdx(i)));
        for (int j = i; j < nnodT; j++)    {
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

/*    vec sourceDiams = zeros<vec>(nnodT);
    vec treeFlow = zeros<vec>(nnodT);
    for (int i = 0; i < (int) sourceIdx.n_elem; i++)    {
        for (int jnodbc = 0; jnodbc < discreteNet.getNnodbc(); jnodbc++) {
            if (sourceBCtyp(jnodbc) == -1 && sourceTree(sourceIdx(i)) == sourceTree(jnodbc))  {
                treeFlow(i) = discreteNet.BCflow(jnodbc);
                for (int iseg = 0; iseg < discreteNet.getNseg(); iseg++)    {
                    if (discreteNet.bcnod(jnodbc) == discreteNet.ista(iseg) || discreteNet.bcnod(jnodbc) == discreteNet.iend(iseg)) {
                        sourceDiams(i) = discreteNet.diam(iseg);
                        iseg = discreteNet.getNseg();
                    }
                }
                jnodbc = discreteNet.getNnodbc();
            }
        }
    }*/

    //r0 = r0 % (1./abs(treeFlow));
    //cout<<r0<<endl;

    //r0.fill(min(r0));
    //r0.fill(0.5 * rnod(find(rnod > 0.)).min());

    // Pack source spheres to touching distance
    //packSpheres();

    vec closestNeighbour = zeros<vec>(nnodT);
    for (int i = 0; i < nnodT; i++) {
        vec tmp = rnod.col(i);
        tmp = tmp(find(tmp > 0.));
        closestNeighbour(i) = min(tmp);
        //r0(i) = sourceDiams(i);
        //if (discreteNet.BCgeo(sourceIdx(i)) == 1)    {r0(i) = 0.5 * closestNeighbour(i);}
        r0(i) = 0.5 * closestNeighbour(i);

    }

/*    for (int i = 0; i < nnodT; i++)   {
        if (r0(i) > 0.5 * closestNeighbour(i))    {r0(i) = 0.5 * closestNeighbour(i);}
    }*/

    printText("Source radii = "+to_string(mean(r0))+" ± "+to_string(stddev(r0))+" mm",0);

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
double DiscreteContinuum::evalTiss(const double &r, const double &sourceR0, const double &lambda)   {

    double lr0 = lambda * sourceR0;
    double i0 = SPHI(0, lr0);
    double i1 = SPHI(1, lr0);
    double k0 = SPHK(0, lr0);
    double k1 = SPHK(1, lr0);
    double C1 = i1 / (i1*k0 + i0*k1);
    double C2 = k1 / (i1*k0 + i0*k1);

    double Mij = 0.;
    double lamR = lambda * r;
    if (r < sourceR0) {Mij = (1.-C2*SPHI(0,lamR)) * (3./(4.*M_PI*pow(sourceR0,3)));}
    else    {Mij = C1*SPHK(0,lamR) * (3./(4.*M_PI*pow(sourceR0,3)));}

    return Mij;

}


void DiscreteContinuum::capHomogenisation(double &iBeta, double &lambda) {

    tissueMat(lambda);
    NINV = inv((Mnet + Mtiss/iBeta));

    // Calculating the flow values at the source nodes
    qout = NINV * (Pa_Pv - capPress) * alpha;
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
    logbeta = -1.e-2;

    qsum = 0.0;
    int nitmax = 1e2;
    int iter = 1;
    double iBeta{};
    double oldlogbeta{};
    vec oldPout = zeros<vec>(nnodT);
    vec oldQout = zeros<vec>(nnodT);
    NewtRaphExplosion = false;
    double maxflow{},rmsPErr{};
    if (artIn > abs(venOut))    {maxflow = artIn;}
    else {maxflow = venOut;}
    mat optLog = zeros<mat>(nitmax,2); // Store beta and flow diff. data
    while(iter == 1 || (abs(qsum-qact) >  0.001 * abs(maxflow) && iter <= nitmax))   {

        printText(to_string(iter)+"/"+to_string(nitmax)+": ",1, -1);

        oldPout = pout;
        oldQout = qout;

        iBeta = pow(10.,logbeta);
        oldlogbeta = logbeta;
        lambda = sqrt(iBeta/optKappa);

        capHomogenisation(iBeta, lambda);

        rmsPErr = sqrt(mean(pow(continuumPress-pout,2)));
        if (rmsPErr > 1.e-8)   {printText("Mtiss/Mnet. -> RMS pressure error = " + to_string(rmsPErr),4);}

        if (isnan(qsum - qact) == 0)    {
            optLog(iter-1,0) = qsum - qact;
            optLog(iter-1,1) = logbeta;
        }
        printText( "Flow error = "+to_string(qsum-qact)+", lambda = "+to_string(lambda)+", beta =  "+to_string(pow(10.,logbeta)),1,0);

        dQdB_INT = Mtiss*qout;
        dQdB = NINV*dQdB_INT;
        dQdB_SUM = (1./pow(iBeta,2)) * accu(dQdB);

        dbeta = (1./(iBeta *  log(10.))) * gamma * (qsum-qact) / dQdB_SUM;
        logbeta -= dbeta;

        if ((logbetaLOW-logbeta)*(logbeta-logbetaHIGH) < 1.e-6)   {printText("Newton-Raphson Method -> logbeta = " + to_string(logbeta),4);}
        if (logbeta < -10.)  {
            iter = nitmax + 1;
            terminateNM = true;
        }

        iter += 1;

    }
    int idx{};
    double tmp{1e4};
    if (terminateNM) {
        printText("Floating-point error -> Newton method terminated",4);
        double val{};
        for (int i = 0; i < nitmax; i++)  {
            val = abs(optLog(i,0));
            if (val < tmp && val != 0.0 && isnan(val) == 0.)  {
                tmp = val;
                idx = i;
            }
        }
        optBeta = pow(10,optLog(idx,1));
        //optBeta = pow(10.,oldlogbeta);
        NewtRaphExplosion = true;
        //terminateNM = false;
        pout = oldPout;
        qout = oldQout;
    }
    else {
        optBeta = pow(10.,logbeta);
        if (iter > nitmax)  {
            printText("Conservation of flow not achieved",4);
            terminateNM = true;
            for (int i = 0; i < nitmax; i++)  {
                if (optLog(i,0) < tmp)  {
                    tmp = abs(optLog(i,0));
                    idx = i;
                }
            }
            optBeta = pow(10,optLog(idx,1));
        }
    }
    optLambda = sqrt(optBeta/optKappa);
    printText( "Final: lambda (1/mm) = "+to_string(optLambda)+", beta (mm.s/kg) = "+to_string(optBeta)+", kappa (mm3.s/kg) =  "+to_string(optKappa),1,1);

}

double DiscreteContinuum::evalTissPress(vec &x) {

    int cntr{};
    double p{capPress},dist{},pgreen{};
    vec y;
    if (cell.aniScaleY > 0.)    {x(1) /= cell.aniScaleY;}
    if (cell.aniScaleZ > 0.)    {x(2) /= cell.aniScaleZ;}
    for (int i = 0; i < (int) sourceIdx.n_elem; i++)    {
        y = discreteNet.cnode.col(discreteNet.bcnod(sourceIdx(i)));
        if (cell.aniScaleY > 0.)    {y(1) /= cell.aniScaleY;}
        if (cell.aniScaleZ > 0.)    {y(2) /= cell.aniScaleZ;}
        dist = eucDistance(x, y);
        pgreen = evalTiss(dist, r0(i), optLambda) / optBeta;
        if (pgreen < 0.)    {cntr +=  1;}
        p += pgreen * qout(i) / alpha; // mmHg scaling
    }
    if (cntr > 0)   {printText("Negative Green's at "+to_string(cntr)+"source locations when computing blood pressure",5);}

    return p;

}


double DiscreteContinuum::evalTissVel(vec &x) {

    double v{},dist{},vgreen{};
    vec y;
    if (cell.aniScaleY > 0.)    {x(1) /= cell.aniScaleY;}
    if (cell.aniScaleZ > 0.)    {x(2) /= cell.aniScaleZ;}

    double lr0{},i0{},i1{},k0{},k1{},C1{},C2{},ldist{};
    for (int i = 0; i < nnodT; i++) {

        lr0 = optLambda*r0(i);
        i0 = SPHI(0, lr0);
        i1 = SPHI(1, lr0);
        k0 = SPHK(0, lr0);
        k1 = SPHK(1, lr0);
        C1 = i1/(i1*k0+i0*k1);
        C2 = k1/(i1*k0+i0*k1);

        y = discreteNet.cnode.col(discreteNet.bcnod(sourceIdx(i)));
        if (cell.aniScaleY > 0.)    {y(1) /= cell.aniScaleY;}
        if (cell.aniScaleZ > 0.)    {y(2) /= cell.aniScaleZ;}
        dist = eucDistance(x, y);

        ldist = optLambda*dist;
        if (dist < r0(i))   {vgreen = -C2*(3/(4*M_PI*pow(r0(i),3)*optBeta))*(ldist*cosh(ldist) - sinh(ldist)) / (optLambda*pow(dist,2));}
        else {vgreen = C1*optLambda*SPHK(0,ldist)*(3/(4*M_PI*pow(r0(i),3)*optBeta));}
        v += vgreen*qout(i); // mm/s
    }
    //if (v < 0.)    {printText("Negative speed calculated in Green's",5);}

    return log(v);
}
