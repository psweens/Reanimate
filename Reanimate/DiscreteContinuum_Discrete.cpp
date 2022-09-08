//
// Created by Paul Sweeney on 27/02/2021.
//

#include "DiscreteContinuum.hpp"

using namespace reanimate;

void DiscreteContinuum::computeDiscrete()  {

    printText("Generating discrete network matrix");

    nnodD = 0;
    nnodB = 0;
    for (int inodbc = 0; inodbc < discreteNet.getNnodbc(); inodbc++) {
        if (sourceBCtyp(inodbc) == 2)    {nnodB += 1;}
        else if (sourceBCtyp(inodbc) == 4)   {nnodD += 1;}
    }
    nnodT = nnodB + nnodD;

    printText("No. of terminal nodes = "+to_string(nnodB),0,0);
    printText("No. of dummy nodes = "+to_string(nnodD),0,0);

    qout_act = zeros<vec>(nnodT);
    pout_act = zeros<vec>(nnodT);
    Pa_Pv = zeros<vec>(nnodT);
    sourceGeoTyp = zeros<ivec>(nnodT);
    for (int i = 0; i < (int) sourceIdx.n_elem; i++)    {
        qout_act(i) = discreteNet.BCflow(sourceIdx(i));
        pout_act(i) = discreteNet.BCpress(sourceIdx(i));
        for (int jnodbc = 0; jnodbc < discreteNet.getNnodbc(); jnodbc++) {
            if (sourceBCtyp(jnodbc) == -1 && sourceTree(sourceIdx(i)) == sourceTree(jnodbc))  {
                Pa_Pv(i) = discreteNet.BCpress(jnodbc);
                if (Pa_Pv(i) < pout_act(i) && discreteNet.BCgeo(jnodbc) == 1)  {
                    printText("Arteriolar branch : base pressure > 1", 4);
                }
                else if (Pa_Pv(i) > pout_act(i) && discreteNet.BCgeo(jnodbc) == 3)    {
                    printText("Venular branch : base pressure < 1", 4);
                }
                jnodbc = discreteNet.getNnodbc();
            }
        }
        sourceGeoTyp(i) = discreteNet.BCgeo(sourceIdx(i));
    }
    if (any(sourceGeoTyp == 2)) {printText("Source classified as capillary",4);}

    double inflow = accu(discreteNet.BCflow(find(discreteNet.BCflow > 0.)));

    // No flow even at terminal branches in full network
    for (int inodbc = 0; inodbc < discreteNet.getNnodbc(); inodbc++) {
        if (sourceBCtyp(inodbc) == -1)  {
            discreteNet.bctyp(inodbc) = 0;
            discreteNet.bcprfl(inodbc) = discreteNet.BCpress(inodbc);
        }
        else if (sourceBCtyp(inodbc) == -2) {
            discreteNet.bctyp(inodbc) = 1;
            discreteNet.bcprfl(inodbc) = 0.;//discreteNet.BCflow(inodbc);
        }
        else    {
            discreteNet.bctyp(inodbc) = 1;
            discreteNet.bcprfl(inodbc) = 0.;
        }
    }

    // Populate Mnet matrix - sequential solver
    discreteNet.varviscosity = true;
    discreteNet.phaseseparation = false;
    ivec storeBCtyp = discreteNet.bctyp;
    storeBC = discreteNet.bcprfl;
    vec storeBCflow = discreteNet.BCflow;
    vec storeBCpress = discreteNet.BCpress;
    vec storeSegpress = discreteNet.segpress;
    vec storeqq = discreteNet.qq;
    Mnet = zeros<mat>(nnodT,nnodT);

    discreteNet.rheolParams();
    discreteNet.printNetwork("Branch_Network.txt");
    discreteNet.scaleNetwork(1.e-3);
    for (int i = 0; i < nnodT; i++)    {

        if (Pa_Pv(i) - storeBCpress(sourceIdx(i)) > 0.) {discreteNet.bcprfl(sourceIdx(i)) = -1.;}
        else {discreteNet.bcprfl(sourceIdx(i)) = 1.;}

        /*if (discreteNet.BCgeo(sourceIdx(i)) == 1)   {discreteNet.bcprfl(sourceIdx(i)) = -1.;}
        else if (discreteNet.BCgeo(sourceIdx(i)) == 3)   {discreteNet.bcprfl(sourceIdx(i)) = 1.;}*/
        discreteNet.bctyp(sourceIdx(i)) = 1;

        for (int j = 0; j < nnodT; j++)    {
            if (i != j) {
                discreteNet.bcprfl(sourceIdx(j)) = 0.;
                discreteNet.bctyp(sourceIdx(j)) = 1;
            }
        }

        discreteNet.setup_flowArrays();
        discreteNet.splitHD(&Network::fullSolver, graph);
        for (int inodbc = 0; inodbc < discreteNet.getNnodbc(); inodbc++)    {
            discreteNet.BCpress(inodbc) = discreteNet.nodpress(discreteNet.bcnod(inodbc));
        }

        for (int j = 0; j < nnodT; j++)     {Mnet(i,j) = abs(Pa_Pv(j) - discreteNet.BCpress(sourceIdx(j)));}

        discreteNet.bcprfl = storeBC;
        discreteNet.bctyp = storeBCtyp;
        discreteNet.BCflow = storeBCflow;
        discreteNet.BCpress = storeBCpress;

    }

    Mnet *= (alpha/gamma); // mmHg / nl/min -> Kg / mm4 s

    discreteNet.segpress = storeSegpress;
    discreteNet.qq = storeqq;

}