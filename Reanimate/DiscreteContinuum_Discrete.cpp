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
    for (int i = 0; i < (int) sourceIdx.n_elem; i++)    {
        qout_act(i) = discreteNet.BCflow(sourceIdx(i));
        pout_act(i) = discreteNet.BCpress(sourceIdx(i));
        for (int jnodbc = 0; jnodbc < discreteNet.getNnodbc(); jnodbc++) {
            if (sourceBCtyp(jnodbc) == -1 && sourceTree(sourceIdx(i)) == sourceTree(jnodbc))  {
                Pa_Pv(i) = discreteNet.BCpress(jnodbc);
                jnodbc = discreteNet.getNnodbc();
            }
        }
    }


    // No flow even at terminal branches in full network
    for (int inodbc = 0; inodbc < discreteNet.getNnodbc(); inodbc++) {
        if (sourceBCtyp(inodbc) == -1)  {
            discreteNet.bctyp(inodbc) = 0;
            discreteNet.bcprfl(inodbc) = discreteNet.BCpress(inodbc);
        }
        else if (sourceBCtyp(inodbc) == -2) {
            discreteNet.bctyp(inodbc) = 1;
            discreteNet.bcprfl(inodbc) = 0.;//-discreteNet.BCflow(inodbc);
            /*for (int jnodbc = 0; jnodbc < discreteNet.getNnodbc(); jnodbc++)    {
                if (discreteNet.bcnodname(jnodbc) == sourceTree(inodbc))    {
                    discreteNet.bcprfl(inodbc) = discreteNet.BCflow(inodbc) / discreteNet.BCflow(jnodbc);
                    jnodbc = discreteNet.getNnodbc();
                }
            }*/

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
    int n{};
    for (int i = 0; i < nnodT; i++)    {

        discreteNet.bcprfl(sourceIdx(i)) = discreteNet.BCflow(sourceIdx(i)) / abs(discreteNet.BCflow(sourceIdx(i)));
        discreteNet.bctyp(sourceIdx(i)) = 1;

        for (int j = 0; j < nnodT; j++)    {
            if (i != j) {
                discreteNet.bcprfl(sourceIdx(j)) = 0.;
                discreteNet.bctyp(sourceIdx(j)) = 1;
            }
        }

        if (discreteNet.BCgeo(sourceIdx(i)) == 1) {
            n = 0;
            discreteNet.bcprfl(sourceIdx(i)) = -1.;
            if (discreteNet.bcprfl(sourceIdx(i)) == 1.)  {n = 1;}
        }
        else if (discreteNet.BCgeo(sourceIdx(i)) == 3)    {
            n = 1;
            discreteNet.bcprfl(sourceIdx(i)) = 1.;
            if (discreteNet.bcprfl(sourceIdx(i)) == -1.)  {n = 0;}
        }

        discreteNet.setup_flowArrays();
        discreteNet.splitHD(&Network::fullSolver, graph);
        for (int inodbc = 0; inodbc < discreteNet.getNnodbc(); inodbc++)    {
            discreteNet.BCpress(inodbc) = discreteNet.nodpress(discreteNet.bcnod(inodbc));
        }

        Mnet.col(i) = pow(-1.,n)*(Pa_Pv(i) - discreteNet.BCpress(sourceIdx));

        discreteNet.bcprfl = storeBC;
        discreteNet.bctyp = storeBCtyp;
        discreteNet.BCflow = storeBCflow;
        discreteNet.BCpress = storeBCpress;

    }

    for (int i = 0; i < nnodT; i++)    {
        for (int j = 0; j < nnodT; j++)    {
            if (sourceTree(sourceIdx(i)) != sourceTree(sourceIdx(j)))   {
                Mnet(i,j) = 0.;
            }
        }
    }

    Mnet *= (alpha/gamma); // mmHg / nl/min -> Kg / mm4 s

    discreteNet.segpress = storeSegpress;
    discreteNet.qq = storeqq;

}