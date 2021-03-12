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
        if (sourceTyp(discreteNet.bcnod(inodbc)) == 2)    {nnodB += 1;}
        else if (sourceTyp(discreteNet.bcnod(inodbc)) == 4)   {nnodD += 1;}
    }
    nnodT = nnodB + nnodD;


    int cntr = 0;
    qout_act = zeros<vec>(nnodT);
    pout_act = zeros<vec>(nnodT);
    Pa_Pv = zeros<vec>(nnodT);
    for (int inodbc = 0; inodbc < discreteNet.getNnodbc(); inodbc++) {
        if (sourceBCtyp(inodbc) == 2 || sourceBCtyp(inodbc) == 4)   {
            qout_act(inodbc-cntr) = discreteNet.BCflow(inodbc);
            pout_act(inodbc-cntr) = discreteNet.BCpress(inodbc);
            for (int jnodbc = 0; jnodbc < discreteNet.getNnodbc(); jnodbc++) {
                if (sourceBCtyp(jnodbc) == -1 && sourceTree(inodbc) == sourceTree(jnodbc) && inodbc != jnodbc)  {
                    Pa_Pv(inodbc-cntr) = discreteNet.BCpress(jnodbc);
                }
            }
        }
        else {
            cntr += 1;
        }
    }
    cntr = 0;
    for (int inodbc = 0; inodbc < discreteNet.getNnodbc(); inodbc++)   {
        if (sourceBCtyp(inodbc) == -1)  {
            double basePress = discreteNet.BCpress(inodbc);
            for (int jnodbc = 0; jnodbc < discreteNet.getNnodbc(); jnodbc++)    {
                if (sourceBCtyp(jnodbc) == 2 || sourceBCtyp(jnodbc) == 4)   {
                    if (sourceTree(inodbc) == sourceTree(jnodbc))   {Pa_Pv(jnodbc-cntr) = basePress;}
                }
                else {cntr += 1;}
            }
        }
        cntr = 0;
        cout<<discreteNet.bcnodname(inodbc)<<"\t"<<sourceBCtyp(inodbc)<<"\t"<<discreteNet.BCpress(inodbc)<<endl;
    }
    cntr = 0;



    for (int inodbc = 0; inodbc < discreteNet.getNnodbc(); inodbc++) {
        if (sourceBCtyp(inodbc) == -1)  {
            discreteNet.bctyp(inodbc) = 0;
            discreteNet.bcprfl(inodbc) = discreteNet.BCpress(inodbc);
        }
        else if (sourceBCtyp(inodbc) == 1)  {
            discreteNet.bctyp(inodbc) = 1;
            for (int jnodbc = 0; jnodbc < discreteNet.getNnodbc(); jnodbc++) {
                if (sourceBCtyp(jnodbc) == -1 && sourceTree(inodbc) == sourceTree(jnodbc))  {
                    discreteNet.bcprfl(inodbc-cntr) = discreteNet.BCflow(inodbc) / fabs(discreteNet.BCflow(jnodbc));
                }
            }
        }
        else    {
            discreteNet.bctyp(inodbc) = 1;
            discreteNet.bcprfl(inodbc) = 0.;
        }
    }



    // Populate Mnet matrix - sequential solver
    discreteNet.phaseseparation = false;
    ivec storeBCtyp = discreteNet.bctyp;
    vec storeBC = discreteNet.bcprfl;
    vec storeBCflow = discreteNet.BCflow;
    vec storeSegpress = discreteNet.segpress;
    vec storeqq = discreteNet.qq;
    Mnet = zeros<mat>(nnodT,nnodT);
    mat testBC = zeros<mat>(discreteNet.getNnodbc(),discreteNet.getNnodbc());
    cntr = 0;


    discreteNet.diam *= 1e-3;
    discreteNet.rseg *= 1e-3;
    discreteNet.lseg *= 1e-3;
    discreteNet.hd.fill(0.45);
    discreteNet.rheolParams();
    for (int inodbc = 0; inodbc < discreteNet.getNnodbc(); inodbc++) {

        if (sourceBCtyp(inodbc) == 2 || sourceBCtyp(inodbc) == 4)  {
            discreteNet.bcprfl(inodbc) = discreteNet.BCflow(inodbc)/abs(discreteNet.BCflow(inodbc));
            discreteNet.bctyp(inodbc) = 1;

            for (int jnodbc = 0; jnodbc < discreteNet.getNnodbc(); jnodbc++)    {
                if (inodbc != jnodbc && (sourceBCtyp(jnodbc) == 2 || sourceBCtyp(jnodbc) == 4))   {
                    discreteNet.bcprfl(jnodbc) = 0.0;
                    discreteNet.bctyp(inodbc) = 1;
                }
            }
        }

        testBC.row(inodbc) = discreteNet.bcprfl.t();

        int n = -1;
        if (discreteNet.BCgeo(inodbc) == 1) {
            n = 0;
            if (discreteNet.bcprfl(inodbc) == 1)    {n = 1;}
        }
        else if (discreteNet.BCgeo(inodbc) == 3)    {
            n = 1;
            if (discreteNet.bcprfl(inodbc) == -1)   {n = 0;}
        }


        if (sourceBCtyp(inodbc) == -1 || sourceBCtyp(inodbc) == 1)   {
            cntr += 1;
        }
        else {

            discreteNet.setup_flowArrays();
            discreteNet.splitHD(&Network::fullSolver, graph);
            Mnet.row(inodbc-cntr) = pow(-1,n)*(Pa_Pv(inodbc-cntr) - discreteNet.BCpress(find(sourceBCtyp == 2 || sourceBCtyp == 4)).t());

        }

        discreteNet.bcprfl = storeBC;
        discreteNet.bctyp = storeBCtyp;
        discreteNet.BCflow = storeBCflow;
    }

    cntr = 0;
    for (int inodbc = 0; inodbc < discreteNet.getNnodbc(); inodbc++) {
        if (sourceBCtyp(inodbc) == 2 || sourceBCtyp(inodbc) == 4)   {

            int cntr2 = 0;
            for (int jnodbc = 0; jnodbc < discreteNet.getNnodbc(); jnodbc++) {
                if (sourceBCtyp(jnodbc) == 2 || sourceBCtyp(jnodbc) == 4)   {
                    if (sourceTree(inodbc) != sourceTree(jnodbc))   {
                        Mnet(inodbc-cntr,jnodbc-cntr2) = 0.;
                    }
                }
                else    {
                    cntr2 += 1;
                }
            }

        }
        else    {
            cntr += 1;
        }
    }

    Mnet *= (alpha/gamma); // mmHg / nl/min -> Kg / mm4 s

    discreteNet.segpress = storeSegpress;
    discreteNet.qq = storeqq;

    //cout<<testBC<<endl;

}