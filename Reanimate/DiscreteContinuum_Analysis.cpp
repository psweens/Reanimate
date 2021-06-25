//
// Created by Paul Sweeney on 21/02/2021.
//

#include "DiscreteContinuum.hpp"

using namespace reanimate;

void DiscreteContinuum::analyseBranches()  {

    printText("Analysing branching structures");

    // Store mean capillary data
    capPress = mean(discreteNet.segpress(find(discreteNet.vesstyp == 2)));
    capLseg = mean(discreteNet.lseg(find(discreteNet.vesstyp == 2)));
    capDiam = mean(discreteNet.diam(find(discreteNet.vesstyp == 2)));
    artDiam = min(discreteNet.diam(find(discreteNet.vesstyp == 1)));
    venDiam = min(discreteNet.diam(find(discreteNet.vesstyp == 3)));

    // Remove capillary segments
    ivec idx = zeros<ivec>(discreteNet.getNseg());
    idx(find(discreteNet.vesstyp == 2)).fill(1);
    DiscreteContinuum clone = *this;
    discreteNet.subNetwork(idx, false, true);
    graph.mapClassification(discreteNet, true); // Rerun tree classification

    // Plot discrete network
    discreteNet.pictureNetwork("Hybrid_BranchingNetwork.ps",conv_to<vec>::from(discreteNet.vesstyp));

    // Setup arrays
    setup_hybridArrays();

    // Source labels: (1) branch trunks; (2) new terminal nodes; (3) intermediate node (for dummies); (4) dummy nodes' (-2) branch boundary nodes
    for (int inod = 0; inod < discreteNet.getNnod(); inod++)    {
        uvec jdx = find(discreteNet.nodname(inod) == clone.discreteNet.nodname);
        if (discreteNet.nodtyp(inod) == 1 && clone.discreteNet.nodtyp(jdx(0)) == 1)  {sourceTyp(inod) = -2;}
        else if (discreteNet.nodtyp(inod) == 1) {sourceTyp(inod) = 2;}
        else if (discreteNet.nodtyp(inod) != clone.discreteNet.nodtyp(jdx(0))) {sourceTyp(inod) = 3;}
    }


    // Add dummy segments to intermediate nodes & index new connectivity
    addDummies();
    discreteNet.setup_networkArrays();
    discreteNet.rseg = 0.5 * discreteNet.diam;
    discreteNet.indexSegmentConnectivity();
    discreteNet.indexNodeConnectivity();
    discreteNet.indexBCconnectivity();

    discreteNet.pictureNetwork("Hybrid_DummyNetwork.ps",discreteNet.diam);

    // Define boundary nodes types
    sourceBCtyp = zeros<ivec>(discreteNet.getNnodbc());
    for (int inodbc = 0; inodbc < discreteNet.getNnodbc(); inodbc++) {
        sourceBCtyp(inodbc) = sourceTyp(discreteNet.bcnod(inodbc));
        for (int i = 0; i < graph.InOutlets.n_rows; i++)    {
            if (discreteNet.bcnodname(inodbc) == graph.nodname(graph.InOutlets(i,3)))    {sourceBCtyp(inodbc) = -1;}
        }
    }


    // Define geometry of dummy nodes & tree allegiance
    sourceTree = zeros<ivec>(discreteNet.getNnodbc());
    sourceTree.fill(-1);
    for (int inodbc = 0; inodbc < discreteNet.getNnodbc(); inodbc++) {
        if (sourceBCtyp(inodbc) == 4)   {
            discreteNet.BCgeo(inodbc) = discreteNet.vesstyp(discreteNet.nodseg(0,discreteNet.bcnod(inodbc)));
        }
        sourceTree(inodbc) = discreteNet.flagTree(discreteNet.nodseg(0,discreteNet.bcnod(inodbc)));
    }
    if (any(sourceTree == -1)) {printText("Boundary node(s) have not been allocated to a tree",4);}

    // Index sources types of (2) and (4)
    sourceIdx = find(sourceBCtyp == 2 || sourceBCtyp == 4);

}


void DiscreteContinuum::addDummies()    {

    printText("Adding dummy segments");

    int jnod = discreteNet.getNnod();
    int jseg = discreteNet.getNseg();
    for (int inod = 0; inod < (int) discreteNet.getNnod(); inod++) {
        if (sourceTyp(inod) == 3) {
            discreteNet.nodpress.resize(jnod+1);
            discreteNet.nodname.resize(jnod+1);
            discreteNet.nodname(jnod) = discreteNet.nodname.max() + 1;
            discreteNet.cnode.resize(3,jnod+1);
            discreteNet.cnode.col(jnod) = discreteNet.cnode.col(inod);
            sourceTyp.resize(jnod+1);
            sourceTyp(jnod) = 4;
            for (int iseg = 0; iseg < (int) discreteNet.getNseg(); iseg++) {
                if (discreteNet.ista(iseg) == inod || discreteNet.iend(iseg) == inod) {
                    discreteNet.segname.resize(jseg+1);
                    discreteNet.segname(jseg) = discreteNet.segname.max() + 1;
                    discreteNet.vesstyp.resize(jseg+1);
                    discreteNet.vesstyp(jseg) = discreteNet.vesstyp(iseg);
                    discreteNet.segpress.resize(jseg+1);
                    discreteNet.segpress(jseg) = discreteNet.segpress(iseg);
                    discreteNet.segnodname.resize(2,jseg+1);
                    discreteNet.segnodname(0,jseg) = discreteNet.nodname(inod);
                    discreteNet.segnodname(1,jseg) = discreteNet.nodname(jnod);
                    discreteNet.diam.resize(jseg+1);
                    //discreteNet.diam(jseg) = capDiam;
                    if (discreteNet.vesstyp(iseg) == 1) {discreteNet.diam(jseg) = artDiam;}
                    else if (discreteNet.vesstyp(iseg) == 3)    {discreteNet.diam(jseg) = venDiam;}
                    discreteNet.rseg.resize(jseg+1);
                    discreteNet.rseg(jseg) = 0.5 * discreteNet.diam(jseg);
                    discreteNet.lseg.resize(jseg+1);
                    discreteNet.lseg(jseg) = 5.0;//capLseg;
                    discreteNet.q.resize(jseg+1);
                    discreteNet.qq.resize(jseg+1);
                    discreteNet.qq(jseg) = discreteNet.qq(iseg);
                    discreteNet.hd.resize(jseg+1);
                    discreteNet.hd(jseg) = discreteNet.hd(iseg);
                    discreteNet.flagTree.resize(jseg+1);
                    discreteNet.flagTree(jseg) = discreteNet.flagTree(iseg);
                    jseg += 1;
                    iseg = discreteNet.getNseg();
                }
            }
            jnod += 1;
        }
    }
    discreteNet.setNnod((int) discreteNet.nodname.n_elem);
    discreteNet.setNseg((int) discreteNet.segname.n_elem);

    int jnodbc = (int) discreteNet.getNnodbc();
    for (int inod = 0; inod < (int) discreteNet.getNnod(); inod++) {
        if (sourceTyp(inod) == 4) {
            discreteNet.bcnodname.resize(jnodbc+1);
            discreteNet.bcnodname(jnodbc) = discreteNet.nodname(inod);
            discreteNet.bctyp.resize(jnodbc+1);
            discreteNet.bctyp(jnodbc) = 1;
            discreteNet.bcprfl.resize(jnodbc+1);
            discreteNet.bcprfl(jnodbc) = 0.;
            discreteNet.bchd.resize(jnodbc+1);
            discreteNet.bchd(jnodbc) = discreteNet.consthd;
            discreteNet.BCgeo.resize(jnodbc+1);
            discreteNet.BCgeo(jnodbc) = 2;
            jnodbc += 1;
        }
    }
    discreteNet.setNnodbc((int) discreteNet.bcnodname.n_elem);

}


void DiscreteContinuum::bridgeFlow()    {

    printText("Calculating dummy flow");

    discreteNet.BCpress = zeros<vec>(discreteNet.getNnodbc());
    discreteNet.BCflow = zeros<vec>(discreteNet.getNnodbc());
    for (int inodbc = 0; inodbc < discreteNet.getNnodbc(); inodbc++) {
        if (sourceBCtyp(inodbc) != 4)    {
            discreteNet.BCpress(inodbc) = discreteNet.nodpress(discreteNet.bcnod(inodbc));
            for (int iseg = 0; iseg < discreteNet.getNseg(); iseg++) {
                if (discreteNet.bcnod(inodbc) == discreteNet.ista(iseg)) {
                    if (discreteNet.nodpress(discreteNet.ista(iseg)) > discreteNet.nodpress(discreteNet.iend(iseg)))  {discreteNet.BCflow(inodbc) = discreteNet.qq(iseg);}
                    else  {discreteNet.BCflow(inodbc) = -discreteNet.qq(iseg);}
                }
                else if (discreteNet.bcnod(inodbc) == discreteNet.iend(iseg))   {
                    if (discreteNet.nodpress(discreteNet.ista(iseg)) < discreteNet.nodpress(discreteNet.iend(iseg)))  {discreteNet.BCflow(inodbc) = discreteNet.qq(iseg);}
                    else  {discreteNet.BCflow(inodbc) = -discreteNet.qq(iseg);}
                }
            }
        }
        else    {
            discreteNet.BCflow(inodbc) = 0.;
            for (int iseg = 0; iseg < discreteNet.getNseg(); iseg++) {
                if (discreteNet.bcnod(inodbc) == discreteNet.ista(iseg))    {
                    discreteNet.BCpress(inodbc) = discreteNet.nodpress(discreteNet.iend(iseg));
                    for (int kseg = 0; kseg < discreteNet.getNseg(); kseg++) {
                        if (iseg != kseg && (discreteNet.iend(iseg) == discreteNet.ista(kseg) || discreteNet.iend(iseg) == discreteNet.iend(kseg)))  {discreteNet.BCflow(inodbc) += discreteNet.q(kseg);}
                    }
                }
                else if (discreteNet.bcnod(inodbc) == discreteNet.iend(iseg))    {
                    discreteNet.BCpress(inodbc) = discreteNet.nodpress(discreteNet.ista(iseg));
                    for (int kseg = 0; kseg < discreteNet.getNseg(); kseg++) {
                        if (iseg != kseg && (discreteNet.ista(iseg) == discreteNet.ista(kseg) || discreteNet.ista(iseg) == discreteNet.iend(kseg)))  {discreteNet.BCflow(inodbc) += discreteNet.q(kseg);}
                    }
                }
            }
            // Adjust for dummies
            if (discreteNet.BCgeo(inodbc) == 1) {discreteNet.BCflow(inodbc) = -abs(discreteNet.BCflow(inodbc));}
            else if (discreteNet.BCgeo(inodbc) == 3)   {discreteNet.BCflow(inodbc) = abs(discreteNet.BCflow(inodbc));}
        }

    }

}
