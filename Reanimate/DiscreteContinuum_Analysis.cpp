//
// Created by sweene01 on 29/06/2021.
//

#include "DiscreteContinuum.hpp"

using namespace reanimate;

void DiscreteContinuum::mapContinuum(Vasculature &network)  {

    printText("Mapping continuum solution to discrete network");

    ivec nodTag = zeros<ivec>(network.getNnod());
    for (int inod = 0; inod < network.getNnod(); inod++)    {
        bool lvess = false;
        bool cvess = false;
        int nod = (int) network.nodtyp(inod);
        for (int jnod = 0; jnod < nod; jnod++)   {
            if (network.vesstyp(network.nodseg(jnod, inod)) == 2)   {cvess = true;}
            else {lvess = true;}
        }
        if (cvess && lvess) {nodTag(inod) = 1;}
    }

    int nod1{},nod2{};
    double press1{},press2{};
    vec cnod1,cnod2;
    segpressPred = zeros<vec>(network.getNseg());
    qPred = zeros<vec>(network.getNseg());
    nodpressPred = zeros<vec>(network.getNnod());
    network.conductance = zeros<vec>(network.getNseg());

    network.scaleNetwork(1.e-3);

    network.rheolParams();
    network.computeConductance();
    ivec flagnod = zeros<ivec>(network.getNnod());
    ivec branchtyp = zeros<ivec>(network.getNnod());
    for (int iseg = 0; iseg < network.getNseg(); iseg++)    {
        if (network.vesstyp(iseg) == 1 || network.vesstyp(iseg) == 3)  {
            segpressPred(iseg) = network.segpress(iseg);
            qPred(iseg) = network.q(iseg);
        }
        else {
            nod1 = (int) network.ista(iseg);
            nod2 = (int) network.iend(iseg);
            cnod1 = network.cnode.col(nod1);
            cnod2 = network.cnode.col(nod2);
            if (nodTag(nod1) == 1)    {press1 = network.nodpress(nod1);}
            else {press1 = evalTissPress(cnod1);}
            if (nodTag(nod2) == 1)    {press2 = network.nodpress(nod2);}
            else {press2 = evalTissPress(cnod2);}
            segpressPred(iseg) = 0.5 * (press1 + press2);
            qPred(iseg) = network.conductance(iseg) * (press1 - press2) * (alpha/gamma);
            flagnod(nod1) = 1;
            flagnod(nod2) = 1;
        }
        nodpressPred(nod1) = press1;
        nodpressPred(nod2) = press2;
    }

    uvec idx = find(network.vesstyp == 2);
    printText("Mean Predicted Capillary Pressure (mmHg) = "+to_string(mean(segpressPred(idx)))+" ± "+to_string(stddev(segpressPred(idx))),0);
    printText("Min / Max Predicted Capillary Pressure (mmHg) = "+to_string(min(segpressPred(idx)))+" / "+to_string(max(segpressPred(idx))),0,0);

    field<string> headers(4,1);
    headers(0,0) = "Discrete_Press";
    headers(1,0) = "Darcy_Press";
    headers(2,0) = "Discrete_Flow";
    headers(3,0) = "Hybrid_Flow";

    mat data = zeros<mat>(idx.n_elem,4);
    data.col(0) = network.segpress(idx);
    data.col(1) = segpressPred(idx);
    data.col(2) = network.qq(idx);
    data.col(3) = abs(qPred(idx));
    printRawData("Hybrid_SegmentData.txt", data, headers);

    field<string> headers2(2,1);
    headers(0,0) = "Discrete_Press";
    headers(1,0) = "Darcy_Press";
    idx = find(flagnod == 1);
    data = zeros<mat>(idx.n_elem,4);
    data.col(0) = network.nodpress(idx);
    data.col(1) = nodpressPred(idx);
    printRawData("Hybrid_NodeData.txt", data, headers);

    field<string> headers3(3,1);
    headers3(0,0) = "Discrete_Press";
    headers3(1,0) = "Hybrid_SourcePress";
    headers3(2,0) = "Node_Classification";
    idx = find(nodTag == 1);
    data = zeros<mat>(idx.n_elem,3);
    data.col(0) = network.nodpress(idx);
    for (int inod = 0; inod < (int) idx.n_elem; inod++)    {
        cnod1 = network.cnode.col(idx(inod));
        data(inod,1) = evalTissPress(cnod1);
        for (int iseg = 0; iseg < network.getNseg(); iseg++)    {
            if ((network.ista(iseg) == idx(inod) || network.iend(iseg) == idx(inod)) && network.vesstyp(iseg) != 2)  {
                data(inod,2) = network.vesstyp(iseg);
                iseg = network.getNseg();
            }
        }
    }
    printRawData("Hybrid_SourceData.txt", data, headers3);

    network.pictureNetwork("Hybrid_DiscretePressure.ps", segpressPred);
    network.pictureNetwork("Hybrid_DiscreteFlow.ps", abs(qPred), false);

    const char *aheaders[] = {"Pressure"};
    data = zeros<mat>(network.getNseg(), 1);
    data.col(0) = segpressPred;
    network.printAmira("Hybrid_DiscretePressure.am", data, true, aheaders);

}
