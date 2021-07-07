//
// Created by sweene01 on 29/06/2021.
//

#include "DiscreteContinuum.hpp"

using namespace reanimate;

void DiscreteContinuum::mapContinuum(Vasculature &network)  {

    ivec nodTag = zeros<ivec>(network.getNnod());
    for (int inod = 0; inod < network.getNnod(); inod++)    {
        bool lvess = false;
        bool cvess = false;
        int nod = network.nodtyp(inod);
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
    network.computeConductance();
    for (int iseg = 0; iseg < network.getNseg(); iseg++)    {
        if (network.vesstyp(iseg) == 1 || network.vesstyp(iseg) == 3)  {
            segpressPred(iseg) = network.segpress(iseg);
            qPred(iseg) = network.q(iseg);
        }
        else {
            nod1 = network.ista(iseg);
            nod2 = network.iend(iseg);
            cnod1 = network.cnode.col(nod1);
            cnod2 = network.cnode.col(nod2);
            if (nodTag(nod1) == 1)    {press1 = network.nodpress(nod1);}
            else {press1 = evalTissPress(cnod1);}
            if (nodTag(nod2) == 1)    {press2 = network.nodpress(nod2);}
            else {press2 = evalTissPress(cnod2);}
            segpressPred(iseg) = 0.5 * (press1 + press2);
            qPred(iseg) = network.conductance(iseg) * (press1 - press2) * (alpha/gamma);
        }
        nodpressPred(nod1) = press1;
        nodpressPred(nod2) = press2;
    }

    field<string> headers(4,1);
    headers(0,0) = "Discrete_Press";
    headers(1,0) = "Hybrid_Press";
    headers(2,0) = "Discrete_Flow";
    headers(3,0) = "Hybrid_Flow";
    uvec idx = find(network.vesstyp == 2);
    mat data = zeros<mat>(idx.n_elem,4);
    data.col(0) = network.segpress(idx);
    data.col(1) = segpressPred(idx);
    data.col(2) = network.qq(idx);
    data.col(3) = abs(qPred(idx));
    printRawData("Hybrid_PressFlowData.txt", data, headers);

    field<string> headers2(2,1);
    headers2(0,0) = "Discrete_Press";
    headers2(1,0) = "Hybrid_SourcePress";
    idx = find(nodTag == 1);
    data = zeros<mat>(idx.n_elem,2);
    data.col(0) = network.nodpress(idx);
    for (int inod = 0; inod < (int) idx.n_elem; inod++)    {
        cnod1 = network.cnode.col(idx(inod));
        data(inod,1) = evalTissPress(cnod1);
    }
    printRawData("Hybrid_SourceData.txt", data, headers2);

    network.pictureNetwork("Hybrid_DiscretePressure.ps", segpressPred);
    network.pictureNetwork("Hybrid_DiscreteFlow.ps", abs(qPred));

}
