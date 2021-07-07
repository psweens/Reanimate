//
// Created by Paul Sweeney on 21/02/2021.
//

#include "Vasculature.hpp"
#include "MicroCell.hpp"

#ifndef DiscreteContinuum_hpp
#define DiscreteContinuum_hpp


namespace reanimate   {

    class DiscreteContinuum : public Vasculature {

    public:

        MicroCell cell;
        Vasculature discreteNet;
        spatGraph graph;
        void runHybrid();

        DiscreteContinuum();
        ~DiscreteContinuum();

    protected:

        int nart{},ncap{},nven{},nnodD{},nnodB{},nnodT{},nIterLambda{};
        double qact{},qsum{},capFlow{},artIn{},venOut{},capPress{},capLseg{},artDiam{},capDiam{},venDiam{},optKappa{},optLambda{},optBeta{},lambda{};
        ivec sourceTyp,sourceBCtyp,sourceTree;
        uvec sourceIdx;
        vec qout_act,pout_act,Pa_Pv,qout,pout,continuumPress,r0,segpressPred,nodpressPred,qPred,Mpress_pred,sigmaFlow_pred,sigmaPress_pred,R2_flow,R2_press;
        mat rnod,Mnet,Mtiss,snode,lambdaRange,NINV;

        double evalTiss(const double &r, const double &r0, const double &lambda);
        double evalTissPress(vec &x);
        double evalTissVel(vec &x);

        void setup_hybridArrays();
        void setup_continuumArrays();
        void analyseBranches();
        void addDummies();
        void bridgeFlow();
        void computeDiscrete();
        void computeContinuum();
        void tissueMat(const double &lambda);
        void distribSource();
        void capHomogenisation(double &iBeta, double &lambda);
        void NewtRaph();
        void mapContinuum(Vasculature &network);
        void printDataAnalysis();

    };

}

#endif /* DiscreteContinuum.hpp */
