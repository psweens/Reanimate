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

        void setBuildPath(const string buildPath, bool deleteFiles=true) override;
        void setLoadPath(const string loadPath) override;

        DiscreteContinuum();
        ~DiscreteContinuum();

    protected:

        bool terminateNM{false};
        int nart{},ncap{},nven{},nnodD{},nnodB{},nnodT{},nIterLambda{};
        double qact{},qsum{},capFlow{},artIn{},venOut{},capPress{},capLseg{},artDiam{},capDiam{},venDiam{},optKappa{},optLambda{},optBeta{};
        ivec sourceTyp,sourceBCtyp,sourceTree,sourceGeoTyp;
        uvec sourceIdx;
        vec qout_act,pout_act,Pa_Pv,qout,pout,continuumPress,r0,segpressPred,nodpressPred,qPred,Mpress_pred,sigmaFlow_pred,sigmaPress_pred,R2_flow,R2_press;
        mat rnod,Mnet,Mtiss,snode,NINV;

        double evalTiss(const double &r, const double &r0, const double &lambda);
        double evalTissPress(vec &x) override;
        double evalTissVel(vec &x) override;

        void setup_hybridArrays();
        void setup_continuumArrays();
        void analyseBranches();
        void addDummies();
        void bridgeFlow();
        void packSpheres(int idx=-1, bool repack=false);
        void computeDiscrete();
        void computeContinuum();
        void tissueMat(const double &lambda);
        void distribSource();
        void capHomogenisation(double &iBeta, double &lambda);
        void NewtRaph();
        void mapContinuum(Vasculature &network);
        void printDataAnalysis();


    private:

        bool NewtRaphExplosion{};

    };

}

#endif /* DiscreteContinuum.hpp */
