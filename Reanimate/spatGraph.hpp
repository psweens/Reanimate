//
// Created by sweene01 on 04/07/2020.
//

#include "Network.hpp"

#ifndef spatGraph_hpp
#define spatGraph_hpp

namespace reanimate {

    class spatGraph : public Network {

    public:

        int nTrees{},narterioles,nvenules;
        double Dr{},Dmin{};
        ivec geometry,isBridge,isBridgehead,visited;
        imat InOutlets;
        //template <class CallNetwork>
        void generate(Network network, bool setNoflow=false);
        void defineTrunk();
        void loadTrunks(const string &filepath);
        void analyseTopology(imat predefinedInput=NULL);
        void dfsBasic(int v, int tag);
        void dfsBridge(int v, int p = -1);
        void findBridges();
        void findBridgeheads();
        void traverseGraph();

        spatGraph();
        ~spatGraph();

    private:

        int timer{};
        ivec Pa,Rs,feedNod,dFeedNod,drainNod,flagTree,daughter,Rn,tin,low;
        imat segnod;

        void setup_graphArrays();
        void classifyNetwork(imat &InOutlets, ivec &geometry);
        void internalClassificationLoop(const int &init_seg, const int &classify, ivec &geometry, const int &algo_2, const int &tree, ivec &fTree);

    };

}



#endif /* spatGraph_hpp */
