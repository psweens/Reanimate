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
        ivec geometry;
        imat InOutlets;
        void generate(Network &network);
        void defineTrunk();
        void loadTrunks(const string &filepath);
        void analyseTopology(imat predefinedInput=NULL);

        spatGraph();
        ~spatGraph();

    private:

        ivec Pa,Rs,feedNod,dFeedNod,drainNod,flagTree,daughter,Rn;

        void setup_graphArrays();
        void classifyNetwork(imat &InOutlets, ivec &geometry);
        void internalClassificationLoop(const int &init_seg, const int &classify, ivec &geometry, const int &algo_2, const int &tree, ivec &fTree);

    };

}



#endif /* spatGraph_hpp */
