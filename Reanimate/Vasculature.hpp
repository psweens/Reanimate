#include "Network.hpp"
#include "spatGraph.hpp"

#ifndef Vasculature_hpp
#define Vasculature_hpp

namespace reanimate {

    class Vasculature : public Network {

    public:

        bool varviscosity,memoryeffects,updateBoundaryHD,loadDeadEnds;
        spatGraph graph;

        void assignBoundaryHD();
        void rheolParams();
        template <typename Call>
        void splitHD(Call solver, spatGraph &hdGraph);
        void bloodFlow(bool varviscosity=true, bool phaseseparation=false, bool memoryeffects=false, bool updateBoundaryHD=false);
        void printSummary();
        void printVisuals();

        Vasculature();
        ~Vasculature();

    protected:

        double viscor(const double &d, const double &hd);
        void computeConductance();
        void empiricalWSS();

    private:

        double vplas{};
        vec viscpar,cpar,bifpar;


        double recovfn(double &len, double &dp);
        void dishem(bool &memoryeffects, Network &graph);
        void woMemory(int &nout, int &nodt, int &segfltot, ivec &segs, vec &flow, Network &graph);
        void wMemory(int &nout, int &nodt, int &segfltot, ivec &segs, vec &flow, Network &graph);
        void iterateFlowDir(spatGraph &hdGraph);
        vec computeFlowError(double &relax);
        void relaxBoundaryHD();
        void mapFlow(Vasculature &Network);
        void analyseVascularFlow();

    };

}




#endif /* Vasculature_hpp */
