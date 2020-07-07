#include "Network.hpp"

#ifndef Vasculature_hpp
#define Vasculature_hpp

namespace reanimate {

    class Vasculature : public Network {

    public:

        bool varviscosity,phaseseparation,memoryeffects;

        void bloodFlow(bool varviscosity=true, bool phaseseparation=false, bool memoryeffects=false);

        Vasculature();
        ~Vasculature();

    private:

        double vplas{};
        vec viscpar,cpar,bifpar;

        double viscor(const double &d, const double &hd);
        double recovfn(double len, double dp);
        void dishem(bool &memoryeffects);
        void woMemory(int &nout, int &nodt, int &segfltot, ivec &segs, vec &flow);
        void wMemory(int &nout, int &nodt, int &segfltot, ivec &segs, vec &flow);
        void rheolParams();
        template <typename Call>
        void splitHD(Call solver);
        void iterateFlowDir();

    };

}




#endif /* Vasculature_hpp */
