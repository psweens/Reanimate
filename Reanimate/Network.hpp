//
//  Network.hpp
//  TEEPEE
//
//  Created by Paul Sweeney on 16/06/2020.
//  Copyright © 2020 Paul Sweeney. All rights reserved.
//

#include <string>
#include <armadillo>

using namespace arma;
using namespace std;

#ifndef Network_hpp
#define Network_hpp

namespace reanimate {

    class Network {

    public:

        string networkName,networkPath,root,loadroot;
        bool unknownBCs;
        int mxx{},myy{},mzz{},nodsegm{},nsol{},nnodfl{},nedge{};
        float alx{},aly{},alz{},lb{},maxl{};
        double targPress{},targStress{},tissperfusion{},inflow{};
        ivec ista,iend,segname,vesstyp,nodname,bcnodname,bctyp,nodtyp,bcnod,BCgeo,zeroflow,edgeLabels;
        vec diam,rseg,lseg,q,qq,hd,bcprfl,bchd,nodpress,BCflow,BCpress,tau,segpress,elseg,ediam;
        imat segnodname,nodnod,nodseg;
        mat cnode,bcp;

        void loadNetwork(const string& filepath);
        void analyse_network();
        void subNetwork(ivec &index);
        void edgeNetwork();
        void pictureNetwork(const string &filename, vec vector, bool logdist = false, int nl=20, bool nodes=false, bool segs=false);
        void fullSolver();
        void estimationSolver();
        void putrank();
        void printNetwork(const string& filename, bool resetDomain = true);
        void printHistogram(const string &filename, mat &data, const field<string> &headers);

        // 'Getter' functions
        int getNseg();
        int getNnod();
        int getNnodbc();

        // 'Setter' functions
        void setNseg(int nseg);
        void setNnod(int nnod);
        void setNnodbc(int nnodbc);

        Network(); // Constructor
        ~Network(); // Destructor

    protected:

        // Dimensional constants - millimetre scaling
        double gamma,alpha,beta,xi;

        // Network
        int nseg{},nnod{},nnodbc{},computeLseg{};

        // Flow parameters
        int estimationarraysize{},nIBnod{},nunknown{};
        double constvisc{},consthd{},mcv{},hdtol{},qtol{},ktau{},oldktau{},kp{},targetpress{};
        ivec nodout,nodrank,nk,unknownnod,storeBCtyp;
        vec conductance,c,qold,hdold,flowsign,oldFlowsign,tau0,oldTau,Qo,B,p0,storeBC,storeBChd,storeHD,oldHd,oldNodpress,oldq;
        sp_mat M,L,K,A,H1,H2,W;

        void setup_networkArrays();
        void setup_flowArrays();
        void setup_estimationArrays();

    };

}



#endif /* Network_hpp */
