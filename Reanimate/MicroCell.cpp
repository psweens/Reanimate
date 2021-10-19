//
// Created by Paul Sweeney on 14/02/2021.
//

#include "MicroCell.hpp"
#include <stdio.h>

using namespace reanimate;

MicroCell::MicroCell() = default;
MicroCell::~MicroCell() = default;

void MicroCell::loadMicroCell() {}

void MicroCell::setup_microCellArrays()  {

    nodsegm = 10;
    computeLseg = 1;

    setup_networkArrays();

    segname = zeros<ivec>(nseg);
    nodname = zeros<ivec>(nnod);
    diam = zeros<vec>(nseg);
    lseg = zeros<vec>(nseg);
    q = zeros<vec>(nseg);
    hd = ones<vec>(nseg) * consthd;
    nodpress = zeros<vec>(nnod);

    cnode = zeros<mat>(3,nnod);
    cellSegpress = zeros<mat>(3,nseg);

}


void MicroCell::setup_mcFlowArrays() {

    conductance = zeros<vec>(nseg);
    c = zeros<vec>(nseg);

    cA = zeros<mat>(nseg,nnod);
    cC = zeros<mat>(int(nnodbc/2.),nnod);
    cE = zeros<mat>(int(nnodbc/2.),nseg);
    cF = zeros<mat>(nseg,nnod);
    conductivity = zeros<mat>(3,3);

}

void MicroCell::analyseMicroCell()  {

    // (microns)

    // Network volume
    netVol = sum(lseg % (diam % diam))*M_PI/4;

    // Vascular density
    vascDens = 100*netVol/(alx*aly*alz);

    // Vascular Length Density
    lsegDens = sum(lseg)/(alx*aly*alz)*1e6; // 1/mm2

    // Vascular surface density
    surfDens = sum(diam % lseg)*M_PI/(alx*aly*alz)*1e3; // 1/mm

    // Vascular surface area to volume ratio
    surfVolRatio = sum(diam % lseg)*M_PI/netVol*1e3; // 1/mm

    // Maximum extravascular diffusion distance
    R = 1/sqrt(M_PI*lsegDens*1e-6);
    
}

void MicroCell::printCellAnalysis(string filename) {

    FILE *ofp;

    filename = buildPath + filename;

    ofp = fopen(filename.c_str(),"w");
    fprintf(ofp,"Example statistic of micro-cell");
    fprintf(ofp,"Network Dimensions (um): %.1f x %.1f x %.1f\n",alx,aly,alz);
    fprintf(ofp,"Vascular Density = %.3f %%\n",vascDens);
    fprintf(ofp,"Vascular Length Density = %.3f mm^-2\n",lsegDens);
    fprintf(ofp,"Vascular Surface Density = %.3f mm^-1\n",surfDens);
    fprintf(ofp,"Vascular Surface Area/Vascular Volume = %.3f mm^-1\n",surfVolRatio);
    fprintf(ofp,"Max. Extravascular Diffusion Distance = %.3f um^-1\n",R);
    fprintf(ofp,"Mean Diam. = %.3f ± %.3f um\n",mean(diam),stddev(diam));
    fprintf(ofp,"Min. Diam. = %.3f um\n",diam.min());
    fprintf(ofp,"Max. Diam. = %.3f um\n",diam.max());
    fprintf(ofp,"Mean Length = %.3f ± %.3f um\n",mean(lseg),stddev(lseg));
    fprintf(ofp,"Min. Length = %.3f um\n",lseg.min());
    fprintf(ofp,"Max. Length = %.3f um\n",lseg.max());
    fprintf(ofp,"Min. P1 = %.3f um\n",min(cellSegpress.row(0)));
    fprintf(ofp,"Max. P1 = %.3f um\n",max(cellSegpress.row(0)));
    fprintf(ofp,"Mean. P1 = %.3f ± %.3f um\n",mean(cellSegpress.row(0)),stddev(cellSegpress.row(0)));
    fprintf(ofp,"Min. P2 = %.3f um\n",min(cellSegpress.row(1)));
    fprintf(ofp,"Max. P2 = %.3f um\n",max(cellSegpress.row(1)));
    fprintf(ofp,"Mean. P2 = %.3f ± %.3f um\n",mean(cellSegpress.row(1)),stddev(cellSegpress.row(1)));
    fprintf(ofp,"Min. P3 = %.3f um\n",min(cellSegpress.row(2)));
    fprintf(ofp,"Max. P3 = %.3f um\n",max(cellSegpress.row(2)));
    fprintf(ofp,"Mean. P3 = %.3f ± %.3f um\n",mean(cellSegpress.row(2)),stddev(cellSegpress.row(2)));
    fprintf(ofp,"Cell Conductivity (mm3 s / kg) x 10^-4\n");
    fprintf(ofp,"K11 = %.4e\n",conductivity(0,0)*1e4);
    fprintf(ofp,"K21 = %.4e\n",conductivity(1,0)*1e4);
    fprintf(ofp,"K31 = %.4e\n",conductivity(2,0)*1e4);
    fprintf(ofp,"K12 = %.4e\n",conductivity(0,1)*1e4);
    fprintf(ofp,"K22 = %.4e\n",conductivity(1,1)*1e4);
    fprintf(ofp,"K32 = %.4e\n",conductivity(2,1)*1e4);
    fprintf(ofp,"K13 = %.4e\n",conductivity(0,2)*1e4);
    fprintf(ofp,"K23 = %.4e\n",conductivity(1,2)*1e4);
    fprintf(ofp,"K33 = %.4e\n",conductivity(2,2)*1e4);
    fprintf(ofp,"Y-Scale = %.4f\n",aniScaleY);
    fprintf(ofp,"Z-Scale = %.4f\n",aniScaleZ);

    fclose(ofp);
    
}

void MicroCell::setEdgeDiamDistrib(const vec d) {diamDistrib = d;}
void MicroCell::setEdgeLengthDistrib(const vec l) {lengthDistrib = l;}
void MicroCell::setEucLengthDistrib(const vec l)   {eucLengths = l;}
void MicroCell::setRotAngle(const double angle) {rotationAngle = angle;}