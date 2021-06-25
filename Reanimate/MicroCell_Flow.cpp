//
// Created by Paul Sweeney on 14/02/2021.
//

#include "MicroCell.hpp"

using namespace reanimate;

void MicroCell::computeConductivity(const string cellType, const int iterations, const string filename)  {

    printText("Micro-Cell Module", 3);

    rheolParams();
    cell2D = false;

    if (cellType == "hexCell2D")    {printText("Applying 2D hexagonal cell");}
    else if (cellType == "crossCell3D") {printText("Applying 3D cross cell");}
    else if (cellType == "load")    {printText("Loading custom cell");}

    mat normConductivity = zeros<mat>(3, 3);
    for (int i = 0; i < iterations; i++)   {
        if (cellType == "hexCell2D")    {hexCell2D();}
        else if (cellType == "crossCell2D") {crossCell2D();}
        else if (cellType == "crossCell3D") {crossCell3D();}
        else if (cellType == "load")    {loadMicroCell();}
        flowMicroCell();
        normConductivity += conductivity;
    }
    conductivity = normConductivity / double(iterations);
    conductivity.print();

    // mm3.s/kg
    kappa = conductivity(0,0);
    normConductivity = conductivity / kappa;

    aniScaleY = sqrt(normConductivity(1, 1));
    aniScaleZ = sqrt(normConductivity(2, 2));

    printNum("Conductivity (scalar, mm3.s/kg) =", kappa);
    printNum("Anisotropic y-scaling =", aniScaleY);
    printNum("Anisotropic z-scaling =", aniScaleZ);

    analyseMicroCell();
    printCellAnalysis("MicroCell_Analysis.txt");

}

void MicroCell::flowMicroCell()    {

    // Millimetre scaling
    lseg *= 1e-3;
    diam *= 1e-3;

    // Allocate arrays
    setup_mcFlowArrays();


    // Construct A matrix
    for (int iseg = 0; iseg < nseg; iseg++) {
        for (int inod = 0; inod < nnod; inod++)    {
            if (inod == ista(iseg)) {cA(iseg,inod) = 1.0 / lseg(iseg);}
            else if (inod == iend(iseg))    {cA(iseg,inod) = -1.0 / lseg(iseg);}
        }
    }


    // Construct B matrix
    uvec idx = find(nodtyp != 1);
    int Innod = (int) idx.n_elem;
    cB = zeros<mat>(Innod,nseg);
    for (int inod = 0; inod < Innod; inod++)    {
        for (int iseg = 0; iseg < nseg; iseg++) {
            if ((int) idx(inod) == ista(iseg)) {cB(inod,iseg) = -1.;}
            else if ((int) idx(inod) == iend(iseg))    {cB(inod,iseg) = 1.;}
        }
    }


    // Construct C matrix
    idx = find(Bin == 1);
    for (int inodbc = 0; inodbc < (int) idx.n_elem; inodbc++) {
        int nod = idx(inodbc);
        cC(inodbc,bcnod(nod)) = 1;
        for (int jnodbc = 0; jnodbc < nnodbc; jnodbc++) {
            if (bcPairs(nod) == bcPairs(jnodbc) && nod != jnodbc)    {cC(inodbc,bcnod(jnodbc)) = -1;}
        }
    }


    // Construct E matrix
    for (int inodbc = 0; inodbc < (int) idx.n_elem; inodbc++) {
        int nod = idx(inodbc);
        for (int iseg = 0; iseg < nseg; iseg++)   {
            if (bcnod(nod) == ista(iseg) || bcnod(nod) == iend(iseg)) {cE(inodbc,iseg) = 1;}
        }
        for (int jnodbc = 0; jnodbc < nnodbc; jnodbc++) {
            if (bcPairs(nod) == bcPairs(jnodbc) && nod != jnodbc)    {
                for (int iseg = 0; iseg < nseg; iseg++)   {
                    if (bcnod(jnodbc) == ista(iseg) || bcnod(jnodbc) == iend(iseg)) {cE(inodbc,iseg) = -1;}
                }
            }
        }

    }


    // Construct F matrix
    for (int iseg = 0; iseg < nseg; iseg++)   {
        cF(iseg,ista(iseg)) = 1;
        cF(iseg,iend(iseg)) = 1;
    }


    // Define conductance
    computeConductance();
    conductance = conductance % lseg; // Modified conductance

    // Conservation of internal cell flux
    mat BA = cB * (conductance % cA.each_col());


    // Conservation of boundary cell flux
    mat EA = cE * (conductance % cA.each_col());


    // Volume pressure average
    double omega = alx*aly*alz*1e-9;
    cF = cF.each_col() % ((M_PI/(omega*8)) * (lseg % pow(diam,2)));

    // Flag vessels aligned with axes
    ivec flag = zeros<ivec>(nseg);
    for (int iseg = 0; iseg < nseg; iseg++) {
        if (cnode(1,ista(iseg)) == cnode(1,iend(iseg)) && cnode(2,ista(iseg)) == cnode(2,iend(iseg))) {
            flag(iseg) = 1; // Parallel to x-axis
        }
        else if (cnode(0,ista(iseg)) == cnode(0,iend(iseg)) && cnode(2,ista(iseg)) == cnode(2,iend(iseg)))   {
            flag(iseg) = 2; // Parallel to y-axis
        }
        else if (cnode(0,ista(iseg)) == cnode(0,iend(iseg)) && cnode(1,ista(iseg)) == cnode(1,iend(iseg)))   {
            flag(iseg) = 3; // Parallel to z-axis
        }
        else {
            flag(iseg) = 4; // Diagonal vessels
        }
    }


    int dim{};
    if (cell2D) {dim = 2;}
    else {dim = 3;}
    conductivity.zeros();
    for (int i = 0; i < dim; i++) {

        vec tissForce = zeros<vec>(3);
        tissForce(i) = 1;

        vec nod1, nod2;
        mat unitCL = zeros<mat>(nseg, 3);
        for (int iseg = 0; iseg < nseg; iseg++) {
            nod1 = cnode.col(ista(iseg));
            nod2 = cnode.col(iend(iseg));
            unitCL.row(iseg) = (abs(nod1-nod2) / eucDistance(nod1,nod2)).t();
        }
        vec forcing = unitCL * tissForce;

        // Conservation of internal cell flux
        vec v1 = -cB * (conductance % forcing);

        // Conservation of boundary flux
        vec v3 = -cE * (conductance % forcing);

        // Conservation of boundary pressures
        vec v2 = zeros<vec>(int(nnodbc/2.));

        // Volume pressure average
        vec v4 = zeros<vec>(nseg);

        mat matrix = join_cols(join_cols(join_cols(BA,cC),EA),cF);
        vec vector = join_cols(join_cols(join_cols(v1,v2),v3),v4);
        nodpress = solve(matrix,vector);

        // Pressure - mm
        for (int iseg = 0; iseg < nseg; iseg++) {
            cellSegpress(i,iseg) = (nodpress(ista(iseg)) + nodpress(iend(iseg)))/2.;
        }

        // Conductivity - mm3.s/kg

        vec p_grad = cA*nodpress;
        for (int iseg = 0; iseg < nseg; iseg++) {
            conductivity.col(i) -= (conductance(iseg) * lseg(iseg) / omega) * (p_grad(iseg) - forcing(iseg)) * unitCL.row(iseg).t();
        }


        cellSegpress.row(i) *= 1e3; // mm to um
        char file1[32];
        char file2[32];
        char file3[32];
        sprintf(file1,"MicroCell_Pressure %i.ps",i);
        sprintf(file2,"Micro-Cell/Amira_Pressure %i.txt",i);
        sprintf(file3,"Micro-Cell/Amira_Gradient %i.txt",i);

        diam *= 1e3;
        rseg *= 1e3;
        pictureNetwork(file1, cellSegpress.row(i).t());
        diam *= 1e-3;
        rseg *= 1e-3;

        //net2amira(file2,"Pressure", nnod, nseg, cnode, ista, iend, 1e3*diam/2, temp_p);
        //net2amira(file3,"Pressure_Gradient", nnod, nseg, cnode, ista, iend, 1e3*diam/2, p_grad);

    }


    // Needs to be a velocity matrix (same with flow) for as they are vectors
    for (int iseg = 0; iseg < nseg; iseg++)   {
        //vel(iseg) = (qq(iseg))/(M_PI*pow((0.5*diam(iseg)),2));
    }

    lseg *= 1e3;
    diam *= 1e3;
}
