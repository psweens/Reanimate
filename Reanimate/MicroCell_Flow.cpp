//
// Created by Paul Sweeney on 14/02/2021.
//

#include "MicroCell.hpp"

using namespace reanimate;

void MicroCell::computeConductivity(const string cellType, const int iterations)  {

    printText("Micro-Cell Module", 3);

    rheolParams();
    cell2D = false;

    if (cellType == "hexCell2D")    {printText("Applying 2D hexagonal cell");}
    else if (cellType == "crossCell3D") {
        printText("Generating 3D cross cell");
        crossCell3D();
    }
    else if (cellType == "load")    {printText("Loading custom cell");}

    setup_mcFlowArrays();

    printText("Computing hydraulic conductivity");
    mat normConductivity = zeros<mat>(3, 3);
    progressBar.reset();
    progressBar.set_niter(iterations);
    for (int i = 0; i < iterations; i++)   {
        vascDens = 1.e2;
        crossCell3D();
        setup_mcFlowArrays();
        assignGeometry();
        while (vascDens > 3.25 || vascDens < 2.75)  {
            if (cellType == "hexCell2D")    {hexCell2D();}
            else if (cellType == "crossCell2D") {crossCell2D();}
            else if (cellType == "crossCell3D") {assignGeometry();}
            else if (cellType == "load")    {loadMicroCell();}
        }
        flowMicroCell();
        normConductivity += conductivity;
        progressBar.update();
    }
    cout<<endl;
    conductivity = normConductivity / double(iterations);
    conductivity.print();

    const char *headers[] = {"Diam", "P1", "P2", "P3"};
    mat data = zeros<mat>(nseg, 4);
    data.col(0) = diam;
    data.col(1) = cellSegpress.row(0).t();
    data.col(2) = cellSegpress.row(1).t();
    data.col(3) = cellSegpress.row(2).t();
    printAmira("MicroCell_amiraBloodPress.am", data, false, headers);

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

    lseg *= 1e-3;
    diam *= 1e-3;

    // Construct A matrix
    for (int iseg = 0; iseg < nseg; iseg++) {
        for (int inod = 0; inod < nnod; inod++)    {
            if (inod == ista(iseg)) {cA(iseg,inod) = 1.0 / lseg(iseg);}
            else if (inod == iend(iseg))    {cA(iseg,inod) = -1.0 / lseg(iseg);}
        }
    }

    // Define vessel conductance
    computeConductance();
    conductance = conductance % lseg; // Modified conductance

    // Conservation of internal cell flux
    BA = cB * (conductance % cA.each_col());

    // Conservation of boundary cell flux
    EA = cE * (conductance % cA.each_col());

    // Volume pressure average
    omega = alx*aly*alz * 1e-9; // mm3 scaling
    cF = cF.each_col() % ((M_PI/(omega*8)) * (lseg % pow(diam,2)));

    int dim{};
    if (cell2D) {dim = 2;}
    else {dim = 3;}
    conductivity.zeros();
    vec nods1, nods2, v1, v2, v3, v4;
    // Conservation of boundary pressures
    v2 = zeros<vec>(int(nnodbc/2.));
    // Volume pressure average
    v4 = zeros<vec>(nseg);
    for (int iseg = 0; iseg < nseg; iseg++) {
        nods1 = cnode.col(ista(iseg));
        nods2 = cnode.col(iend(iseg));
        unitCL.row(iseg) = ((nods1-nods2) / eucDistance(nods1,nods2)).t();
    }
    vec flow = zeros<vec>(nseg);
    for (int i = 0; i < dim; i++) {

        tissForce.zeros();
        tissForce(i) = 1.;

        vec forcing = unitCL * tissForce;

        // Conservation of internal cell flux
        v1 = cB * (conductance % forcing);

        // Conservation of boundary flux
        v3 = cE * (conductance % forcing);

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
            flow(iseg) = conductance(iseg) * (p_grad(iseg) - forcing(iseg));
            conductivity.col(i) -= (conductance(iseg) * lseg(iseg) / omega) * (p_grad(iseg) - forcing(iseg)) * unitCL.row(iseg).t();
        }

        // Check conservation conditions
/*        vec checkB = cB*flow;
        vec checkC = cC * nodpress;
        cout<<"checkB"<<endl;
        cout<<checkB<<endl;
        cout<<"checkC"<<endl;
        cout<<checkC<<endl;*/


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
