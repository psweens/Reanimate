#include "Network.hpp"
#include "spatGraph.hpp"

using namespace reanimate;

int Network::readAmira(const string &filename, const string &networkname, bool stubs) {

    const char* errormess = "Something went wrong while reading the binary data section.\nPremature end of file?\n";

    FILE* fp = fopen((loadPath + filename).c_str(), "rb");
    if(!fp){printf("Could not find %s\n", filename.c_str());}

    printf("Reading %s ...\n", filename.c_str());

    // We read the first 2k bytes into memory to parse the header.
    // The fixed buffer size looks a bit like a hack, and it is one, but it gets the job done.
    char buffer[2048];
    fread(buffer, sizeof(char), 2047, fp);
    buffer[2047] = '\0'; //The following string routines prefer null-terminated strings

    if(!strstr(buffer, "# AmiraMesh BINARY-LITTLE-ENDIAN 2.1")){
        printf("*** Error: Not a AmiraMesh BINARY-LITTLE-ENDIAN 2.1 file ***\n");
        fclose(fp);
    }

    sscanf(FindAndJump(buffer, "VERTEX"), "%d %d %d", &nvertex);
    printf("\tVertices = : %d\n", nvertex);
    sscanf(FindAndJump(buffer, "EDGE"), "%d %d %d", &nedge);
    printf("\tEdges = %d\n", nedge);
    sscanf(FindAndJump(buffer, "POINT"), "%d %d %d", &npoint);
    printf("\tPoints = %d\n", npoint);

    //set up arrays
    vertexCoordinates = zeros<mat>(nvertex,3);
    edgeConnectivity = zeros<imat>(nedge,2);
    nedgePoints = zeros<ivec>(nedge);
    edgePointCoordinates = zeros<mat>(npoint,3);
    thickness = zeros<vec>(npoint);

    //Find the beginning of the data section
    const long idxStartData = strstr(buffer, "# Data section follows") - buffer;
    if (idxStartData > 0){
        fseek(fp, idxStartData, SEEK_SET);  //Set the file pointer to the beginning of "# Data section follows"
        fgets(buffer, 128, fp);	//Consume this line, which is "# Data section follows"
        fgets(buffer, 128, fp);	//Consume the next line, which is "@1"
        printf("Should be @1: %s", buffer);

        //Read the data
        const size_t NumToRead = nvertex*3;
        float* pData = new float[NumToRead];    // - prepare memory; use malloc() if you're using pure C
        if (pData){
            const size_t ActRead = fread((void*)pData, sizeof(float), NumToRead, fp);		// - do it
            if (NumToRead != ActRead){		// - ok?
                printf(errormess);
                delete[] pData;
                fclose(fp);
                return 0;
            }
            printf("Putting vertex coordinates in an array ... ");
            for (int k = 0; k < nvertex; k++) {
                for (int c = 0; c < 3; c++) {
                    vertexCoordinates(k,c) = pData[k*3 + c];
                }
            }
            delete[] pData;
        }
        fgets(buffer, 128, fp);
        fgets(buffer, 128, fp);
        printf("Should be @2: %s", buffer);

        //Read the data
        const size_t NumToRead1 = nedge*2;
        int* pData1 = new int[NumToRead1];
        if (pData1)  {
            const size_t ActRead1 = fread((void*)pData1, sizeof(int), NumToRead1, fp);
            if (NumToRead1 != ActRead1) {
                printf(errormess);
                delete[] pData1;
                fclose(fp);
                return 0;
            }
            printf("Putting edge connectivity in an array ... ");
            for (int k = 0; k < nedge; k++) {
                for (int c = 0; c < 2; c++) {
                    edgeConnectivity(k,c) = pData1[k*2 + c];
                }
            }
            delete[] pData1;
        }
        fgets(buffer, 128, fp);
        fgets(buffer, 128, fp);
        printf("Should be @3: %s", buffer);

        //Read the data
        const size_t NumToRead2 = nedge;
        int* pData2 = new int[NumToRead2];
        if (pData2) {
            const size_t ActRead2 = fread((void*)pData2, sizeof(int), NumToRead2, fp);
            if (NumToRead2 != ActRead2){
                printf(errormess);
                delete[] pData2;
                fclose(fp);
                return 0;
            }
            printf("Putting the number of edge points in an array ... ");
            for (int k = 0; k < nedge; k++) {
                nedgePoints(k) = pData2[k];
            }
            delete[] pData2;
        }
        fgets(buffer, 128, fp);
        fgets(buffer, 128, fp);
        printf("Should be @4: %s", buffer);

        //Read the data
        const size_t NumToRead3 = npoint*3;
        float* pData3 = new float[NumToRead3];
        if (pData3) {
            const size_t ActRead3 = fread((void*)pData3, sizeof(float), NumToRead3, fp);
            if(NumToRead3 != ActRead3){
                printf(errormess);
                delete[] pData3;
                fclose(fp);
                return 0;
            }
            printf("Putting edge point coordinates in an array ... ");
            for (int k = 0; k < npoint; k++) {
                for (int c = 0; c < 3; c++) {
                    edgePointCoordinates(k,c) = pData3[k*3 + c];
                }
            }
            delete[] pData3;
        }
        fgets(buffer, 128, fp);
        fgets(buffer, 128, fp);
        printf("Should be @5: %s", buffer);

        //Read the data
        const size_t NumToRead4 = npoint;
        float* pData4 = new float[NumToRead4];
        if (pData4){
            const size_t ActRead4 = fread((void*)pData4, sizeof(float), NumToRead4, fp);
            if (NumToRead4 != ActRead4){
                printf(errormess);
                delete[] pData4;
                fclose(fp);
                return 0;
            }
            printf("Putting thickness in an array ...\n");
            for (int k = 0; k < npoint; k++) thickness(k) = pData4[k];
            delete[] pData4;
        }
    }
    fclose(fp);

    processAmira(stubs);
    printNamira(networkname+".txt", networkname);
    printReducedAmira("reducedAmira.am");

    return 0;

}
const char* Network::FindAndJump(const char* buffer, const char* SearchString)
{
    const char* FoundLoc = strstr(buffer, SearchString);
    if (FoundLoc) {return FoundLoc + strlen(SearchString);}
    return buffer;
}

void Network::processAmira(const bool &stubs) {

    printf("Processing Amira spatial graph ...\n");

    // Perform desired modifications to spatial graph
    imat edgeConnectivityNew = zeros<imat>(nedge,2);
    ivec nedgePointsNew = zeros<ivec>(nedge);
    ivec vertexUsed = zeros<ivec>(nvertex);
    ivec edgeUsed = zeros<ivec>(nedge);
    ivec vertexEdge1 = zeros<ivec>(nvertex);
    mat vertexCoordinatesNew = zeros<mat>(nvertex,3);
    mat edgePointCoordinatesNew = zeros<mat>(npoint,3);
    vec thicknessNew = zeros<vec>(npoint);
    vec thicknessMean = zeros<vec>(nedge);

    int nvertexNew{},nedgeNew{},npointNew{};

    int k{};	// Index of point
    for (int iEdge = 0; iEdge < nedge; iEdge++){
        edgeUsed(iEdge) = 1;
        for (int j = 0; j < nedgePoints(iEdge); j++){
            thicknessMean(iEdge) += thickness(k);
            k++;
        }
        thicknessMean(iEdge) = thicknessMean(iEdge) / nedgePoints(iEdge);
        if (thicknessMean(iEdge) < 1.5) {edgeUsed(iEdge) = 0;}	// Reject edge if average < 0.75 micron thickness
        if (edgeConnectivity(iEdge,0) == edgeConnectivity(iEdge,1)) {edgeUsed(iEdge) = 0;}	// Reject edge if it is a closed loop
        k -= nedgePoints(iEdge);  // Reset k
        for (int j = 0; j < nedgePoints(iEdge); j++){
            thickness(k) = thicknessMean(iEdge);		// Set all segments to average thickness
            k++;
        }
    }
    k -= 1;

    //for(i=1; i<=nEDGE; i++)	if(EdgeUsed[i] == 1) for(c=1; c<=2; c++) VertexUsed[EdgeConnectivity[i][c]+1] = 1;
    for (int iEdge = 0; iEdge < nedge; iEdge++)	{
        if (edgeUsed(iEdge) == 1) {
            for(int c = 0; c < 2; c++){	// This version stores the first edge number associated with the vertex
                int iVertex = edgeConnectivity(iEdge,c);
                vertexUsed(iVertex)++;											// and counts the number of associated edges
                if (vertexUsed(iVertex) == 1) {vertexEdge1(iVertex) = iEdge;}
            }
        }
    }

    if (stubs)  {
        for (int iVertex = 0; iVertex < nvertex; iVertex++) {
            if (vertexUsed(iVertex) == 1)    {	// If a vertex has only one edge, compute its length
                int iEdge = vertexEdge1(iVertex);
                double lcumul = 0.;
                for (int j = 0; j < nedgePoints(iEdge); j++)    {
                    double length = 0.;
                    for(int c = 0; c < 3; c++) {length += pow(edgePointCoordinates(k,c) - edgePointCoordinates(k-1,c),2);}
                    length = sqrt(length);
                    lcumul += length;
                }
                if (lcumul / thicknessMean[iEdge] < 1. || lcumul / lthresh < 2.){		// If the length is too short relative to diameter, remove the edge and the vertex
                    edgeUsed(iEdge) = 0;
                    vertexUsed(iVertex) = 0;
                }
            }
        }
    }

    //for(i=1; i<=nVERTEX; i++) if(VertexUsed[i] == 1){
    for (int iVertex = 0; iVertex < nvertex; iVertex++) {
        if (vertexUsed(iVertex)){
            vertexUsed(iVertex) = nvertexNew;	// Create a mapping of old to new vertex numbers
            for (int c = 0; c < 3; c++) {vertexCoordinatesNew(nvertexNew,c) = vertexCoordinates(iVertex,c);}
            nvertexNew++;
        }
    }
    k = -1;	// Index of point
    nedgeNew = -1;
    npointNew = -1;
    for (int iEdge = 0; iEdge < nedge; iEdge++){	// Index of edge
        if (edgeUsed(iEdge) == 1)   {
            nedgeNew++;
            for (int c = 0; c < 2; c++) {edgeConnectivityNew(nedgeNew,c) = vertexUsed(edgeConnectivity(iEdge,c)) - 1;}	// -1? Use mapping
            nedgePointsNew(nedgeNew) = nedgePoints(iEdge);
            for (int j = 0; j < nedgePoints(iEdge); j++){
                k++;
                npointNew++;
                for (int c = 0; c < 3; c++) {edgePointCoordinatesNew(npointNew,c) = edgePointCoordinates(k,c);}
                thicknessNew(npointNew) = thickness(k);
            }
        }
        else {k += nedgePoints(iEdge);}
    }
    nedgeNew += 1;
    npointNew += 1;

    if (nedgeNew <= 0)  {printf("*** Error: no edges ***\n");}


    // Copy new structure to original arrays
    nedge = nedgeNew;
    npoint = npointNew;
    nvertex = nvertexNew;

    vertexCoordinates = vertexCoordinatesNew;
    edgeConnectivity = edgeConnectivityNew;
    nedgePoints = nedgePointsNew;
    edgePointCoordinates = edgePointCoordinatesNew;
    thickness = thicknessNew;



    // Combine segments up to length > lthresh
    npointNew = -1;	// Index of new point
    k = -1;	//index of point
    nedgePointsNew.zeros();
    for (int iEdge = 0; iEdge < nedge; iEdge++){	// Index of edge
        double lcumul = 0.;
        for (int j = 1; j <= nedgePoints(iEdge); j++)    {
            k++;
            double length = 0.;
            if (j > 1)  {
                for (int c = 0; c < 3; c++) {
                    length += pow(edgePointCoordinates(k-1,c) - edgePointCoordinates(k,c),2);
                }
                length = sqrt(length);
                lcumul += length;
            }
            if (j == 1 || j == (nedgePoints(iEdge)) || lcumul > lthresh)  {	// Keep this one
                npointNew++;
                nedgePointsNew(iEdge)++;
                for (int c = 0; c < 3; c++) {edgePointCoordinatesNew(npointNew,c) = edgePointCoordinates(k,c);}
                thicknessNew(npointNew) = thickness(k);
                lcumul = 0.;
            }
        }
    }
    npointNew += 1;
    // Copy new structure to original arrays
    npoint = npointNew;
    nedgePoints = nedgePointsNew;
    edgePointCoordinates = edgePointCoordinatesNew;
    thickness = thicknessNew;

    printf("Reduced vertices = %i\n", nvertex);
    printf("Reduced edges = %i\n", nedge);
    printf("Reduced points = %i\n", npoint);

}

void Network::printNamira(const string &filename, const string &networkname) {

    FILE *ofp;

    nseg = npoint - nedge;
    nnod = nvertex + npoint - 2*nedge;

    double maxEdge = max(edgePointCoordinates.col(0));
    double maxVertex = max(vertexCoordinates.col(0));
    maxEdge > maxVertex ? alx = maxEdge : alx = maxVertex;
    maxEdge = max(edgePointCoordinates.col(1));
    maxVertex = max(vertexCoordinates.col(1));
    maxEdge > maxVertex ? aly = maxEdge : aly = maxVertex;
    maxEdge = max(edgePointCoordinates.col(2));
    maxVertex = max(vertexCoordinates.col(2));
    maxEdge > maxVertex ? alz = maxEdge : alz = maxVertex;

    printf("Writing network file\n");
    ofp = fopen((buildPath + filename).c_str(), "w");
    fprintf(ofp,"%s network derived from Amira Spatial Graph\n",networkname.c_str());
    fprintf(ofp,"%f %f %f box dimensions in microns - adjust as needed\n",alx,aly,alz);
    fprintf(ofp,"32 32 32 number of tissue points in x,y,z directions - adjust as needed\n");
    fprintf(ofp,"10	outer bound distance - adjust as needed\n");
    fprintf(ofp,"100	max. segment length - adjust as needed\n");
    fprintf(ofp,"30	maximum number of segments per node - adjust as needed\n");
    fprintf(ofp,"%i	total number of segments\n",nseg);
    fprintf(ofp,"SegName Type StartNode EndNode Diam   Flow[nl/min]    Hd\n");
    int k{0},segnodname1{},segnodname2{},cnt{},idiammax{};
    double diammax{};

    imat segnodname = zeros<imat>(2,nseg);
    for (int i = 0; i < nedge; i++){
        for (int j = 1; j < nedgePoints(i); j++){
            k++;
            if (j == 1) {segnodname1 = edgeConnectivity(i,0) + 2;}	// node from nvertex
            else {segnodname1 = nvertex + k;}	// node from npoint
            if (j == nedgePoints(i) - 1) {segnodname2 = edgeConnectivity(i,1) +  2;}	// node from nvertex
            else {segnodname2 = nvertex + k + 1;}	// node from npoint
            double diam = fmax(2.*thickness(k),4.);
            if (diam > diammax){
                diammax = diam;
                idiammax = k;
            }
            fprintf(ofp,"%i %i %i %i %f %f %f\n",k,3,segnodname1,segnodname2,diam,0.0,0.45);
            segnodname(0,cnt) = segnodname1;
            segnodname(1,cnt) = segnodname2;
            cnt++;
        }
        k++;
    }

    segnodname = segnodname.submat(0,0,1,cnt-1);
    printf("Max diameter = %f, idx = %i\n",diammax,idiammax);

    if (cnt != nseg) {printf("*** Error: incorrect number of segments ***\n");}

    fprintf(ofp,"%i   number of nodes\n",nnod);
    fprintf(ofp,"Name    x       y       z\n");
    cnt = 0;
    for (int i = 0; i < nvertex; i++)    {	//nodes from nvertex
        fprintf(ofp,"%i %f %f %f\n",i+1,vertexCoordinates(i,0),vertexCoordinates(i,1),vertexCoordinates(i,2));
        cnt++;
    }

    k = -1;
    int inod = nvertex;
    for (int i = 0; i < nedge; i++){	//nodes from npoint
        for (int j = 1; j <= nedgePoints(i); j++){
            k++;
            inod++;
            if (j > 1 && j < nedgePoints(i))    {
                fprintf(ofp,"%i %f %f %f\n",inod,edgePointCoordinates(k,0),
                        edgePointCoordinates(k,1),edgePointCoordinates(k,2));
                cnt++;
            }
        }
    }
    if (cnt != nnod) {printf("*** Error: incorrect number of nodes\n");}

    ista = zeros<ivec>(nseg);
    iend = zeros<ivec>(nseg);
    indexSegmentConnectivity();

    nodtyp = zeros<ivec>(nnod);
    indexNodeConnectivity();

    nnodbc = 0;
    for (int inod = 0; inod < nnod; inod++) {
        if (nodtyp(inod) == 1) {nnodbc += 1;}
    }
    fprintf(ofp,"%i   number of boundary nodes\n",nnodbc);
    fprintf(ofp,"Name    bctyp     bcprfl     bchd\n");
    for (int inod = 0; inod < nnod; inod++) {
        if (nodtyp(inod) == 1)  {
            fprintf(ofp,"%i %i %f %f\n",inod,3,0.0,consthd);
        }
    }

    fclose(ofp);

}


void Network::printReducedAmira(const string &filename) {

    FILE *ofp;

    ofp = fopen((buildPath + filename).c_str(), "wb");
    fprintf(ofp,"# AmiraMesh BINARY-LITTLE-ENDIAN 2.1\n");
    fprintf(ofp,"\n\n");
    fprintf(ofp,"define VERTEX %i\n",nvertex);
    fprintf(ofp,"define EDGE %i\n",nedge);
    fprintf(ofp,"define POINT %i\n",npoint);
    fprintf(ofp,"\n");
    fprintf(ofp,"Parameters {\n");
    fprintf(ofp,"    ContentType \"HxSpatialGraph\"\n");
    fprintf(ofp,"}\n\n");
    fprintf(ofp,"VERTEX { float[3] VertexCoordinates } @1\n");
    fprintf(ofp,"EDGE { int[2] EdgeConnectivity } @2\n");
    fprintf(ofp,"EDGE { int NumEdgePoints } @3\n");
    fprintf(ofp,"POINT { float[3] EdgePointCoordinates } @4\n");
    fprintf(ofp,"POINT { float thickness } @5\n");
    fprintf(ofp,"\n# Data section follows");

    //Write the data
    fprintf(ofp,"\n@1\n");
    const size_t NumToWrite5 = nvertex*3;
    float* pData5 = new float[NumToWrite5];
    for (int k = 0; k < nvertex; k++) {
        for (int c = 0; c < 3; c++) {
            pData5[k*3 + c] = vertexCoordinates(k,c);
        }
    }
    printf("Writing VertexCoordinates in file ...\n");
    fwrite((void*)pData5, sizeof(float), NumToWrite5, ofp);
    delete[] pData5;

    fprintf(ofp,"\n@2\n");
    const size_t NumToWrite6 = nedge*2;
    int* pData6 = new int[NumToWrite6];
    for (int k = 0; k < nedge; k++) {
        for (int c = 0; c < 2; c++) {
            pData6[k*2 + c] = edgeConnectivity(k,c);
        }
    }
    printf("Writing EdgeConnectivity in file ...\n");
    fwrite((void*)pData6, sizeof(int), NumToWrite6, ofp);
    delete[] pData6;

    fprintf(ofp,"\n@3\n");
    const size_t NumToWrite7 = nedge;
    int* pData7 = new int[NumToWrite7];
    for (int k = 0; k < nedge; k++) {
        pData7[k] = nedgePoints(k);
    }
    printf("Writing NumEdgePoints in file ...\n");
    fwrite((void*)pData7, sizeof(int), NumToWrite7, ofp);
    delete[] pData7;

    fprintf(ofp,"\n@4\n");
    const size_t NumToWrite8 = npoint*3;
    float* pData8 = new float[NumToWrite8];
    for (int k = 0; k < npoint; k++) {
        for (int c = 0; c < 3; c++) {
            pData8[k*3 + c] = edgePointCoordinates(k,c);
        }
    }
    printf("Writing EdgePointCoordinates in file ...\n");
    fwrite((void*)pData8, sizeof(float), NumToWrite8, ofp);
    delete[] pData8;

    fprintf(ofp,"\n@5\n");
    const size_t NumToWrite9 = npoint;
    float* pData9 = new float[NumToWrite9];
    for(int k = 0; k < npoint; k++) {
        pData9[k] = thickness(k);
    }
    printf("Writing thickness in file ...\n");
    fwrite((void*)pData9, sizeof(float), NumToWrite9, ofp);
    delete[] pData9;

    fprintf(ofp,"\n");
    fclose(ofp);

}


void Network::printAmira(const string &filename, const mat &extraData, bool smooth, const char *headers[]) {

    if (smooth)     {

        spatGraph graph;
        graph.generate(*this, true, true);

        uvec idx;
        ivec pntsPerEdge = zeros<ivec>(graph.getNseg());
        for (int iseg = 0; iseg < graph.getNseg(); iseg++)  {
            idx = find(graph.segname(iseg) == edgeLabels);
            //if (graph.ista(iseg) == graph.iend(iseg))   {pntsPerEdge(iseg) = (int) idx.n_elem;}
            pntsPerEdge(iseg) = (int) idx.n_elem + 1;
        }
        npoint = accu(pntsPerEdge);


        FILE *ofp1;

        ofp1 = fopen((buildPath + filename).c_str(), "w");
        fprintf(ofp1,"# AmiraMesh 3D ASCII 2.0 \n \n ");
        fprintf(ofp1,"\n define VERTEX %i",graph.getNnod());
        fprintf(ofp1,"\n define EDGE %i",graph.getNseg());
        fprintf(ofp1,"\n define POINT %i",npoint);
        fprintf(ofp1,"\n \n Parameters { \n \t ContentType \"HxSpatialGraph\" \n } \n");

        fprintf(ofp1,"\n VERTEX { float[3] VertexCoordinates } @1");
        fprintf(ofp1,"\n EDGE { int[2] EdgeConnectivity } @2");
        fprintf(ofp1,"\n EDGE { int NumEdgePoints } @3");
        fprintf(ofp1,"\n POINT { float [3] EdgePointCoordinates } @4");
        fprintf(ofp1,"\n POINT { float Radii } @5");
        if ((int) extraData.n_cols == 1)   {
            for (int j = 0; j < (int) extraData.n_cols; j++)  {
                fprintf(ofp1,"\n POINT { float Extra-%i } @%i",j,6+j);
            }
        }
        else {
            for (int j = 0; j < (int) extraData.n_cols; j++)  {
                fprintf(ofp1,"\n POINT { float %s } @%i",headers[j],6+j);
            }
        }
        fprintf(ofp1,"\n\n");


        // Vertex coordinates
        fprintf(ofp1,"@1\n");
        for(int inod = 0; inod < graph.getNnod(); inod++)  {
            fprintf(ofp1,"%.15e %.15e %.15e\n",graph.cnode(0,inod),graph.cnode(1,inod),graph.cnode(2,inod));
        }

        // Connecting vertices
        fprintf(ofp1,"\n@2\n");
        for(int iseg = 0; iseg < graph.getNseg(); iseg++)  {
            fprintf(ofp1,"%i %i\n",(int) graph.ista(iseg),(int) graph.iend(iseg));
        }

        // Number of points per edge
        fprintf(ofp1,"\n@3\n");
        for (int iseg = 0; iseg < graph.getNseg(); iseg++)  {
            fprintf(ofp1,"%i \n",pntsPerEdge(iseg));
        }

        // Coordinates of points - in order of start/end nodes for a segment (as read in network data file)
        fprintf(ofp1,"\n@4\n");
        int nod1{}, nod2{-1}, pnt{};
        uvec snod, enod;
        ivec tag, pntIdx=zeros<ivec>(npoint);
        for (int iseg = 0; iseg < graph.getNseg(); iseg++)  {
            nod2 = -1;
            idx = find(graph.segname(iseg) == edgeLabels);
            tag = zeros<ivec>((int) idx.n_elem);
            snod = find(graph.segnodname(0,iseg) == nodname);
            enod = find(graph.segnodname(1,iseg) == nodname);
            if ((int) idx.n_elem > 1) {
                fprintf(ofp1, "%.15e %.15e %.15e\n", cnode(0, snod(0)), cnode(1, snod(0)), cnode(2, snod(0)));
                nod1 = (int) snod(0);
                pntIdx(pnt) = snod(0);
                pnt += 1;
                while (nod2 != (int) enod(0)) {
                    for (int iseg = 0; iseg < (int) idx.n_elem; iseg++) {
                        if (ista(idx(iseg)) == nod1 && tag(iseg) == 0) {
                            nod2 = iend(idx(iseg));
                            fprintf(ofp1, "%.15e %.15e %.15e\n", cnode(0, nod2), cnode(1, nod2), cnode(2, nod2));
                            nod1 = nod2;
                            tag(iseg) = 1;
                            pntIdx(pnt) = nod2;
                            pnt += 1;
                            iseg = nseg;
                        }
                        else if (iend(idx(iseg)) == nod1 && tag(iseg) == 0) {
                            nod2 = ista(idx(iseg));
                            fprintf(ofp1, "%.15e %.15e %.15e\n", cnode(0, nod2), cnode(1, nod2), cnode(2, nod2));
                            nod1 = nod2;
                            tag(iseg) = 1;
                            pntIdx(pnt) = nod2;
                            pnt += 1;
                            iseg = nseg;
                        }
                    }
                }
            }
            else {
                fprintf(ofp1, "%.15e %.15e %.15e\n", cnode(0, snod(0)), cnode(1, snod(0)), cnode(2, snod(0)));
                pntIdx(pnt) = snod(0);
                pnt += 1;
                fprintf(ofp1, "%.15e %.15e %.15e\n", cnode(0, enod(0)), cnode(1, enod(0)), cnode(2, enod(0)));
                pntIdx(pnt) = enod(0);
                pnt += 1;
            }
        }


        // Segment radii
        fprintf(ofp1,"\n@5\n");
        for(int i = 0; i < npoint; i++)    {fprintf(ofp1,"%.15e\n",pointAverage(i, pntIdx, rseg));}

        // Extra data
        for (int j = 0; j < (int) extraData.n_cols; j++)  {
            fprintf(ofp1,"\n@%i\n",6+j);
            for(int i = 0; i < npoint; i++)    {fprintf(ofp1,"%.15e\n",pointAverage(i, pntIdx, extraData.col(j)));}
        }

        fclose(ofp1);

    }
    else {

        FILE *ofp1;

        ofp1 = fopen((buildPath + filename).c_str(), "w");
        fprintf(ofp1,"# AmiraMesh 3D ASCII 2.0 \n \n ");
        fprintf(ofp1,"\n define VERTEX %i",nnod);
        fprintf(ofp1,"\n define EDGE %i",nseg);
        fprintf(ofp1,"\n define POINT %i",(2*nseg));
        fprintf(ofp1,"\n \n Parameters { \n \t ContentType \"HxSpatialGraph\" \n } \n");

        fprintf(ofp1,"\n VERTEX { float[3] VertexCoordinates } @1");
        fprintf(ofp1,"\n EDGE { int[2] EdgeConnectivity } @2");
        fprintf(ofp1,"\n EDGE { int NumEdgePoints } @3");
        fprintf(ofp1,"\n POINT { float [3] EdgePointCoordinates } @4");
        fprintf(ofp1,"\n POINT { float Radii } @5");
        for (int j = 0; j < (int) extraData.n_cols; j++)  {
            fprintf(ofp1,"\n POINT { float Extra-%i } @%i\n\n",j,6+j);
        }

        // Vertex coordinates
        fprintf(ofp1,"@1\n");
        for(int inod = 0; inod < nnod; inod++)  {
            fprintf(ofp1,"%.15e %.15e %.15e\n",cnode(0,inod),cnode(1,inod),cnode(2,inod));
        }

        // Connecting vertices
        fprintf(ofp1,"\n@2\n");
        for(int iseg = 0; iseg < nseg; iseg++)  {
            fprintf(ofp1,"%i %i\n",(int) ista(iseg),(int) iend(iseg));
        }

        // Number of points per edge
        fprintf(ofp1,"\n@3\n");
        for(int iseg = 0; iseg < nseg; iseg++)  {
            fprintf(ofp1,"%i \n",2);
        }

        // Coordinates of points - in order of start/end nodes for a segment (as read in network data file)
        fprintf(ofp1,"\n@4\n");
        for(int iseg = 0; iseg < nseg; iseg++)  {
            fprintf(ofp1,"%.15e %.15e %.15e\n",cnode(0,ista(iseg)),cnode(1,ista(iseg)),cnode(2,ista(iseg)));
            fprintf(ofp1,"%.15e %.15e %.15e\n",cnode(0,iend(iseg)),cnode(1,iend(iseg)),cnode(2,iend(iseg)));
        }

        // Segment r0
        fprintf(ofp1,"\n@5\n");
        for(int iseg = 0; iseg < nseg; iseg++)    {
            fprintf(ofp1,"%.15e\n",rseg(iseg));
            fprintf(ofp1,"%.15e\n",rseg(iseg));
        }

        // Extra data
        for (int j = 0; j < (int) extraData.n_cols; j++)  {
            fprintf(ofp1,"\n@%i\n",6+j);
            for(int iseg = 0; iseg < nseg; iseg++)    {
                fprintf(ofp1,"%.15e\n",extraData(iseg,0));
                fprintf(ofp1,"%.15e\n",extraData(iseg,0));
            }
        }


        fclose(ofp1);

    }

}


double Network::pointAverage(const int &pnt, const ivec &pntIdx, const vec &param)    {

    double var{};
    for (int iseg= 0; iseg < nseg; iseg++)  {
        if (ista(iseg) == pntIdx(pnt) || iend(iseg) == pntIdx(pnt)) {
            var += param(iseg);
        }
    }

    return var /= nodtyp(pntIdx(pnt));
}
