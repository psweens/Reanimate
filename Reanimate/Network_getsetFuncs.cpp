#include "Network.hpp"

using namespace reanimate;

int Network::getNseg()    {return nseg;}
int Network::getNnod()    {return nnod;}
int Network::getNnodbc()  {return nnodbc;}

void Network::setNseg(int nseg)  {this->nseg = nseg;}
void Network::setNnod(int nnod)  {this->nnod = nnod;}
void Network::setNnodbc(int nnodbc)  {this->nnodbc = nnodbc;}