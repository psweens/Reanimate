#include "Network.hpp"

using namespace reanimate;

int Network::getNseg()    {return nseg;}
int Network::getNnod()    {return nnod;}
int Network::getNnodbc()  {return nnodbc;}

void Network::setNseg(int nseg)  {this->nseg = nseg;}
void Network::setNnod(int nnod)  {this->nnod = nnod;}
void Network::setNnodbc(int nnodbc)  {this->nnodbc = nnodbc;}

void Network::setStackSize(int stackSize) {

    // Default stack size set to 16 MB minimum
    const rlim_t kStackSize = stackSize;
    struct rlimit rl;
    int result = getrlimit(RLIMIT_STACK, &rl);
    if (result == 0)    {
        if (rl.rlim_cur < kStackSize)   {
            rl.rlim_cur = kStackSize;
            result = setrlimit(RLIMIT_STACK, &rl);
            if (result != 0)    {fprintf(stderr, "setrlimit returned result = %d\n", result);}
        }
    }

}
