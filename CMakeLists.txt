cmake_minimum_required(VERSION 3.15)
project(Reanimate)

set(CMAKE_CXX_STANDARD 14)

set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -std=c++14 -Wall -O3 -fopenmp -DARMA_USE_SUPERLU -DARMA_DONT_USE_HDF5 -DLLVM_TARGETS_TO_BUILD=all")

include_directories(Reanimate)
include_directories(Reanimate/Misc)

#include_directories(/home/linuxbrew/.linuxbrew/include)
#link_directories(/home/linuxbrew/.linuxbrew/lib)

#include_directories(/usr/local/include)
#link_directories(/usr/local/lib)


link_libraries(-llapack)
link_libraries(-lopenblas)
link_libraries(-lsuperlu)

add_executable(Reanimate_Update
        Reanimate/main.cpp
        Reanimate/Network.cpp
        Reanimate/Network.hpp
        Reanimate/Vasculature_BloodFlow.cpp
        Reanimate/Network_Flow_Solvers.cpp
        Reanimate/Vasculature_Rheology.cpp
        Reanimate/Network_Analysis.cpp
        Reanimate/Network_Visuals.cpp
        Reanimate/spatGraph_Classification.cpp
        Reanimate/spatGraph.hpp
        Reanimate/spatGraph.cpp
        Reanimate/Network_getsetFuncs.cpp
        Reanimate/Vasculature.cpp
        Reanimate/Vasculature.hpp
        Reanimate/Network_Flow.cpp
        Reanimate/Network_Aux.cpp
        Reanimate/Network_setupArrays.cpp Reanimate/Network_Print.cpp Reanimate/Tissue.cpp Reanimate/Tissue.hpp Reanimate/Tissue_Solvers.cpp Reanimate/Solute.cpp Reanimate/Solute.hpp Reanimate/Network_Amira.cpp Reanimate/spatGraph_GraphAnalysis.cpp)
