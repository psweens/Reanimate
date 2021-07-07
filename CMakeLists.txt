cmake_minimum_required(VERSION 3.15)
project(Reanimate)

set(CMAKE_CXX_STANDARD 14)

set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS}  -O3 -std=c++14 -fopenmp -Wall -DARMA_USE_SUPERLU -DARMA_DONT_USE_WRAPPER -DARMA_DONT_USE_HDF5 -DLLVM_TARGETS_TO_BUILD=all -DARMA_OPENMP_THREADS=42")
#set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -std=c++14 -Wall -O3 -fopenmp -DARMA_USE_SUPERLU -DARMA_DONT_USE_WRAPPER")

include_directories(Reanimate)
include_directories(Reanimate/Misc)

#include_directories(/home/linuxbrew/.linuxbrew/include)
#link_directories(/home/linuxbrew/.linuxbrew/lib)

#include_directories(/usr/local/include)
#link_directories(/usr/local/lib)

#link_libraries(-larmadillo)
link_libraries(-llapack)
link_libraries(-lopenblas)
link_libraries(-lsuperlu)

add_executable(Reanimate_Update
        Reanimate/main.cpp
        Reanimate/Network.cpp
        Reanimate/Network.hpp
        Reanimate/Network_Flow_Solvers.cpp
        Reanimate/Network_Analysis.cpp
        Reanimate/Network_Visuals.cpp
        Reanimate/Network_getsetFuncs.cpp
        Reanimate/Network_Flow.cpp
        Reanimate/Network_Aux.cpp
        Reanimate/Network_setupArrays.cpp
        Reanimate/Network_Print.cpp
        Reanimate/Network_Amira.cpp
        Reanimate/Network_Topology.cpp
        Reanimate/Network_Math.cpp
        Reanimate/spatGraph.hpp
        Reanimate/spatGraph.cpp
        Reanimate/spatGraph_Topology.cpp
        Reanimate/spatGraph_Classification.cpp
        Reanimate/Vasculature.cpp
        Reanimate/Vasculature.hpp
        Reanimate/Vasculature_BloodFlow.cpp
        Reanimate/Vasculature_Rheology.cpp
        Reanimate/Vasculature_BoundaryConditions.cpp
        Reanimate/Vasculature_Analysis.cpp
        Reanimate/DiscreteContinuum.cpp
        Reanimate/DiscreteContinuum.hpp
        Reanimate/DiscreteContinuum_Setup.cpp
        Reanimate/DiscreteContinuum_Continuum.cpp
        Reanimate/DiscreteContinuum_Discrete.cpp
        Reanimate/MicroCell.cpp
        Reanimate/MicroCell.hpp
        Reanimate/MicroCell_Designs.cpp
        Reanimate/MicroCell_Flow.cpp
        Reanimate/Tissue.cpp
        Reanimate/Tissue.hpp
        Reanimate/Tissue_Solvers.cpp
        Reanimate/Solute.cpp
        Reanimate/Solute.hpp Reanimate/DiscreteContinuum_Analysis.cpp)
