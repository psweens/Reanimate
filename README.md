# Reanimate: Microvascular Fluid & Mass Transport Solver

Reanimate is a C++ library to simulate fluid and mass transport in microvascular networks. 
The main goal of this project is to provide a simple and flexible framework to simulate blood flow, interstitial flow, and tracer or drug delivery in large microvascular networks. 
These networks are imported as weighted, undirected graphs which represent a network(s) of blood vessels generated either synthetically or by segmenting and skeletonising biomedical images. 

Reanimate implements several mathematical models in research literature including:
* [Poiseuille's Law](https://www.annualreviews.org/doi/10.1146/annurev.fl.25.010193.000245) for axisymmetric, laminar pipe flow in 1D networks with known boundary conditions.
* [Empirical blood viscosity laws](https://journals.physiology.org/doi/full/10.1152/ajpheart.00297.2005) to compute bulk blood viscosity as a function of vessel diameter and haematocrit, thereby capturing the Fåhraeus-Lindqvist effect.
* [Empirical phase separation law](https://www.ahajournals.org/doi/10.1161/01.res.67.4.826?url_ver=Z39.88-2003&rfr_id=ori:rid:crossref.org&rfr_dat=cr_pub%20%200pubmed) to calculate the disproportion distribution of haematocrit at microvascular bifurcations.
* [Numerical update to phase separation law](https://www.pnas.org/content/117/45/27811) which accounts for cell-free layer disruption and recovery.
* [Blood flow estimation model](https://onlinelibrary.wiley.com/doi/10.1111/j.1549-8719.2012.00184.x) to simulate blood flow in a microvascular network with incomplete boundary conditions.
* [Interstitial flow model](https://journals.plos.org/ploscompbiol/article/comments?id=10.1371/journal.pcbi.1006751) which uses a Green's function method to simulate transvascular fluid transport.
