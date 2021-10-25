# Reanimate: Microvascular Fluid & Mass Transport Solver

Reanimate is a C++ library to simulate fluid and mass transport in microvascular tissue. 
The main goal of this project is to provide a simple and flexible framework to simulate blood flow, interstitial flow, and tracer or drug delivery in large microvascular networks. 
These networks are imported as weighted, undirected graphs which represent a network(s) of blood vessels generated either synthetically or by segmenting and skeletonising biomedical images. 

Reanimate implements several mathematical models in research literature including:
* [Poiseuille's Law](https://www.annualreviews.org/doi/10.1146/annurev.fl.25.010193.000245) for axisymmetric, laminar pipe flow in 1D networks with known boundary conditions.
* [Empirical blood viscosity laws](https://journals.physiology.org/doi/full/10.1152/ajpheart.00297.2005) to compute bulk blood viscosity as a function of vessel diameter and haematocrit, thereby capturing the Fåhraeus-Lindqvist effect.
* [Empirical phase separation law](https://www.ahajournals.org/doi/10.1161/01.res.67.4.826?url_ver=Z39.88-2003&rfr_id=ori:rid:crossref.org&rfr_dat=cr_pub%20%200pubmed) to calculate the disproportion distribution of haematocrit at microvascular bifurcations.
* [Numerical update to phase separation law](https://www.pnas.org/content/117/45/27811) which accounts for cell-free layer disruption and recovery.
* [Blood flow estimation model](https://onlinelibrary.wiley.com/doi/10.1111/j.1549-8719.2012.00184.x) to simulate blood flow in a microvascular network with incomplete boundary conditions.
* [Interstitial flow model](https://journals.plos.org/ploscompbiol/article/comments?id=10.1371/journal.pcbi.1006751) which uses a Green's function method to simulate transvascular fluid transport.

The Reanimate C++ library is the 2nd generation, user-friendly version of the [code](https://zenodo.org/record/1414160#.YXbN7y1Q1bV) which forms the basis of the REANIMATE (**Rea**listic **N**umerical **I**mage-based **M**odelling of Biologic**a**l **T**issue Substrat**e**s) framework published [here](http://www.nature.com/articles/s41551-018-0306-y) and associate mathematical methods [here](https://journals.plos.org/ploscompbiol/article/comments?id=10.1371/journal.pcbi.1006751).

If you wish to cite the Reanimate library in your work, please use the following references:
> [Computational fluid dynamics with imaging of cleared tissue and of in vivo perfusion predicts drug uptake and treatment responses in tumours](http://www.nature.com/articles/s41551-018-0306-y)<br>
> Angela d'Esposito & Paul W. Sweeney et al.

> [Modelling the transport of fluid through heterogeneous, whole tumours in silico](https://journals.plos.org/ploscompbiol/article/comments?id=10.1371/journal.pcbi.1006751)<br>
> Paul W. Sweeney et al.

## Installation

## Contributing
The Reanimate C++ library is an open-source project started by [Dr Paul Sweeney](www.psweeney.co.uk) during his PhD at University College London under the supervision of [Prof. Rebecca Shipley](https://mecheng.ucl.ac.uk/people/profile/dr-rebecca-shipley/) and [Prof. Simon Walker-Samuel](http://simonwalkersamuel.com). The package is maintained by Dr Paul Sweeney on Github. All contributions of all types are most welcome!

Contribution guidelines are available [here](https://github.com/psweens/Reanimate/blob/master/CONTRIBUTION.md) and a list of feature requests can be found [here](https://github.com/psweens/Reanimate/projects).
