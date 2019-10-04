bakery
===========================

Generic adaptive resolution approach to reverse mapping of polymers

Structure:

 - ``doc/``  - the sphinx documentation
 - ``src/``  - the main code

   - start_backmapping.py - run the MD simulation with backmapping
   - prepare_files.py - prepare hybrid files

 - ``examples/`` - the set of examples

Requirements:

 - `networkx >= 2.0`
 - modified `espressopp` - https://github.com/bakery-cg2at/espressopp
 - `numpy`

If you like to use the method in your work, please cite following materials:

```
@article{doi:10.1021/acs.jctc.6b00595,
  author = {Krajniak, Jakub and Pandiyan, Sudharsan and Nies, Eric and Samaey, Giovanni}
  title = {Generic Adaptive Resolution Method for Reverse Mapping of Polymers from Coarse-Grained to Atomistic Descriptions},
  journal = {Journal of Chemical Theory and Computation},
  doi = {10.1021/acs.jctc.6b00595},
  note ={PMID: 27685340},
  URL = { http://dx.doi.org/10.1021/acs.jctc.6b00595},
}
```
If in your work you willing to reverse map complex polymer structures like polymer networks, please cite following article:
```
@article {JCC:JCC25129,
    author = {Krajniak, Jakub and Zhang, Zidan and Pandiyan, Sudharsan and Nies, Eric and Samaey, Giovanni},
    title = {Reverse mapping method for complex polymer systems},
    journal = {Journal of Computational Chemistry},
    issn = {1096-987X},
    url = {http://dx.doi.org/10.1002/jcc.25129},
    doi = {10.1002/jcc.25129},
    pages = {n/a--n/a},
    keywords = {molecular dynamics, simulation, soft matter, polymer networks, reverse mapping, epoxy, melamine, dendritic, polyester},
}
```

and the software package:

```
@misc{jakub_krajniak_2017_819783,
  author       = {Jakub Krajniak},
  title        = {{bakery-cg2at/bakery: bakery: reverse mapping tool 
                   v2.0}},
  month        = jun,
  year         = 2017,
  doi          = {10.5281/zenodo.819783},
  url          = {https://doi.org/10.5281/zenodo.819783}
}
```

Please note that if you use in your published work any modified version of this software package then you are obliged to
publish your modified code.

[![DOI](https://zenodo.org/badge/DOI/10.5281/zenodo.819783.svg)](https://doi.org/10.5281/zenodo.819783)
[![Build Status](https://travis-ci.org/bakery-cg2at/bakery.svg?branch=devel)](https://travis-ci.org/bakery-cg2at/bakery)
