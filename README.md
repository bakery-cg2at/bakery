bakery
===========================

[![DOI](https://zenodo.org/badge/DOI/10.5281/zenodo.819783.svg)](https://doi.org/10.5281/zenodo.819783)

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

# Citation
If you would like to use the method in your work, please cite the following materials:

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
If in your work you are willing to reverse map complex polymer structures like polymer networks, please cite the following article:
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

# Recent publications

- Zhang, Z.; Krajniak, J.; Keith, J. R.; Ganesan, V. Mechanisms of Ion Transport in Block Copolymeric Polymerized Ionic Liquids. ACS Macro Lett. 2019, 8 (9), 1096-1101. DOI: 10.1021/acsmacrolett.9b00478
- Zhang, Z.; Nasrabadi, A. T.; Aryal, D.; Ganesan, V. Mechanisms of Ion Transport in Lithium Salt-Doped Polymeric Ionic Liquid Electrolytes. Macromolecules 2020, 53 (16), 6995-7008. DOI: 10.1021/acs.macromol.0c01444
- Zhang, Z.; Krajniak, J.; Ganesan, V. A Multiscale Simulation Study of Influence of Morphology on Ion Transport in Block Copolymeric Ionic Liquids. Macromolecules 2021, 54 (11), 4997-5010. DOI: 10.1021/acs.macromol.1c00025
- Zhang, Z.; Sass, J.; Krajniak, J.; Ganesan, V. Ion Correlations and Partial Ionicities in the Lamellar Phases of Block Copolymeric Ionic Liquids. ACS Macro Lett. 2022, 11 (11), 1265-1271. DOI: 10.1021/acsmacrolett.2c00401
- Zhang, Z.; Krajniak, J.; Sass, J.; Sachar, H. S.; Marioni, N.; Duncan, T. J.; Ganesan, V. Conductivity and Transference Numbers in Lithium Salt-Doped Block Copolymeric Ionic Liquid Electrolytes. Macromolecules 2023. DOI: 10.1021/acs.macromol.3c01791
