bakery
===========================

Generic adaptive resolution approach to reverse mapping of polymers

Structure:

 - ``doc/``  - the sphinx documentation
 - ``src/``  - the main code

   - start_backmapping.py - run the MD simulation with backmapping
   - start_simulation.py - run normal MD simulation
   - prepare_files.py - prepare hybrid files

 - ``examples/`` - the set of examples

Requirments:

 - `networkx >= 1.11`
 - modified `espressopp` - https://github.com/MrTheodor/espressopp/tree/Adress-newdd-rebased
 - `numpy`

If you like to use the method in your project, modify and reuse code of `bakery` or the ingredients from `espressopp`
that are related to the reverse mapping, please cite following materials:

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

```
@misc{jakub_krajniak_2016_61233,
  author       = {Jakub Krajniak},
  title        = {bakery: reverse mapping tool v1.3},
  month        = aug,
  year         = 2016,
  doi          = {10.5281/zenodo.61233},
  url          = {http://dx.doi.org/10.5281/zenodo.61233}
}
```

[![DOI](https://zenodo.org/badge/20122/MrTheodor/bakery.svg)](https://zenodo.org/badge/latestdoi/20122/MrTheodor/bakery)
[![Build Status](https://travis-ci.org/bakery-cg2at/bakery.svg?branch=master)](https://travis-ci.org/bakery-cg2at/bakery)
