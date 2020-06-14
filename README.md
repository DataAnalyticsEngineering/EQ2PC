# EQ2PC
Generation of root and child structures with identical 2-point correlation.

[![Actions Status](https://github.com/mauricio-fernandez-l/EQ2PC/workflows/Build/badge.svg)](https://github.com/mauricio-fernandez-l/EQ2PC/actions)

Contact:
* mauricio.fernandez.lb@gmail.com
* fritzen@mechbau.uni-stuttgart.de

Research group website: https://www.mib.uni-stuttgart.de/en/emma/

## Related work
Manuscript

"On the generation of periodic discrete structures with identical two-point correlation"

by Mauricio Fernández and Felix Fritzen.

## Description

The present repository offers a Python 3.7 implementation for the generation of discrete periodic structures with identical 2-point correlation as described in the work cited above. The source files are contained in the `src` folder, where `eq2pc.py` is the main module. The functionalities of the routines are demonstrated in the notebooks:

* [Demo 1: Search for 2PC-equivalent root structures](demo1_root_structures.ipynb)
* [Demo 2: Generation of 2PC-equivalent child structures](demo2_child_structures.ipynb)

Additionally, routines for the computation of the Voigt, Reuss and HS bounds for linear conductivity and elasticity are provided for given structures. This is demonstrated in:

* [Demo 3: Computation of bounds for conductivity and elasticity](demo3_bounds.ipynb)

Finally, limitations of the work of 

* Niezgoda, S.R., Fullwood, D.T., and Kalidindi, S.R. (2008): Delineation of the space of 2-point correlations in a composite material system. Acta Materialia 56, 5285 - 5292 

with respect to the determination of 2PC are discussed in the related work and demonstrated in the notebook

* [Demo 4: Limitations of Niezgoda et al (2008)](demo4_niezgoda_2008.ipynb)

**For notebook interaction:** please either download the repository or open it in 
[![Binder](https://mybinder.org/badge_logo.svg)](https://mybinder.org/v2/gh/mauricio-fernandez-l/EQ2PC/master)