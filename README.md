# EQ2PC
Generation of root and child periodic discrete structures with identical 2-point correlation.

[![Actions Status](https://github.com/mauricio-fernandez-l/EQ2PC/workflows/Build/badge.svg)](https://github.com/mauricio-fernandez-l/EQ2PC/actions)

Contact:
* Mauricio Fernández: [email](mailto:fernandez@cps.tu-darmstadt.de), [website](https://www.maschinenbau.tu-darmstadt.de/cps/department_cps/team_1/team_detail_184000.en.jsp)
* Felix Fritzen: [email](mailto:fritzen@mechbau.uni-stuttgart.de), [website](https://www.simtech.uni-stuttgart.de/detail/mitarbeiter/Fritzen-00001/)

Research group website: http://www.mib.uni-stuttgart.de/dae

Last update: 2020-06-17

## Related work
Manuscript

"On the generation of periodic discrete structures with identical two-point correlation"

by Mauricio Fernández and Felix Fritzen. Submitted to *Proceedings of the Royal Society A* (under review). Pre-print available at https://arxiv.org/abs/2002.01234 .

## Description

The present repository offers a Python 3.7 implementation for the generation of discrete periodic structures with identical 2-point correlation as described in the work cited above. The source files are contained in the `src` folder, where `eq2pc.py` is the main module. The functionalities of the routines are demonstrated in the notebooks:

* [Demo 1: Search for 2PC-equivalent root structures](demo1_root_structures.ipynb)
* [Demo 2: Generation of 2PC-equivalent child structures](demo2_child_structures.ipynb)

Additionally, routines for the computation of the Voigt, Reuss and HS bounds for linear conductivity and elasticity are provided for given structures. This is demonstrated in:

* [Demo 3: Computation of bounds for conductivity and elasticity](demo3_bounds.ipynb)

Finally, limitations of the work of 

* Niezgoda, S.R., Fullwood, D.T., and Kalidindi, S.R. (2008): Delineation of the space of 2-point correlations in a composite material system. Acta Materialia 56, 5285 - 5292 

with respect to the determination of 2PC are discussed in the related work and demonstrated in the notebook

* [Demo 4: Formal limitations of Niezgoda et al (2008)](demo4_niezgoda_2008.ipynb)

**For notebook interaction:** please either download the repository or open it in 
[![Binder](https://mybinder.org/badge_logo.svg)](https://mybinder.org/v2/gh/mauricio-fernandez-l/EQ2PC/master)