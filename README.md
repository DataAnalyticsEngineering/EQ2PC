# EQ2PC
Generation of structures with equivalent 2-point correlation

Contact:
* mauricio.fernandez.lb@gmail.com
* fritzen@mechbau.uni-stuttgart.de

Research group website: https://www.mib.uni-stuttgart.de/en/emma/

## Related work
Manuscript

"On the generation of periodic discrete structures with identical two-point correlation"

by Mauricio Fernández and Felix Fritzen.

## Description

The present repository offers a Python 3 implementation for the generation of discrete periodic structures with identical 2-point correlation as described in the work cited above. The source files are contained in the 'src' folder, where 'eq2pc.py' is the main module. The functionalities of the routines are demonstrated in the notebooks:

* [Demo 1: Generation of root structures](demo1_root_structures.ipynb)
* [Demo 2: Generation of child structures](demo2_child_structures.ipynb)

Additionally, routines for the computation of the Voigt, Reuss and HS bounds for linear conductivity and elasticity are provided for given structures. This is demonstrated in:

* [Demo 3: Computation of bounds](demo3_bounds.ipynb)

Finally, limitations of the work of Niezgoda, Fullwood and Kalidindi (2008) with respect to the determination of 2PC are discussed in the related work and demonstrated in the notebooks

* [Demo 4: Limitations of NFK](demo4_niezgoda_2008.ipynb)