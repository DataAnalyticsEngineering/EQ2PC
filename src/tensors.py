#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
@author: Mauricio Fernández

Module with auxiliary tensor routines.
"""

import numpy as np

# %% Utilities


def av(c, A):
    '''
    Average with phase concentrations c = [c1,c2,...] and matrices A = [A1,A2,...]
    '''
    return np.sum(np.transpose(A, axes=(1, 2, 0)) * c, axis=-1)

# %% Algebra


def sp(a, b):
    '''
    sp(a,b) returns the scalar product (full contraction) of tensors a and b of identical
    tensor.
    '''
    return np.sum(a * b)
