#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
@author: Mauricio Fernández

Module with routines for elastictiy problem.
"""

# Computation
import numpy as np

# Plotting
import matplotlib.pyplot as plt

# My stuff
from . import tensors as ten

# %% Overall convention
'''
Normalized Voigt notation, i.e., 6x6 matrices with sqrt(2) blocks for indices
(1,1), (2,2), (3,3), (2,3), (1,3), (1,2).
'''

# %% Isotropic projectors

P1 = np.zeros([6, 6])
P1[:3, :3] = 1 / 3

IS = np.eye(6)

P2 = IS - P1

# %% Stiffnesses


def iso(k):
    '''
    Return isotropic stiffness based on given eigenvalues.
    '''
    return k[0] * P1 + k[1] * P2


def iso_ev(K):
    '''
    Extract eigenvalues of given isotropic stiffness.
    '''
    return np.array([
        ten.sp(K, P1) / ten.sp(P1, P1), ten.sp(K, P2) / ten.sp(P2, P2)
    ])

# %% Bounds


def voigt(v, K):
    '''
    Voigt bound.
    '''
    return ten.av(v, K)


def reuss(v, K):
    '''
    Reuss average.
    '''
    return np.linalg.inv(ten.av(v, np.linalg.inv(K)))


def rP0(K0):
    '''
    Polarization tensor from Willis 1977 for hs in 3-dimensional elasticity.
    (Willis 1977, (5.11), i.e., isotropic 2-point statistics, no long-range order, isotropic C0,...)
    '''
    k = iso_ev(K0)
    p = [1 / (k[0] + 2 * k[1]), 2 / (5 * k[1]) *
         (k[0] + 3 * k[1]) / (k[0] + 2 * k[1])]
    return iso(p)


def hs(v, K, K0):
    '''
    Hashin-Shtrikman bound for comparison medium C0 for anisotropic 3d case.
    (Willis 1977, i.e., isotropic 2-point statistics, no long-range order, isotropic C0,...)
    (Algebraic form : Walpole 1966, see Lobos Fernández and Böhlke 2018)
    CHS = C0 - P0^{-1} + <L>^{-1}
    L = [(C-C0) + P0^{-1}]^{-1}
    '''
    P0 = rP0(K0)
    L = [np.linalg.inv(KK - K0 + np.linalg.inv(P0)) for KK in K]
    Lav = ten.av(v, L)
    return K0 - np.linalg.inv(P0) + np.linalg.inv(Lav)

# %% Evaluate structure


def S_iso_bounds(S, K):
    '''
    Evaluate structure for isotropic stiffnesses.
    '''
    nph = np.max(S)
    ntotal = np.prod(S.shape)
    v = np.array([np.sum(S == i) / ntotal for i in range(1, nph + 1)])

    # Get optimal isotropic comparison media
    # (for generalization to anisotropic stiffness list K, use Lobos and Böhlke (2016))
    evs = np.array([iso_ev(KK) for KK in K])
    K0low = iso(np.min(evs, axis=0))
    K0upp = iso(np.max(evs, axis=0))

    # Copute HS bounds
    hslow = hs(v, K, K0low)
    hsupp = hs(v, K, K0upp)

    # List of bounds
    out = np.array([
        reuss(v, K), hslow, hsupp, voigt(v, K)
    ])

    # Plot for 2 phases
    if nph == 2:
        plot_2ph_iso(K, vl=v[0])

    # Return bounds
    return out

# %% Plots


def plot_2ph_iso(K, n=20, vl=0):
    '''
    Plot for 2-phase isotropic 3d case.
    '''
    v1 = np.linspace(0, 1, n)

    vb = np.array([
        voigt([v, 1 - v], K)
        for v in v1
    ])
    rb = np.array([
        reuss([v, 1 - v], K)
        for v in v1
    ])

    evs = np.array([iso_ev(KK) for KK in K])
    K0low = iso(np.min(evs, axis=0))
    K0upp = iso(np.max(evs, axis=0))

    hslow = np.array([
        hs([v, 1 - v], K, K0low)
        for v in v1
    ])
    hsupp = np.array([
        hs([v, 1 - v], K, K0upp)
        for v in v1
    ])

    for i in [0, 3]:
        plt.figure()
        plt.plot(v1, vb[:, i, i], color='blue', linestyle='-', label='Voigt', linewidth=2)
        plt.plot(v1, hsupp[:, i, i], color='lightgray',
                 linestyle='-', label='HS_upp', linewidth=2)
        plt.plot(v1, hslow[:, i, i], color='lightgray',
                 linestyle='--', label='HS_low', linewidth=2)
        plt.plot(v1, rb[:, i, i], color='blue', linestyle='--', label='Reuss', linewidth=2)
        if vl != 0:
            plt.axvline(vl, color='gray', linestyle='--', linewidth=2)
        plt.legend()
        plt.xlabel('$v_1$')
        plt.ylabel('Stiffness component $K_{%i%i}$' % (i + 1, i + 1))
        plt.show()
