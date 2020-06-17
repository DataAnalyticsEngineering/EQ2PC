#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
@author: Mauricio Fernández

Module with routines for conductivity problem.
"""

# Computation
import numpy as np

# Plotting
import matplotlib.pyplot as plt

# My stuff
from . import tensors as ten

# %% Bounds


def voigt_iso(v, K):
    '''
    Voigt bound for isotropic phases with volume fractions v and conductivities
    K in dimension d in {2,3}. Returns a scalar.
    Example:
        v = [0.1,0.2,0.7]
        K = [10,2,5]
    '''
    return np.sum(np.array(v) * np.array(K))


def voigt(v, K):
    '''
    Voigt bound with volume fractions v and conductivities
    K for anisotorpic 3d problems. Returns a matrix.
    Example:
        v = [0.1,0.2,0.7]
        K = [np.array([[10,2,0],[2,5,0],[0,0,8]]),2*np.eye(3),5*np.eye(3)]
    '''
    return ten.av(v, K)


def reuss_iso(v, K):
    '''
    Reuss bound for isotropic phases.
    Returns a scalar.
    '''
    return 1 / np.sum(np.array(v) / np.array(K))


def reuss(v, K):
    '''
    Reuss for anisotropic phases.
    Returns a matrix.
    '''
    return np.linalg.inv(ten.av(v, np.linalg.inv(K)))


def hs_iso(v, K, d=2):
    '''
    Hashin-Shtrikman bounds for isotropic phases in dimension d.
    Torquato (2002), equation (21.23).
    '''
    Kmin = np.min(K)
    Kmins = (d - 1) * Kmin
    Kmax = np.max(K)
    Kmaxs = (d - 1) * Kmax
    low = 1 / np.sum(v * (1 / (K + Kmins))) - Kmins
    upp = 1 / np.sum(v * (1 / (K + Kmaxs))) - Kmaxs
    return np.array([low, upp])


def rP0(K0):
    '''
    Return polarization tensor from Willis 1977 for hs in 2/3-dimensional thermal conductivity.
    (Willis 1977, i.e., isotropic 2-point statistics, no long-range order, isotropic C0,...)
    '''
    k0 = K0[0, 0]
    d = len(K0)
    return np.eye(d) / (d * k0)


def hs(v, K, K0):
    '''
    Hashin-Shtrikman bound for comparison medium K0 for anisotropic 3d case.
    (Willis 1977, i.e., isotropic 2-point statistics, no long-range order, isotropic K0,...)
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
    Evaluate structure.
    '''
    nph = np.max(S)
    ntotal = np.prod(S.shape)
    d = np.ndim(S)
    print('...computing for dimension %i' % d)
    v = np.array([np.sum(S == i) / ntotal for i in range(1, nph + 1)])
    print('...volume fractions')
    print(v)
    out = np.array([
        reuss_iso(v, K), hs_iso(v, K, d=d)[0], hs_iso(
            v, K, d=d)[1], voigt_iso(v, K)
    ])
    if nph == 2:
        plot_2ph_iso(K, d=d, vl=v[0])
    return out

# %% Plots


def plot_2ph_iso(
        K,
        d=2,
        n=20,
        vl=0,
        x_label='$v_1$',
        plot_title='',
        point=[],
        point_label='',
        point2=[],
        point2_label='',
        save_path=''):
    '''
    Plot for 2-phase isotropic case in dimension d.
    '''
    v1 = np.linspace(0, 1, n)
    vb = np.array([
        voigt_iso([v, 1 - v], K)
        for v in v1
    ])
    rb = np.array([
        reuss_iso([v, 1 - v], K)
        for v in v1
    ])
    hs = np.array([
        hs_iso([v, 1 - v], K, d=d)
        for v in v1
    ])

    plt.figure()
    plt.plot(v1, vb, color='blue', linestyle='-', label='Voigt', linewidth=2)
    plt.plot(v1, hs[:, 1], color='lightgray', linestyle='-', label='HS_upper', linewidth=2)
    plt.plot(v1, hs[:, 0], color='lightgray', linestyle='--', label='HS_lower', linewidth=2)
    plt.plot(v1, rb, color='blue', linestyle='--', label='Reuss', linewidth=2)
    if vl != 0:
        plt.axvline(vl, color='gray', linestyle='--', linewidth=2)
    if len(point) > 0:
        plt.scatter([point[0]], [point[1]], label=point_label)
    if len(point2) > 0:
        plt.scatter([point2[0]], [point2[1]], label=point2_label)
    plt.legend()
    plt.title(plot_title)
    plt.xlabel(x_label)
    plt.ylabel('$K_{11,22}$')
    if len(save_path) > 0:
        plt.tight_layout()
        plt.savefig(save_path)
    plt.show()
