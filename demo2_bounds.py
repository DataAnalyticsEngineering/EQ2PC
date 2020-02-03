import numpy as np

import source.eq2pc as e2
import source.conductivity as con
import source.elasticity as el

#%% Conductivity

print('*********************************************************')
print('Bounds in conductivity')
print('*********************************************************')

K = [10,1]
con.plot2phiso(K,d=2)

#%% Elasticity

print('*********************************************************')
print('Bounds in elasticity')
print('*********************************************************')

K = np.array([el.iso([10,2]),el.iso([1,5])])
el.plot2phiso(K)

#%% Load an example

print('*********************************************************')
print('Available examples')
print('*********************************************************')

ex = e2.Ceqexamples()

S1,S2 = e2.Ceqload(ex[2])

#%% Evaluation of structures

print('*********************************************************')
print('Load examples and evaluate bounds')
print('*********************************************************')

# 2d conductivity, 2 phases, isotropic
print('-------------------------------------------------')
print('2d conductivity, 2 phases, isotropic')
S1,S2 = e2.Ceqload(ex[2])
K = np.array([10,1])
print(con.Siso(S1,K))

# 3d conductivity, 2 phases, isotropic
print('-------------------------------------------------')
print('3d conductivity, 2 phases, isotropic')
S1,S2 = e2.Ceqload(ex[3])
K = np.array([10,1])
print(con.Siso(S1,K))

# 3d elasticity, 2 phases, isotropic
print('-------------------------------------------------')
print('3d elasticity, 2 phases, isotropic')
S1,S2 = e2.Ceqload(ex[3])
K = np.array([el.iso([1,5]),el.iso([10,2])])
print(el.Siso(S1,K))
