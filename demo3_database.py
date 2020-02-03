import numpy as np
import nexusformat.nexus as nx

import source.eq2pc as e2
import source.database as db

#%% Create hdf5 database from examples folder

print('**************************************************************')
print('Create hdf5 database from examples folder')
print('**************************************************************')

# Create database
db.extodb()

# Print structure
h = nx.nxload('data/database.hdf5')
print(h.tree)
    
#%% Branch in database

print('**************************************************************')
print('Branch in database')
print('**************************************************************')

# Create branch step by step
# First level
Ks = [
      np.array([[1,1,1],[0,1,0]])
      ,np.array([[1,0,0],[1,1,0]])
      ]
phc = [[2,1]]
db.branch('ex/2d/2ph1',name='b1',Ks=Ks,phc=phc,info=True)
# Second level
Ks = [
      np.array([[1,1,1],[0,1,0]])
      ,np.array([[1,0,0],[1,1,0]])
      ]
phc = []
db.branch('ex/2d/2ph1/b1',name='b1',Ks=Ks,phc=phc)

# Create another branch from root
Ks = [
      np.array([[0,0,0,0,0],[0,0,0,0,0],[0,1,1,1,0],[0,0,0,0,0],[0,0,0,0,0]])
      ,np.array([[0,0,0,0,0],[0,1,0,0,0],[0,0,1,0,0],[0,0,0,1,0],[0,0,0,0,0]])
      ]
phc = [[1,2]]
db.branch('ex/2d/2ph1',name='fibers',Ks=Ks,phc=phc)

# Immediatlely create branch of arbitrary depth
db.branch('ex/2d/2ph1/b2/b1/b1/b1')

# Print structure
h = nx.nxload('data/database.hdf5')
print(h.tree)

#%% Evaluate branches

print('**************************************************************')
print('Evaluate branches')
print('**************************************************************')

# Execute a branch, leaves can be created
out = db.branchev('ex/2d/2ph1/b1/b1',leaves=True,plot=True)
print(out)

# Print structure
h = nx.nxload('data/database.hdf5')
print(h.tree)

# Execute another branch and get generated structure
out = db.branchev('ex/2d/2ph1/fibers',leaves=True,plot=True)

#%% Get structures from database

print('**************************************************************')
print('Get structures from database')
print('**************************************************************')

Ss = db.getstructures('ex/2d/2ph1/fibers')
if len(Ss)>0:
    e2.plotS(Ss[0])
    e2.plotS(Ss[1])