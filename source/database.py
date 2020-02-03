#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
@author: Mauricio Fernández

Module with routines for the generation of an hdf5 database for 2PC-equivalent structures.
"""

# Modules
import numpy as np
import h5py as h5

# My stuff
from . import eq2pc as e2

#%% Migrate examples in normal folder to hdf5 database

def extodb():
    '''
    Create hdf5 database from examples.
    '''
    ex = e2.Ceqexamples()
    
    h = h5.File('data/database.hdf5','w')
    h.create_group('ex/1d')
    h.create_group('ex/2d')
    h.create_group('ex/3d')
    
    for e in ex:
        Ss = e2.Ceqload(e,info=False)
        d = np.ndim(Ss[0])
        nev = e2.rnev(Ss[0])
        keys = h['ex/%id/' % (d)].keys()
        nkeys = len([k for k in keys if '%iph' % len(nev) in k])
        g = h.create_group('ex/%id/%iph%i' % (d,len(nev),nkeys+1))
        for i in range(len(Ss)): g.create_dataset('root%i' % (i+1),data=Ss[i])
    
    h.close()
    
#%% Manipulate hdf5 database
    
def branch(folder,name='b1',Ks=[],phc=[],enforce=False,info=False):
    '''
    Create a branch with kernels and phase concatenation in database.
    '''
    # Access database
    h = h5.File('data/database.hdf5','a')
    
    # Check: if folder does not exists, create
    if not(folder in h):
        if info:
            print('...Folder %s does not exist. Create folder.' % folder)
        h.create_group(folder)
    
    # Branch name
    subfolder = folder+'/'+name
    
    # If branch already exist, print note.
    if subfolder in h and not enforce:
        print('...Branch %s ALREADY exists.' % subfolder)
        print('...Listing existing content.')
        content = h[subfolder].keys()
        print('\tNumber of elements in branch: %i' % len(content))
        for c in content: print('\t%s' % c)
        print('-> Use enforce=True for enforced branch creation.\n')
    
    # Enforce, i.e., delete old existing branch
    if subfolder in h and enforce:
        h.__delitem__(subfolder)
    
    # Branch creation
    if subfolder not in h:
        if info:
            print('...Branch %s does not exist, filling.' % subfolder)
        g = h.create_group(subfolder)
    
        # Write kernels
        for i in range(len(Ks)):
            g.create_dataset('k%i' % (i+1),data=Ks[i])
            
        # Write phase concatenation
        if len(phc)>0:
            temp = np.sum([len(l) for l in phc])
            mat = np.zeros([temp,2],dtype=np.int)
            counter = 0
            for i1 in range(len(phc)):
                for i2 in range(len(phc[i1])):
                    mat[counter] = [i1+1,phc[i1][i2]]
                    counter += 1
            g.create_dataset('phc',data=mat)
            
        # Written content
        if info:
            print('...New content:')
            content = h[subfolder].keys()
            print('\tNumber of elements in branch: %i' % len(content))
            for c in content: print('\t%s' % c)
        
    h.close()

def branchev(branch,leaves=False,plot=False):
    '''
    Evaluate branch.
    '''
    # Access database
    h = h5.File('data/database.hdf5','a')
    
    # Split path
    l = branch.split('/')
    
    # Find root folder
    rootfolder = ''
    for ll in l[:3]:
        rootfolder += '/'+ll
    
    # Get roots
    roots = []
    for r in h[rootfolder].keys():
        if 'root' in r:
            roots.append(h[rootfolder+'/'+r][...])
            if plot:
                e2.plotS(roots[-1])
    
    # Evaluate every level
    depth = len(l[3:])
    for i in range(1,depth+1):
        lev = ''
        for ll in l[:3+i]:
            lev += '/'+ll
        Ks = [h[lev+'/'+k][...] for k in h[lev].keys() if 'k' in k]
        for r in range(len(roots)):
            if len(Ks)>0:        
                roots[r] = e2.se(roots[r],Ks=Ks)
            else:
                roots[r] = e2.se(roots[r])
            if 'phc' in h[lev].keys():
                mat = h[lev+'/phc'][...]
                phc = [list(mat[mat[:,0]==i,1]) for i in range(1,np.max(mat[:,0])+1)]
                roots[r] = e2.phc(roots[r],phc)
            if plot:
                e2.plotS(roots[r])
            if leaves:
                if (lev+'/leaf%i' % (r+1)) in h:
                    h.__delitem__(lev+'/leaf%i' % (r+1))
                h[lev].create_dataset('leaf%i' % (r+1),data=roots[r])
            
    h.close()
    
    return roots

#%% Get items from database

def getstructures(branch):
    '''
    Get leaves in branch.
    '''
    # Access database
    h = h5.File('data/database.hdf5','a')
    
    # Get structures (root or leaves) if branch exists
    if branch not in h:
        print('...Branch does NOT exist in database')
        h.close()
    else:
        if 'root1' in list(h[branch].keys()):
            out = np.array([
                    h[branch+'/root1']
                    ,h[branch+'/root2']
                    ])
        else:
            out = np.array([
                    h[branch+'/leaf1']
                    ,h[branch+'/leaf2']
                    ])
        h.close()
        return out