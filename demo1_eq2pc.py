import source.eq2pc as e2

#%% Search for initial 2PC-equivalent structures

#%% 1d 

print('*****************************************************')
print('1d')
print('*****************************************************')

# 2 phases
out = e2.Ceqsearch([12],[5],info=True)
if len(out)>0:
    S1 = out[0][0]
    S2 = out[0][1]
    print('Are the structures related?\n\t%r' % e2.rel(S1,S2))
    print('Are the structures 2PC-equivalent?\n\t%r' % e2.Ceq(S1,S2))
    e2.plotS(S1)
    e2.plotS(S2)

# 3 phases   
out = e2.Ceqsearch([12],[2,3],info=False)
if len(out)>0:
    S1 = out[0][0]
    S2 = out[0][1]
    print('Are the structures related?\n\t%r' % e2.rel(S1,S2))
    print('Are the structures 2PC-equivalent?\n\t%r' % e2.Ceq(S1,S2))
    e2.plotS(S1)
    e2.plotS(S2)
    
#%% 2d
    
print('*****************************************************')
print('2d')
print('*****************************************************')
  
out = e2.Ceqsearch([4,3],[5],info=False)
if len(out)>0:
    S1 = out[0][0]
    S2 = out[0][1]
    print('Are the structures related?\n\t%r' % e2.rel(S1,S2))
    print('Are the structures 2PC-equivalent?\n\t%r' % e2.Ceq(S1,S2))
    e2.plotS(S1)
    e2.plotS(S2)
    
out = e2.Ceqsearch([4,3],[3,2],info=False)
if len(out)>0:
    S1 = out[0][0]
    S2 = out[0][1]
    print('Are the structures related?\n\t%r' % e2.rel(S1,S2))
    print('Are the structures 2PC-equivalent?\n\t%r' % e2.Ceq(S1,S2))
    e2.plotS(S1)
    e2.plotS(S2)
    
#%% 3d
    
print('*****************************************************')
print('3d')
print('*****************************************************')
    
out = e2.Ceqsearch([2,3,4],[5],info=False)
if len(out)>0:
    S1 = out[0][0]
    S2 = out[0][1]
    print('Are the structures related?\n\t%r' % e2.rel(S1,S2))
    print('Are the structures 2PC-equivalent?\n\t%r' % e2.Ceq(S1,S2))
    e2.plotS(S1)
    e2.plotS(S2)
    
#%% Save to examples folder and load from it

print('*****************************************************')
print('Save to examples folder and load from it')
print('*****************************************************')

# Save
folder = e2.Ceqsave([4,3],[5],info=False)

# Load
out = e2.Ceqload(folder)
print(out)

# Print current examples
out = e2.Ceqexamples()

# Generate and save more examples
e2.Ceqsave([12],[5],info=False)
e2.Ceqsave([4,3],[3,2],info=False)
e2.Ceqsave([2,3,4],[5],info=False)

# Print current examples
out = e2.Ceqexamples()