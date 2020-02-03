import sympy as sym 
from sympy.matrices import matrix_multiply_elementwise

#%% 

print('*****************************************************')
print('Example illustrating issue in Niezgoda et al. (2008)')
print('*****************************************************')

#%%

def rIs(S):
    n = max(S)
    Is = [sym.zeros(*S.shape) for a in range(n)]
    for a in range(n):
        for p1 in range(S.shape[0]):
            for p2 in range(S.shape[1]):
                if S[p1,p2]==a+1: Is[a][p1,p2]=1
    return Is

def dft(a):
    P = a.shape
    F = [
         sym.Matrix([[sym.exp(-sym.I*2*sym.pi*p*q/PP) for q in range(PP)] for p in range(PP)])
         for PP in P
         ]
    return F[0]*a*sym.transpose(F[1])

def idft(a):
    P = a.shape
    F = [
         sym.Matrix([[sym.exp(sym.I*2*sym.pi*p*q/PP)/PP for q in range(PP)] for p in range(PP)])
         for PP in P
         ]
    return F[0]*a*sym.transpose(F[1])

def conv(a,b):
    af = dft(sym.conjugate(a))
    bf = dft(b)
    return idft(matrix_multiply_elementwise(sym.conjugate(af),bf)).expand().simplify()

def rCs(S):
    Is = rIs(S)
    n = len(Is)
    return [[sym.re(conv(Is[a],Is[b]).expand().simplify()) for b in range(n)] for a in range(n)]

#%% Check
    
def reverse(a):
    P = a.shape
    out = sym.zeros(*P)
    for p1 in range(P[0]):
        for p2 in range(P[1]):
            out[p1,p2] = a[-p1,-p2]
    return out

def checkDFTreal(Csf):
    out = []
    n = len(Csf)
    P = Csf[0][0].shape
    z = sym.zeros(*P)
    for a in range(n):
        for b in range(n):
            temp = Csf[a][b]-sym.conjugate(reverse(Csf[a][b]))
            check = temp.expand().simplify().as_real_imag()
            out.append(check[0]==z and check[1]==z)
    return all(out)

def check2PCdef(Csf):
    out = []
    n = len(Csf)
    P = Csf[0][0].shape
    z = sym.zeros(*P)
    for a in range(n):
        for b in range(n):
            temp = Csf[a][b]-sym.conjugate(Csf[b][a])
            check = temp.expand().simplify().as_real_imag()
            out.append(check[0]==z and check[1]==z)
    return all(out)

def check2PCabc(Csf):
    out = []
    n = len(Csf)
    P = Csf[0][0].shape
    z = sym.zeros(*P)
    for a in range(n):
        for b in range(n):
            for c in range(n):
                temp = matrix_multiply_elementwise(Csf[a][c],Csf[c][b]) - matrix_multiply_elementwise(Csf[c][c],Csf[a][b])
                check = temp.expand().simplify().as_real_imag()
                out.append(check[0]==z and check[1]==z)
    return all(out)

def check2PClast(Csf):
    out = []
    n = len(Csf)
    P = Csf[0][0].shape
    z = sym.zeros(*P)
    const = sym.zeros(*P)
    const[0,0] = 1
    for a in range(n):
        temp = sym.zeros(*P)
        for c in range(n):
            temp = temp + Csf[a][c]
        temp = temp - P[0]*P[1]*sym.sqrt(sym.re(Csf[a][a][0,0]))*const
        check = temp.expand().simplify().as_real_imag()
        out.append(check[0]==z and check[1]==z)
    for b in range(n):
        temp = sym.zeros(*P)
        for c in range(n):
            temp = temp + Csf[c][b]
        temp = temp - P[0]*P[1]*sym.sqrt(sym.re(Csf[b][b][0,0]))*const
        check = temp.expand().simplify().as_real_imag()
        out.append(check[0]==z and check[1]==z)
    return all(out)

def checkDFT0(Csf):
    out = []
    n = len(Csf)
    P = Csf[0][0].shape
    for a in range(n):
        for b in range(n):
            temp = 0
            for p1 in range(P[0]):
                for p2 in range(P[1]):
                    temp = temp + Csf[a][b][p1,p2]
            if a==b:
                temp = temp - P[0]*P[1]*sym.sqrt(sym.re(Csf[a][a][0,0]))
            check = temp.expand().simplify().as_real_imag()
            out.append(check[0]==0 and check[1]==0)
    return all(out)

def checkBounds(Csf):
    out = []
    n = len(Csf)
    P = Csf[0][0].shape
    for a in range(n):
        for b in range(n):
            for p1 in range(P[0]):
                for p2 in range(P[1]):
                    temp1 = Csf[a][b][p1,p2].expand().simplify()
                    temp1 = temp1*sym.conjugate(temp1)
                    temp1 = sym.re(temp1.expand().simplify())
                    temp1 = sym.sqrt(temp1)
                    temp2 = sym.re(Csf[a][b][0,0].expand().simplify())
                    check = temp1<=temp2 and temp2<=(P[0]*P[1])**2
                    out.append(check)
    return all(out)

def checkAll(Csf):
    return all([
        checkDFTreal(Csf)
        ,checkDFT0(Csf)
        ,check2PCabc(Csf)
        ,check2PCdef(Csf)
        ,check2PClast(Csf)
        ,checkBounds(Csf)
        ])

#%% Counterexample for Niezgoda et al. (2008)
    
# Structure
S = sym.Matrix([[2,3,1],[2,2,1],[2,3,3],[2,2,3]])

# Extract number of phases, periods
n = max(S)
P = S.shape

# Compute 2PC and DFT
Cs = rCs(S)
Csf = [[dft(Cs[a][b]).expand().simplify() for b in range(n)] for a in range(n)]

# Display structures and independent DFTs
print('\nStructure')
print(S)

print('\nDFT of 2PC - split into real and imaginary parts: (alpha,beta) = (1,1)')
temp = Csf[0][0].as_real_imag() 
print('real')
print(temp[0])
print('imaginary')
print(temp[1])

print('\nDFT of 2PC - split into real and imaginary parts: (alpha,beta) = (1,2)')
temp = Csf[0][1].as_real_imag() 
print('real')
print(temp[0])
print('imaginary')
print(temp[1])

print('\nDFT of 2PC - split into real and imaginary parts: (alpha,beta) = (2,2)')
temp = Csf[1][1].as_real_imag() 
print('real')
print(temp[0])
print('imaginary')
print(temp[1])

# Check all properties
print('\nCheck all properties for computed DFT of 2PC')
print(checkAll(Csf))

#%% Perturbate DFT of 2PC for choice gamma = 1

# Perturbation
z = sym.zeros(*P)
Delta = [[z,z,z],[z,z,z],[z,z,z]]
Delta[1][1] = sym.Matrix([[0,0,0],[0,0,0],[-2,1,1],[0,0,0]])
Delta[2][2] = Delta[1][1]
Delta[1][2] = -Delta[1][1]
Delta[2][1] = Delta[1][2]

# Perturbed DFT of 2PC
Csf_Delta = [[Csf[a][b] + Delta[a][b] for b in range(n)] for a in range(n)]

# Check that all DFT for choice gamma=1 are identical
print('\nCheck that all DFT for choice gamma=1 are identical (+Delta)')
check = all([Csf[0][b]==Csf_Delta[0][b] and Csf[b][0]==Csf_Delta[b][0] for b in range(n)])
print(check)

# Check all properties
print('\nCheck all properties for perturbed DFT of 2PC (+Delta)')
print(checkAll(Csf_Delta))