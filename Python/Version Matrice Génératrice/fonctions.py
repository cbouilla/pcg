import time
import random as r

########Fonctions utiles###########

#Pour afficher en binaire

def dec2bin(d,nb=8):
    """Représentation d'un nombre entier en chaine binaire (nb: nombre de bits du mot)"""
    if d == 0:
        return "0".zfill(nb)
    if d<0:
        d += 1<<nb
    b=""
    while d != 0:
        d, r = divmod(d, 2)
        b = "01"[r] + b
    return b.zfill(nb)

###Calcul modulo (pour méthode LU)###
def euclide(a, b):
    x,y, u,v = 0,1, 1,0
    while a != 0:
        q,r = b//a,b%a; m,n = x-u*q,y-v*q # use x//y for floor "floor division"
        b,a, x,y, u,v = a,r, u,v, m,n
    return b, x, y

def modinv(a, mod): #mod premier
    g, x, y = euclide(a, mod) 
    return x % mod

def FactorisationLUMod(M,mod): #mod premier
    n=len(M)
    L=[[0 for i in range(n)]for i in range(n)]
    U=[[M[i][j] % mod for j in range(n)]for i in range(n)]
    for i in range(n):
        L[i][i]=1
        inv=modinv(U[i][i],mod)
        for j in range(i+1,n):
            L[j][i]=(U[j][i] * inv) % mod
            for k in range(n):
                U[j][k]=(U[j][k] - L[j][i] * U[i][k]) % mod
    return L,U

def RemonteeLUMod(L,U,b,mod): #mod premier
    n=len(b)
    # L*x=b % mod avec x=U*a, on cherche x
    x=b
    for i in range(n):
        for j in range(i):
            b[i]=(b[i] - b[j] * L[i][j]) % mod
    
    # U*a=x % mod on cherche a
    a=x
    for i in range(1,n+1):
        for j in range(1,i):    
            a[n-i]=(a[n-i] - a[n-j] * U[n-i][n-j]) % mod
        a[n-i]=(a[n-i] * modinv(U[n-i][n-i],mod)) % mod
    # ici, vérifier si L*U*x == b
    return a

def prodMatVecMod(M,v,mod): #OK normalement (un test)
    res=[0 for m in M]
    for i in range(len(M)):
        for j in range(len(v)):
            res[i] = (res[i] + M[i][j] * v[j]) % mod
    return res

def prodMatVec(M,v): #OK normalement (un test)
    res=[0 for m in M]
    for i in range(len(M)):
        for j in range(len(v)):
            res[i] += M[i][j] * v[j]
    return res


def prodMatMatMod(M1,M2,mod):
    res=[[0 for m2 in M2[0]] for m1 in M1]
    for i in range(len(M1)):
        for j in range(len(M2[0])):
            for k in range(len(M2[0])):
                res[i][j] = (res[i][j] + M1[i][k] * M2[k][j]) % mod
    return res

def prodMatMat(M1,M2):
    res=[[0 for m2 in M2[0]] for m1 in M1]
    for i in range(len(M1)):
        for j in range(len(M2[0])):
            for k in range(len(M2[0])):
                res[i][j] += M1[i][k] * M2[k][j]
    return res

def TCRVec(X1,X2,p1,p2): #p1 et p2 premiers entre eux
    g, a, b=euclide(p1, p2) #a est l'inverse de x1 mod p2
    p, y1, y2 = p1 * p2, p2 * b, p1 * a
    X=[(X1[i] * y1 + X2[i] * y2) % p for i in range(len(X1))]
    return X
    
###Sorties de pcg sans rotation###
def sortiesGenerateur(k,a,c,nbiter):# ~128bits
    S=[r.randint(0,2**k)]
    for i in range (nbiter-1):
        S.append((S[i]* a + c) % 2**k)
    X=[(s % 2**(k//2))^(s>>(int(k//2))) for s in S]
    return X,S