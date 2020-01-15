import time
import random as r
from parametres import *

######Fonctions utilitaires######

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

def prodMatVec(M,v): #OK normalement (un test)
    res=[0 for m in M]
    for i in range(len(M)):
        for j in range(len(v)):
            res[i] += M[i][j] * v[j]
    return res

def prodMatMat(M1,M2): 
    #Eventuelement pour faire des tests, pas utilisée dans les calculs
    res=[[0 for m2 in M2[0]] for m1 in M1]
    for i in range(len(M1)):
        for j in range(len(M2[0])):
            for k in range(len(M2[0])):
                res[i][j] += M1[i][k] * M2[k][j]
    return res

###### Redéfinition du PCG_128 (avec C aléatoire) ######

def sortiesGenerateur():#OK !
    c = (r.randint(0, 2**(2 * k)) * 2 + 1) % 2**(2 * k) #c est impair
    S=[r.randint(0,2**(k * 2))]
    #c=6364136223846793005*2**64+1442695040888963407#increment par defaut de pcg (connu)
    #S=[8487854484825256858 + 11929896274893053136 * 2**64]
    for i in range (nboutput-1):
        S.append((S[i]* a + c) % 2**(2 * k))
    X=[]
    for i in range (nboutput):
        x=(S[i] % 2**k)^(S[i]>>k)
        rot = S[i] >> (2 * k - known_up)
        X.append((x >> rot) | ((x << (k - rot)) % 2**k))
    return X,S,c

###### Rotations ######
def rotateX(X, rot):#OK !
    rX = [];
    for i in range(nbiter):
        rX.append(((X[i] // 2**rot[i]) + ((X[i] * 2**(k - rot[i])) % 2**(k))))
    return rX

def unrotateX(X, rot):#OK !
    rot2 = []
    for i in range(nbiter):
        rot2.append((k - rot[i]) % k)
    return rotateX(X, rot2)

def unrotate1(Xi):#OK !
    return (Xi >> (k - 1)) | ((Xi << 1) % 2**k)

###### Sous-fonctions de FindDS######
'''
polW = powA * W0
polC = polA * c (idem pour littleC)
'''

def getY(W0, WC, rot, uX):#OK !
    Y = [(((powA[i] * W0 + polA[i] * WC) % 2**known_low) ^ (uX[i] % 2**known_low)) * 2**known_up + (rot[i] ^ (uX[i] // 2**(k - known_up))) for i in range(nbiter)]
    return Y

def getYprim(Y, WC, W0): #OK ! avec erreurs de retenues ~64bits (polC polW)
    Yprim=[(Y[i] - (polA[i] * WC + powA[i]*W0) // 2**(k - known_up)) % 2**(known_up + known_low) for i in range(nbiter)]
    return Yprim

def getDY(Y, WC, W0): #OK ! avec erreurs de retenues ~64bits (polC polW)
    Yprim = getYprim(Y, WC, W0)
    DY=[(Yprim[i+1] - Yprim[i]) % 2**(known_up + known_low) for i in range(nbiter-1)]
    return DY

####### UNIQUEMENT POUR LES TESTS ######
def getDS(polW,S):
    DS=[((S[i+1] - S[i] - polC[i+1] - polW[i+1] + polC[i] + polW[i]) % 2**(k//2 + connus_low))// 2**connus_low for i in range(nbiter-1)]
    return DS

######FINDDS######
def FindDS64(uX, rot, W0,WC): #rajouter rot dans la version non test ? #OK! ~64bits
    #polW = getPolW(W0)
    Y = getY(W0, WC, rot, uX)
    DY = getDY(Y, WC, W0) #OK avec erreurs de retenues!
    tmp = [y * 2**(k - known_up - known_low) for y in DY]#on rajoute les zéros, recentrage impossible à cause des erreurs de retenues
    u = prodMatVec(invG, tmp)
    DS64 = prodMatVec(Greduite, [round(u_) for u_ in u])
    return DS64

######FINDROTI######
#DS64ij = ((polA[j] - polA[i])*DSmod0) % 2**k

def FindRoti(DS640, X, i, Y0, W0,WC):#OK !
    DS640i = (polA[i] * DS640) % 2**k
    DSmod0i = ((DS640i << known_low) + W0 * powA[i] + WC * polA[i] - WC - W0) % 2**(k +known_low)
    print("DS640i")
    print(DS640i)
    # Yi = vraiYi ou vraiYi - 1 à cause de la retenue
    Yi1 = (Y0 + (DSmod0i >> (k - known_up))) % (1 << (known_low + known_up))#avec ou sans retenue
    Yi2 = Yi1 + 1
    Wi = (W0 * powA[i] + WC * polA[i]) % (1 << known_low)
    roti = []
    for i in range(k):
        test1 = (((X ^ (Yi1 >> known_up)) % (1 << known_low)) == Wi) and ((i ^ (X >> (k - known_up))) == Yi1 % (1 << known_up))
        test2 = (((X ^ (Yi2 >> known_up)) % (1 << known_low)) == Wi) and ((i ^ (X >> (k - known_up))) == Yi2 % (1 << known_up))
        
        if (test1 or test2) :
            roti.append(i)
        X = unrotate1(X)
    return roti

def FindRot(DS640,X, Y0, W0, WC, n): #OK !
    rot =[]
    for i in range(n):
        rot.append(FindRoti(DS640, X[i], i, Y0, W0,WC))
        #print(rot[i])
        if(len(rot[i]) == 0):
            return []
    return rot

'''def getUpDS(DSmod0, rot,n):
    DY = []
    for i in range(n):
        DY.append(rot[i][0] - (((powA[i] * DSmod0) % (1<< (2 * k))) >> (2 * k - known_up))'''
