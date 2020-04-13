# coding: utf-8
import time
import random as r
import fpylll as f
import os
from itertools import product

k = 64
known_up  = 6
known_low  = 12
a = 2549297995355413924 * 2^64 + 4865540595714422341
nbiter = 5
nboutput = 40


polA = [0];
powA = [1];
for i in range (1, nboutput):
    polA.append((polA[i-1] + powA [i-1]) % 2^128)
    powA.append((powA[i-1] * a)% 2^128)

def dec2bin(d,nb= 2 * k - known_low + 1):
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

####   MATRIX   ####
def getG(n,mod):
    G = [[1]]
    for i in range(n-1):
        G[0].append(0)
    powA = a % mod
    for i in range(n-1):
        G.append([powA])
        powA = powA * a % mod
        for j in range(n-1):
            if(i == j):
                G[i+1].append(mod)
            else:
                G[i+1].append(0)
    return matrix(G)        

def getGreduite(n, mod):
    G = getG(n, mod)
    Gr = G.transpose().LLL()
    #Glll = f.IntegerMatrix.from_matrix(Gr)
    return Gr

Greduite1 = getGreduite(nbiter - 1, 2^k)

Greduite2 = getGreduite(nboutput - 1, 2^(2 * k - known_low))
print("Vecteur réputé court de taille : 2^{:.1f}".format(float(log(Greduite2[0].norm(), 2).n())))

#### Récupération des données ####
#def recupDonnees():
#   listfiles = os.listdir('results')
#   print(listfiles)
#   listDonnees = []
#   for i in range(listfiles):
#       listDonnees.append[]
#       f = open(listfiles[i], "r")
#       while(1):
#           s = f.readline()
            


def sortiesGenerateur():#OK !
    c = r.randrange(2**128) | 1
    S=[r.randint(0,1<<(k * 2))]
    #c=6364136223846793005 * 2^64 + 1442695040888963407#increment par defaut de pcg (connu)
    #S=[8487854484825256858 + 11929896274893053136 * 2^64]
    for i in range (nboutput-1):
        S.append((S[i] * a + c) % (1<<(2 * k)))
    X=[]
    for i in range(nboutput):
        x = (S[i] % (1<<k)) ^^ (S[i]>>k)
        rot = S[i] >> (2 * k - known_up)
        X.append((x >> rot) | ((x << (k - rot)) % (1 << k)))
    return X,S,c


## Unrotate
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
def getY(W0, WC, rot, uX):#OK !
    Y = [(((powA[i] * W0 + polA[i] * WC) % 2**known_low) ^^ (uX[i] % 2**known_low)) * 2**known_up + (rot[i] ^^ (uX[i] // 2**(k - known_up))) for i in range(nbiter)]
    return Y

def getYprim(Y, WC, W0): #OK ! avec erreurs de retenues ~64bits (polC polW)
    Yprim=[(Y[i] - (polA[i] * WC + powA[i]*W0) // 2**(k - known_up)) % 2**(known_up + known_low) for i in range(nbiter)]
    return Yprim

def getDY(Y, WC, W0): #OK ! avec erreurs de retenues ~64bits (polC polW)
    Yprim = getYprim(Y, WC, W0)
    DY=[(Yprim[i+1] - Yprim[i]) % 2**(known_up + known_low) for i in range(nbiter-1)]
    return DY

######FINDDS######
def FindDS64(uX, rot, W0, WC, Greduite): #rajouter rot dans la version non test ? #OK! ~64bits
    Y = getY(W0, WC, rot, uX)
    DY = getDY(Y, WC, W0) #OK avec erreurs de retenues!
    tmp = vector([y * 1<<(k - known_up - known_low) for y in DY])#on rajoute les zéros, recentrage impossible à cause des erreurs de retenues
    G = f.IntegerMatrix.from_matrix(Greduite)
    DS64 = f.CVP.closest_vector(G,tuple(tmp))
    return DS64, Y[0]

######FINDROTI######
#DS64ij = ((polA[j] - polA[i])*DSmod0) % 2**k

def FindRoti(DS640, X, i, Y0, W0,WC):#OK !
    DS640i = (polA[i] * DS640) % (1<<k)
    DSmod0i = ((DS640i << known_low) + W0 * powA[i] + WC * polA[i] - WC - W0) % (1<<(k +known_low))
    # Yi = vraiYi ou vraiYi - 1 à cause de la retenue
    Yi1 = (Y0 + (DSmod0i >> (k - known_up))) % (1 << (known_low + known_up))#avec ou sans retenue
    Yi2 = Yi1 + 1
    Wi = (W0 * powA[i] + WC * polA[i]) % (1 << known_low)
    roti = []
    for j in range(k):
        test1 = (((X ^^ (Yi1 >> known_up)) % (1 << known_low)) == Wi) and ((j ^^ (X >> (k - known_up))) == Yi1 % (1 << known_up))
        test2 = (((X ^^ (Yi2 >> known_up)) % (1 << known_low)) == Wi) and ((j ^^ (X >> (k - known_up))) == Yi2 % (1 << known_up))
        if (test1 or test2) :
            roti.append(j)
        X = unrotate1(X)
    return roti

def FindRot(DS640,X, Y0, W0, WC): #OK !
    tabrot =[]
    for i in range(nboutput):
        tabrot.append(FindRoti(DS640, X[i], i, Y0, W0,WC))
        if(len(tabrot[i]) == 0):
            return []
    return tabrot

def findDS(rot, W0, WC): #OK!
    rotprim = []
    for i in range(nboutput):
        rotprim.append((rot[i] - ((powA[i] * W0 + polA[i] * WC) >> (2 * k - known_up))) % (1<<known_up))
    tmp = vector([(rotprim[i+1] - rotprim[i]) << (2 * k - known_up - known_low) for i in range(nboutput - 1)])

    '''print("cheat_DS et tmp :")
    for i in range(nboutput - 1):
        print(dec2bin(cheat_DS[i]))
        print(dec2bin(tmp[i]))
        print((cheat_DS[i] - tmp[i])>>(2 * k - known_up - known_low))'''
    
    G = f.IntegerMatrix.from_matrix(Greduite2)
    return f.CVP.closest_vector(G,tuple(tmp))
    '''cvp = f.CVP.closest_vector(G,tuple(shifted_tmp))
    return vector([((cvp[i] << k) + DS64[0] * powA[i]) % (1 << (128 - known_low)) for i in range(nboutput-1)])'''

def findDSdebug(rot, W0, WC): #OK!
    rotprim = []
    for i in range(nboutput):
        rotprim.append((rot[i] - ((powA[i] * W0 + polA[i] * WC) >> (2 * k - known_up))) % (1<<known_up))
    Y = [((S[i] - (powA[i] * W0 + polA[i] * WC)) >> known_low) % (1 << (2*k - known_low)) for i in range(nboutput)]
    #print("rotprim<< et Y+")
    rotprim[1] = (rotprim[1] - 1) % k
    #print("ATTENTION UN ZERO !!!")
    #print(dec2bin(Y[1]))
    #for i in range(nboutput - 1):
    #   print(dec2bin(rotprim[i] << (2 * k - known_up - known_low)))
    #   print(dec2bin(Y[i]))
    #   print('')
    cheat_DS = [(Y[i + 1] - Y[i]) % (1 << (2*k - known_low)) for i in range(nboutput - 1)]
    tmp = vector([(rotprim[i+1] - rotprim[i]) << (2 * k - known_up - known_low) for i in range(nboutput - 1)])

    #print("cheat_DS et tmp :")
    #for i in range(nboutput - 1):
    #   print(dec2bin(cheat_DS[i]))
    #   print(dec2bin(tmp[i]))
    #   print((cheat_DS[i] - tmp[i])>>(2 * k - known_up - known_low))

    G = f.IntegerMatrix.from_matrix(Greduite2)
    return f.CVP.closest_vector(G,tuple(tmp))

cpt = {}
cpt[False, False] = 0
cpt[False, True] = 0
cpt[True,  False] = 0
cpt[True,  True] = 0

cptrotfail = 0
n = 1000
total_output_size = 0
#recupDonnees()
for blabla in range(n):
    X, S, c = sortiesGenerateur()    

    # génère ce qu'on est pas censé connaître, et que le gros calcul doit retrouver
    W0 = S[0] % (1 << known_low)
    WC = c % (1 << known_low)
    rot = []
    for i in range(nboutput):
        rot.append(S[i] >> (2 * k - 6))

    # refait la partie qui trouve DS64.
    Y = [((S[i] - (powA[i] * W0 + polA[i] * WC)) >> known_low) % (1 << (2*k - known_low)) for i in range(nboutput)]
    cheat_DS = [(Y[i + 1] - Y[i]) % (1 << (2*k - known_low)) for i in range(nboutput - 1)]
    assert vector(cheat_DS) in matrix(Greduite2).image()

    uX = unrotateX(X,rot)    
    DS64, Y0 = FindDS64(uX, rot, W0, WC, Greduite1)#OK!
    
    # vérifie que DS64 est correct
    assert cheat_DS[0] % (1 << 64) == DS64[0]

    # calcule toutes les rotations futures à partir de DS64
    tabrot = FindRot(DS64[0], X, Y0, W0, WC) #a l'air OK!
    
    listrot = list(product(*tabrot))
    #for rot in listrot:
    #   for  i in range(len(S)):
    #       print((S[i] >> (2 * k - 6)) - rot[i])
    
    # pour chaque jeu de rotations trouvées, calcule DS complet
    

    test_a = False
    test_b = False

    output = set()

    # vérifie qu'on a bien trouvé un des bon DS complet.
    listDS = [findDS(rot, W0, WC) for rot in listrot]

    output.update(listDS)
    for DS in listDS:
        Sprim = [(S[i] - polA[i] * (c % 1<<known_low) - powA[i] * (S[0] % 2^known_low)) % 2^128 for i in range(nboutput)]
        if(DS[0] == ((Sprim[1] - Sprim[0]) >> known_low)):
            test_a = True
    #if blabla < 3 :
    #   print('###############  NOUVEAU SAMPLE #################')
    #   listDS = [findDSdebug(rot, Greduite2, S, DS64) for rot in listrot]
    #if not test :
        #print('###############  NOUVEAU CAS QUI MARCHE PAS #################')
    listDS = [findDSdebug(rot, W0, WC) for rot in listrot]

    output.update(listDS)
    for DS in listDS:
        Sprim = [(S[i] - polA[i] * (c % 1<<known_low) - powA[i] * (S[0] % 2^known_low)) % 2^128 for i in range(nboutput)]
        if(DS[0] == ((Sprim[1] - Sprim[0]) >> known_low)):
            test_b = True
    
    cpt[test_a, test_b] += 1
    total_output_size += len(output)
print(n)
print("cpt:")
print(cpt)
print("total output size = {}".format(total_output_size))

