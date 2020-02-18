import time
import random as r
import fpylll as f
import os

k = 64
known_up  = 6
known_low  = 13
a = 2549297995355413924 * 2^64 + 4865540595714422341
nbiter = 5
nboutput = 40


polA = [0];
powA = [1];
for i in range (1, nboutput):
    polA.append((polA[i-1] + powA [i-1]) % 2^128)
    powA.append((powA[i-1] * a)% 2^128)


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

def getGreduite(n,mod):
    G = getG(n,mod)
    Gr = G.transpose().LLL()
    Glll = f.IntegerMatrix.from_matrix(Gr)
    return Glll

Greduite1 = getGreduite(nbiter - 1, 2^k)

Greduite2 = getGreduite(nboutput - 1, 2^(2 * k - known_low))
print("Vecteur réputé court de taille : 2^{:.1f}".format(float(log(Greduite2[0].norm(), 2).n())))

#### Récupération des données ####
'''def recupDonnees():
    listfiles = os.listdir('results')
    print(listfiles)
    listDonnees = []
    for i in range(listfiles):
        listDonnees.append[]
        f = open(listfiles[i], "r")
        while(1):
            s = f.readline()'''
            


def sortiesGenerateur():#OK !
    c = (r.randint(0, 1<<(2 * k)) * 2 + 1) % 1<<(2 * k) #c est impair
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
def FindDS64(uX, rot, W0,WC, Greduite): #rajouter rot dans la version non test ? #OK! ~64bits
    Y = getY(W0, WC, rot, uX)
    DY = getDY(Y, WC, W0) #OK avec erreurs de retenues!
    tmp = vector([y * 1<<(k - known_up - known_low) for y in DY])#on rajoute les zéros, recentrage impossible à cause des erreurs de retenues
    DS64 = f.CVP.closest_vector(Greduite,tuple(tmp))
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

def findDS(rot, Greduite, cheat_DS): #OK!
    rotprim = []
    for i in range(nboutput):
        rotprim.append((((rot[i] << 122) - (powA[i] * W0 + polA[i] * WC)) >> known_low) % (1 << (128 - known_low)))
    
    # approximation de cheat_DS
    tmp = vector([(rotprim[i+1] - rotprim[i])  for i in range(nboutput - 1)])
    
    # distance entre approx et vraie solution
    norm = sqrt(sum([(tmp[i] - cheat_DS[i])**2 for i in range(nboutput - 1)]))

    print("distance : 2^{:.1f}".format(float(log(norm, 2).n())))

    return f.CVP.closest_vector(Greduite,tuple(tmp))


def recFindDS(rot, tabrot, Greduite, i, cheat_DS):
    DS = []
    if(i == nboutput):
        DS = [findDS(rot, Greduite, cheat_DS)]
        return(DS)
    for r in tabrot[i]:
        rot.append(r)
        DS += recFindDS(copy(rot), tabrot, Greduite, i+1, cheat_DS)
    return(DS)


cpt = 0
cptrotfail = 0
n = 1000
#recupDonnees()
for blabla in range(n):
    X, S, c = sortiesGenerateur()    

    W0 = S[0] % (1 << known_low)
    WC = c % (1 << known_low)
    rot = []
    for i in range(nboutput):
        rot.append(S[i] >> (2 * k - 6))

    Y = [((S[i] - (powA[i] * W0 + polA[i] * WC)) >> known_low) % (1 << (2*k - known_low)) for i in range(nboutput)]
    cheat_DS = [(Y[i + 1] - Y[i]) % (1 << (2*k - known_low)) for i in range(nboutput - 1)]
    assert vector(cheat_DS) in matrix(Greduite2).image()

    uX = unrotateX(X,rot)    
    DS64, Y0 = FindDS64(uX, rot, W0,WC, Greduite1)#OK!
    tabrot = FindRot(DS64[0],X, Y0, W0, WC)#a l'air OK!
    test = 0
    if(len(tabrot) == 0):
        cptrotfail += 1
    else:
        rot = []
        listDS = recFindDS(rot, tabrot, Greduite2, 0, cheat_DS)
        for DS in listDS:
            Sprim = [(S[i] - polA[i] * (c % 1<<known_low) - powA[i] * (S[0] % 2^known_low)) % 2^128 for i in range(nboutput)]
            if(DS[0] == ((Sprim[1] - Sprim[0]) >> known_low)):
                test = 1
    cpt += test
print(n)
print("cpt:")
print(cpt)
print("cptrotfail:")
print(cptrotfail)
