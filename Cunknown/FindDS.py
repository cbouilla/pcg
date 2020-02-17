import itertools
import math
import time
import random
import fpylll
import os

k = 64
M = 1 << k
N = 1 << (2*k)
known_up  = 6
known_low  = 13
L = 1 << known_low
U = 1 << known_up
a = (2549297995355413924 << 64) + 4865540595714422341
nbiter = 5
nboutput = 14


powA = [1]; # powA[i] == a**i % N
polA = [0]; # polA[i] == a**0 + a**1 + ... + a**(i-1) % N
for i in range(1, nboutput):
    polA.append((polA[i - 1] + powA[i - 1]) % N)
    powA.append((powA[i - 1] * a) % N)


####   MATRIX   ####
def getG(n, mod):
    """
    Renvoie la matrice [n] x [n] dont la première colonne est formée
    des puissances successives de a, avec des mod sur le reste
    de la diagonale.

    Un vecteur de la forme [x, x*a, x*a^2, x*a^3, ....] modulo [mod] 
    s'exprime comme une combinaison linéaire à coefficients sur ZZ des 
    colonnes de cette matrice

    La matrice est renvoyée en tant qu'objet fplll.

    >>> print(getG(4, 2**32))
    [          1          0          0          0 ]
    [ 2681009733 4294967296          0          0 ]
    [  551595673          0 4294967296          0 ]
    [ 3980531005          0          0 4294967296 ]
    """
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
    return fpylll.IntegerMatrix.from_matrix(G)


def getGreduite(n, mod):
    """
    Même chose, mais la matrice est LLL-réduite et transposée.

    Un vecteur de la forme [x, x*a, x*a^2, x*a^3, ....] modulo [mod] 
    s'exprime comme une combinaison linéaire à coefficients sur ZZ des 
    LIGNES de cette matrice

    La matrice est renvoyée en tant qu'objet fplll.
    """
    G = getG(n, mod)
    return fpylll.LLL.reduction(G.transpose())

# these names are REALLY bad
Greduite1 = getGreduite(nbiter - 1, M)
svp = fpylll.SVP.shortest_vector(Greduite1)
norm = sum([x * x for x in svp])
print("(petite matrice) Norme SVP : 2^{:.3f}".format(0.5 * math.log(norm, 2)))

print("Ceci est normalement plus grand que 2^{}".format(58 - known_low))

Greduite2 = getGreduite(nboutput - 1, N // L // M)
svp = fpylll.SVP.shortest_vector(Greduite2, method="fast")
svp2_norm = math.sqrt(sum([x * x for x in svp]))
print("(grosse matrice) Norme SVP : 2^{:.3f}".format(math.log(svp2_norm, 2)))

def sortiesGenerateur():
    """
    Produit des sorties de PCG64, avec graine aléatoire.
    Renvoie la séquence des états internes successifs, les sorties, et l'incrément choisi.
    """
    c = random.randrange(N) | 1       #c est impair
    S = [random.randrange(N)]
    #c=6364136223846793005 * 2^64 + 1442695040888963407#increment par defaut de pcg (connu)
    for i in range(nboutput - 1):
        S.append((S[i] * a + c) % N)
    X = []
    for i in range(nboutput):
        x = (S[i] % M) ^ (S[i] >> k)
        rot = S[i] >> (2 * k - known_up)
        out = (x >> rot) | ((x << (k - rot)) % M)
        X.append(out)
    return X, S, c


## Unrotate
def rotateX(X, rot):
    rX = [];
    for i in range(nbiter):
        rX.append(((X[i] // 2**rot[i]) + ((X[i] * 2**(k - rot[i])) % M)))
    return rX

def unrotateX(X, rot):#OK !
    rot2 = []
    for i in range(nbiter):
        rot2.append((k - rot[i]) % k)
    return rotateX(X, rot2)

def unrotate1(Xi):
    """
    Renvoie Xi <<< 1 (i.e. rotation vers la gauche de 1), modulo 2^k.

    >>> unrotate1(0x0001)
    2
    >>> unrotate1(0x8000000000000000)
    1
    """
    return (Xi >> (k - 1)) | ((Xi << 1) % M)


###### Sous-fonctions de FindDS######
def getY(W0, WC, rot, uX):#OK !
    """
    Y[i] = bits [58:64 + known_low] de X[i], 
    """
    Y = []
    for i in range(nbiter):
        bits_0_known_low = (powA[i] * W0 + polA[i] * WC) % L
        bits_64_known_low = (bits_0_known_low ^ uX[i]) % L
        bits_58_64 = rot[i] ^ (uX[i] >> 58)
        y = (bits_58_64 ^ (bits_64_known_low << 6))
        Y.append(y)
    return Y

def getYprim(Y, WC, W0): #OK ! avec erreurs de retenues ~64bits (polC polW)
    """
    X'[i] == X[i] où les [known_low] bits de poids faibles de X[0] et c auraient été zéro.
             a priori, X'[i] est un multiple de L.
    Y'[i] = bits [58:64 + known_low] de X'[i].
    """
    Yprim = []
    for i in range(nbiter):
        yp = (Y[i] - (polA[i] * WC + powA[i]*W0) // 2**58) % 2**(6 + known_low)
        Yprim.append(yp)
    return Yprim

def getDY(Y, WC, W0): #OK ! avec erreurs de retenues ~64bits (polC polW)
    """
    DY = différence entre deux Y' successifs
    """
    Yprim = getYprim(Y, WC, W0)
    DY = []
    for i in range(nbiter-1):
        d = (Yprim[i + 1] - Yprim[i]) % 2**(6 + known_low)
        DY.append(d)
    return DY

######FINDDS######
def FindDS64(uX, rot, W0,WC, Greduite): #rajouter rot dans la version non test ? #OK! ~64bits
    """
    On a un vecteur de différences approximatives, dont les versions exactes
    appartiennent à l'espace vectoriel engendré par [1, a, a^2, a^3, .....] modulo M.

    Récupère le vecteur des différences modulo 2^64, donc par extension modulo 2^(64+known_low).
    """
    Y = getY(W0, WC, rot, uX)
    DY = getDY(Y, WC, W0) #OK avec erreurs de retenues!
    tmp = [y << (58 - known_low) for y in DY] #on rajoute les zéros, recentrage impossible à cause des erreurs de retenues
    DS64 = fpylll.CVP.closest_vector(Greduite, tmp, method="fast")
    return DS64, Y[0]

######FINDROTI######
#DS64ij = ((polA[j] - polA[i])*DSmod0) % 2**k

def FindRoti(DS640, X, i, Y0, W0, WC):#OK !
    """
    Renvoie la liste des valeurs possibles de la i-ème rotation.
    """
    DS640i = (polA[i] * DS640) % (1<<k)
    DSmod0i = ((DS640i << known_low) + W0 * powA[i] + WC * polA[i] - WC - W0) % (1<<(k +known_low))
    # Yi = vraiYi ou vraiYi - 1 à cause de la retenue
    Yi1 = (Y0 + (DSmod0i >> (k - known_up))) % (1 << (known_low + known_up))#avec ou sans retenue
    Yi2 = Yi1 + 1
    Wi = (W0 * powA[i] + WC * polA[i]) % (1 << known_low)
    roti = []
    for j in range(k):
        test1 = (((X ^ (Yi1 >> known_up)) % (1 << known_low)) == Wi) and ((j ^ (X >> (k - known_up))) == Yi1 % (1 << known_up))
        test2 = (((X ^ (Yi2 >> known_up)) % (1 << known_low)) == Wi) and ((j ^ (X >> (k - known_up))) == Yi2 % (1 << known_up))
        if (test1 or test2) :
            roti.append(j)
        X = unrotate1(X)
    return roti

def FindRot(DS640, X, Y0, W0, WC):
    """
    Détermine les rotations des [nboutput] premiers tours.
    """
    tabrot = []
    for i in range(nboutput):
        tabrot.append(FindRoti(DS640, X[i], i, Y0, W0, WC))
        if len(tabrot[i]) == 0:
            raise ValueError("cannot determine rotations!")
    return tabrot

def findDS(rot, Greduite, W0, WC, DS64, cheat_DS): #OK!
    """
    A partir de toutes les rotations, récupère l'ensemble des différences modulo 2^128
    """
    # check inputs
    for i in range(nboutput - 2):
        assert cheat_DS[i+1] == (a * cheat_DS[i]) % (N // L)
    for i in range(nboutput - 1):
        assert cheat_DS[i] == (pow(a, i, N // L) * cheat_DS[0]) % (N // L)
    dl = DS64[0]
    dh = cheat_DS[0] // M
    for i in range(nboutput - 1):
        assert cheat_DS[i] == (pow(a, i, N // L) * (dl + M * dh)) % (N // L)
        assert (cheat_DS[i] - pow(a, i, N // L) * dl) % (N // L) == (M * pow(a, i, N // L) * dh) % (N // L)
        assert (cheat_DS[i] - pow(a, i, N // L) * dl) % M == 0
        assert ((cheat_DS[i] - pow(a, i, N // L) * dl) // M) % (N // L // M) == (pow(a, i, N // L) * dh) % (N // L // M)

    shifted_DS = [((cheat_DS[i] - pow(a, i, N // L) * dl) // M) % (N // L // M) for i in range(nboutput - 1)]
    for i in range(nboutput - 1):
        assert shifted_DS[i] == (pow(a, i, N // L) * dh) % (N // L // M)

    rotprim = []
    for i in range(nboutput):
        rp = ((((rot[i]) << 122) - (powA[i] * W0 + polA[i] * WC)) // L) % (N // L)
        rotprim.append(rp)
    
    # approx est voisin de cheat_DS
    approx = []
    for i in range(nboutput - 1):
        approx.append(((rotprim[i + 1] - rotprim[i])) % (N // L))

    # shifted est voisin de shifted_DS
    shifted = [((approx[i] - pow(a, i, N // L) * dl) // M) % (N // L // M) for i in range(nboutput - 1)]
 
    # vérifie que la solution apartient bien au réseau ---> c'est toujours le cas
    test = fpylll.fplll.svpcvp.closest_vector(Greduite, shifted_DS, method='proved')
    assert test == tuple(shifted_DS)

    stuff = []
    norm = 0
    for (x, y) in zip(shifted, shifted_DS):
        norm += (x - y)**2
        stuff.append((x - y))
    
    #print("{:32x} --> {:32x}".format(N // L, N // L // 64))
    #for i in range(nboutput-2):
    #    print("{:32x} --> {:32x}".format(approx[i], (approx[i] - cheat_DS[i]) % (N // L)))

    if math.sqrt(norm) < svp2_norm:
        pass
        #print("Approx close enough (distance to solution = 2^{:.1f})".format(0.5 * math.log(norm, 2)))
    else:
        print("ERROR : approx is too far (distance to solution = 2^{:.1f})".format(0.5 * math.log(norm, 2)))
        raise ValueError("It will fail")

    cvp = fpylll.CVP.closest_vector(Greduite, shifted, method='proved')
    
    # check validity
    
    if cvp == tuple(shifted_DS):
        #print("CVP == solution")
        #for i in range(nboutput - 1):
        #    assert (cvp[i] % M) == DS64[i]
        return [cvp[i] * M + DS64[i] for i in range(nboutput - 1)]
    else:
        norm = 0
        for (x, y) in zip(shifted_DS, cvp):
            norm += (x - y)**2
        print("ERROR : distance CVP-solution = 2^{:.1f}".format(0.5 * math.log(norm, 2)))
        raise ValueError("Getting DS failed")



cpt = 0
multirot = 0
n = 1000
for blabla in range(n):
    X, S, c = sortiesGenerateur()
    W0 = S[0] % (1 << known_low)
    WC = c % (1 << known_low)
    rot = []
    for i in range(nboutput):
        rot.append(S[i] >> (2 * k - 6))

    uX = unrotateX(X,rot)
    DS64, Y0 = FindDS64(uX, rot, W0,WC, Greduite1)#OK!
    DS64 = [(powA[i] * DS64[0]) % M for i in range(nboutput - 1)]

    # cheat : vérifie DS64
    cheat_Y = [((S[i] - (polA[i] * WC + powA[i]*W0)) // L) % (N // L) for i in range(nboutput)]
    cheat_DS = [(cheat_Y[i+1] - cheat_Y[i]) % (N // L) for i in range(nboutput-1)]
    for i in range(nbiter - 1):
        assert DS64[i] == cheat_DS[i] % M
    for i in range(nbiter - 2):
        assert DS64[i + 1] == (a * DS64[i]) % M

    # calcule "toutes" les rotations
    # ce nom est MAL choisi !
    tabrot = FindRot(DS64[0], X, Y0, W0, WC)
    rotations = list(itertools.product(*tabrot))
    assert rotations != []
    if len(rotations) > 1:
        multirot += 1

    for rot in rotations:
        try:
            DS = findDS(rot, Greduite2, W0, WC, DS64, cheat_DS)
            Sprim = [(S[i] - (polA[i] * WC + powA[i] * W0)) % N for i in range(nboutput)]
            assert (Sprim[1] - Sprim[0]) % L == 0
            if DS[0] == ((Sprim[1] - Sprim[0]) // L) % (N // L):
                cpt += 1
            else:
                assert DS[0] == cheat_DS[0]
                raise TypeError("FindDS OK but test failed")
        except ValueError:
            pass
        
            #print(DS[0], cheat_DS[0])
    
print("n: {}".format(n))
print("cpt: {}".format(cpt))
print("multirot: {}".format(multirot))
