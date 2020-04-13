import time
import fpylll
import os
from itertools import product
from math import log
import random

k = 64
known_up  = 6
known_low  = 13
a = 2549297995355413924 * 2**64 + 4865540595714422341
nbiter = 5
nboutput = 40
N = 2**128
K = 2**64
LOW = 2**known_low

polA = [0];
powA = [1];
for i in range (1, nboutput):
    polA.append((polA[i-1] + powA [i-1]) % N)
    powA.append((powA[i-1] * a) % N)

def sortiesGenerateur(k): #OK !
    c = random.randrange(N) | 1
    S = [random.randrange(N)]
    for i in range (k - 1):
        S.append((S[i] * a + c) % N)
    X = []
    for i in range(k):
        x = (S[i] % K) ^ (S[i] >> 64)
        rot = S[i] >> 122
        X.append((x >> rot) | ((x << (64 - rot)) % K))
    return X, S, c


## Unrotate
def rotateX(X, rot, k):#OK !
    rX = [];
    for i in range(k):
        rX.append((X[i] >> rot[i]) + ((X[i] << (64 - rot[i])) % K))
    return rX

def unrotateX(X, rot, k):#OK !
    rot2 = []
    for i in range(k):
        rot2.append(64 - rot[i])
    return rotateX(X, rot2, k)

def to_bits(x, k):
    r = []
    for _ in range(k):
        r.append(x & 1)
        x >>= 1
    return r

def carries(x, y):
    """Renvoie les retenues générées lors du calcul de x + y.
    C'est ce qu'il faut XORer sur x^y pour obtenir x+y
    """
    return (x+y) ^ (x ^ y)

def MAJ(a,b,c):
    return (a and b) or (a and c) or (b and c)

if __name__ == '__main__':
    # teste la proba de succès des algorithmes
    n = 1000
    success = 0

    X, S, c = sortiesGenerateur(nboutput)    
    
    # génère ce qu'on est pas censé connaître, et que le gros calcul doit retrouver
    W0 = S[0] % LOW
    WC = c % LOW
    rot = []
    for i in range(nboutput):
        rot.append(S[i] >> 122)
    uX = unrotateX(X, rot, nboutput)
    start_DS = (S[1] - S[0]) % N
    DS = [(polA[i] * start_DS) % N for i in range(nboutput)]
    Carry = [carries(S[0], DS[i]) % N for i in range(nboutput)]

    for i in range(1, nboutput):
        assert S[i] == (S[0] + DS[i]) % N
        assert S[i] == S[0] ^ DS[i] ^ Carry[i]
    
    # Goal : récupérer S[1] et S[0]
    S0 = S[0]
    S1 = S[1]

    S = [to_bits(x, 128) for x in S]
    DS = [to_bits(x, 128) for x in DS]
    X = [to_bits(x, 64) for x in uX]
    Carry = [to_bits(x, 128) for x in Carry]
    rot = [to_bits(x, 6) for x in rot]

    for i in range(nboutput):
        # the rotations give 6 bits
        for j in range(6):
            assert S[i][122+j] == rot[i][j]
            assert S[i][58+j] == S[i][122+j] ^ X[i][58+j]

        for j in range(128):
            assert S[i][j] == S[0][j] ^ DS[i][j] ^ Carry[i][j]
        
        # pas de retenue sur le bit de poids faible
        assert Carry[i][0] == 0
        
        for j in range(64):
            # ceci permet de passer de carry[i][j] à Carry[i][j+64] et vice-versa
            assert (X[i][j] ^ X[0][j]) ^ (DS[i][j] ^ DS[i][j+64]) == Carry[i][j] ^ Carry[i][j+64]

        # en particulier
        assert (X[i][0] ^ X[0][0]) ^ (DS[i][0] ^ DS[i][64]) == Carry[i][64]

        for j in range(1, 128):
            # ceci permet de passer de Carry[i][j-1] et S[0][j-1] à Carry[i][j]
            assert Carry[i][j] == MAJ(S[0][j-1], DS[i][j-1], Carry[i][j-1])

        assert S[i][1] == S[0][1] ^ DS[i][1] ^ Carry[i][1]
        assert S[i][1] == S[0][1] ^ DS[i][1] ^ MAJ(S[0][0], DS[i][0], Carry[i][0])
        assert S[i][1] == S[0][1] ^ DS[i][1] ^ (S[0][0] and DS[i][0])
        
        for j in range(1, 64):
            assert S[i][64+j] == S[0][64+j] ^ DS[i][64+j] ^ MAJ(S[0][63+j], DS[i][63+j], Carry[i][63+j])
            
            # l'équation VRAIMENT utile est la suivante :
            assert (X[i][j] ^ X[0][j]) ^ DS[i][64+j] ^ DS[i][j] == MAJ(S[0][j-1], DS[i][j-1], Carry[i][j-1]) ^ MAJ(S[0][63+j], DS[i][63+j], Carry[i][63+j])

        # bonnes conditions: DS[i][0] == 1,  DS[i][64] == Carry[i][64].

######################################################################################################""

    known_S = [-1] * 128
    known_Carry = []
    for _ in range(nboutput):
        known_Carry.append([-1] * 128)

    for j in range(6):
        known_S[122+j] = rot[0][j]
        known_S[58+j] = known_S[122+j] ^ X[0][58+j]

    for k in range(122, 128):
        assert known_S[k] == S[0][k]
    for k in range(58, 64):
        assert known_S[k] == S[0][k]

    for i in range(nboutput):
        known_Carry[i][0] = 0
        known_Carry[i][64] = (X[i][0] ^ X[0][0]) ^ (DS[i][0] ^ DS[i][64])

    for j in range(58):
        print("j = {}".format(j))
        # on connaît S[0][0:j], Carry[*][0:j+1], S[0][64:64+j], Carry[*][64:64+j+1]

        # vérifie que tous les précédents sont OK
        for k in range(j):
            assert known_S[k] != -1
            assert known_S[k] == S[0][k]
        for k in range(j + 1):
            for i in range(nboutput):
                assert known_Carry[i][k] != -1
                assert known_Carry[i][k] == Carry[i][k]

        # trouve i>0 tq DS[i][j-1] != Carry[i][j-1],  DS[i][64] == Carry[i][64]
        # ça marcherait aussi en échangeant == et !=
        i = None
        for k in range(1, nboutput):
            if DS[k][j] != known_Carry[k][j] and DS[k][64+j] == known_Carry[k][64+j]:
                i = k
                break
        if i is None:
            raise ValueError("no i matching conditions")
        print("OK with i = {}".format(i))

#        assert MAJ(S[0][j], DS[i][j], Carry[i][j]) == S[0][j]
#        assert MAJ(S[0][64+j], DS[i][64+j], Carry[i][64+j]) == DS[i][64+j]
#        assert (X[i][j+1] ^ X[0][j+1]) ^ DS[i][65+j] ^ DS[i][j+1] ^ DS[i][64+j] == S[0][j]
        
        known_S[j] = (X[i][j+1] ^ X[0][j+1]) ^ DS[i][65+j] ^ DS[i][j+1] ^ DS[i][64+j]
        known_S[j+64] = X[0][j] ^ known_S[j]
        
        for i in range(nboutput):
            known_Carry[i][j+1] = MAJ(known_S[j], DS[i][j], known_Carry[i][j])
            known_Carry[i][j+65] = MAJ(known_S[64+j], DS[i][64+j], known_Carry[i][64+j])


    found_S0 = 0
    for i in range(128):
        found_S0 += known_S[i] << i

    assert found_S0 == S0
    print("{:032x}".format(found_S0))

    found_S1 = (found_S0 + start_DS) % N
    assert found_S1 == S1

    assert (found_S1 - a*found_S0) % N == c