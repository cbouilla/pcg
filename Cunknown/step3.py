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


def sortiesGenerateur(S0, c, k):
    S = [S0]
    for i in range (k - 1):
        S.append((S[i] * a + c) % N)
    X = []
    for i in range(k):
        x = (S[i] % K) ^ (S[i] >> 64)
        rot = S[i] >> 122
        X.append((x >> rot) | ((x << (64 - rot)) % K))
    return X


if __name__ == '__main__':
    X = [0] * 48
    X[ 1] = 0xd4166f4c3e02d10a;
    X[ 2] = 0x1d1ceb21e7737101;
    X[ 3] = 0xf8b90f473a5426d3;
    X[ 4] = 0xe3a3b7babb2ad9ca;
    X[ 5] = 0x0077f2c80987dd13;
    X[ 6] = 0xf8ddaf2431548a13;
    X[ 7] = 0x80935e041bbab85a;
    X[ 8] = 0xbe0fde3939201c50;
    X[ 9] = 0xe9604fdf6b2177b7;
    X[10] = 0x95d9cf24a229cedf;
    X[11] = 0x0434a85418759293;
    X[12] = 0x04d230c1debe7999;
    X[13] = 0x83a3cd257d4d04b0;
    X[14] = 0x13990a23037c13c4;
    X[15] = 0xfbfaea2e50411202;
    X[16] = 0x421a394a36baebf8;
    X[17] = 0x5a878c4594ea7221;
    X[18] = 0xbd37307bdd522b9c;
    X[19] = 0x39af06b3e9b3ae10;
    X[20] = 0x1d21f1cee77b8e2e;
    X[21] = 0xe4bddad0aacaf420;
    X[22] = 0x1009ed344cd7f2f4;
    X[23] = 0xff287c5797cceb71;
    X[24] = 0x85e8968b0d8ab49b;
    X[25] = 0x69a9b821830862cc;
    X[26] = 0xf9fe65ed23740aea;
    X[27] = 0x47669184bc43d948;
    X[28] = 0xe3f19b31915ae6d3;
    X[29] = 0x5f3945718b6dac44;
    X[30] = 0x49bfacfe8056b33c;
    X[31] = 0xf2b358ceb722628f;
    X[32] = 0xb4d1bf17f2c57b71;
    X[33] = 0xb4300000ad802deb;
    X[34] = 0xe3125e8022de888d;
    X[35] = 0x2c7f6d404196ed4d;
    X[36] = 0xb10490274ecbe897;
    X[37] = 0xb04da8a406da3814;
    X[38] = 0xbc70124be9a196c9;
    X[39] = 0xf81a244f765141e6;
    X[40] = 0x20413cda8442a149;
    X[41] = 0xf654f8084029c557;
    X[42] = 0xa1677f2f18f8484c;
    X[43] = 0x3ef4f6e355c53e70;
    X[44] = 0x30cceccd8f73c567;
    X[45] = 0xf22a193ded925cb6;
    X[46] = 0xa9c0cba2f6abbd89;
    X[47] = 0xf26c7311f52ededb;
    del X[0]
    original_X = X   

    W_0, W_c, start_DS, rot = 0x01d0, 0x035b, 0x65cca9e25d817ee1460ae35556d5069b, (28, 53, 44, 7, 47, 22, 2, 10, 60, 27, 46, 27, 28, 2, 44, 10, 36, 10, 62, 42, 20, 40, 15, 14, 7, 59, 44, 46, 59, 9, 46, 10, 61, 5, 28, 53, 19, 32, 18, 55)
    # W_0, W_c, start_DS, rot = 0x01d0, 0x035b, 0x65cca9e25d817ee1460ae35556d5069b, (28, 53, 44, 7, 47, 22, 2, 10, 60, 27, 46, 27, 28, 2, 44, 10, 36, 10, 62, 33, 20, 40, 15, 14, 7, 59, 44, 46, 59, 9, 46, 10, 61, 5, 28, 53, 19, 32, 18, 55)
    # W_0, W_c, start_DS, rot = 0x05d0, 0x035b, 0x65cca9e25d817ee1460ae35556d5069b, (28, 53, 44, 7, 47, 22, 2, 10, 60, 27, 46, 27, 28, 2, 44, 10, 36, 10, 62, 42, 20, 40, 15, 14, 7, 59, 44, 46, 59, 9, 46, 10, 61, 5, 28, 53, 19, 32, 18, 55)
    # W_0, W_c, start_DS, rot = 0x05d0, 0x035b, 0x65cca9e25d817ee1460ae35556d5069b, (28, 53, 44, 7, 47, 22, 2, 10, 60, 27, 46, 27, 28, 2, 44, 10, 36, 10, 62, 33, 20, 40, 15, 14, 7, 59, 44, 46, 59, 9, 46, 10, 61, 5, 28, 53, 19, 32, 18, 55)

    uX = unrotateX(X, rot, nboutput)
    DS = [(polA[i] * start_DS) % N for i in range(nboutput)]

    DS = [to_bits(x, 128) for x in DS]
    X = [to_bits(x, 64) for x in uX]
    rot = [to_bits(x, 6) for x in rot]

    known_S = [-1] * 128
    known_Carry = []
    for _ in range(nboutput):
        known_Carry.append([-1] * 128)

    for j in range(6):
        known_S[122+j] = rot[0][j]
        known_S[58+j] = known_S[122+j] ^ X[0][58+j]

    for i in range(nboutput):
        known_Carry[i][0] = 0
        known_Carry[i][64] = (X[i][0] ^ X[0][0]) ^ (DS[i][0] ^ DS[i][64])

    for j in range(58):
        # on connaît S[0][0:j], Carry[*][0:j+1], S[0][64:64+j], Carry[*][64:64+j+1]
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

        known_S[j] = (X[i][j+1] ^ X[0][j+1]) ^ DS[i][65+j] ^ DS[i][j+1] ^ DS[i][64+j]
        known_S[j+64] = X[0][j] ^ known_S[j]
        
        for i in range(nboutput):
            known_Carry[i][j+1] = MAJ(known_S[j], DS[i][j], known_Carry[i][j])
            known_Carry[i][j+65] = MAJ(known_S[64+j], DS[i][64+j], known_Carry[i][64+j])

    found_S0 = 0
    for i in range(128):
        found_S0 += known_S[i] << i
    found_S1 = (found_S0 + start_DS) % N
    found_c = (found_S1 - a*found_S0) % N

    XX = sortiesGenerateur(found_S0, found_c, 64)
    for i, x in enumerate(XX):
        print("X[{:2d}] = 0x{:016x};".format(i, x))
