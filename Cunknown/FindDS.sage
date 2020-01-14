import time
import random as r
k = 64
known_up  = 8
known_low  = 11
a = 2549297995355413924 * 2^64 + 4865540595714422341
nboutput = 30


def sortiesGenerateur():#OK !
    c = (r.randint(0, 2**(2 * k)) * 2 + 1) % 2**(2 * k) #c est impair
    S=[r.randint(0,2**(k * 2))]
    #c=6364136223846793005 * 2^64 + 1442695040888963407#increment par defaut de pcg (connu)
    #S=[8487854484825256858 + 11929896274893053136 * 2^64]
    for i in range (nboutput-1):
        S.append((S[i] * a + c) % (1<<(2 * k)))
    X=[]
    for i in range (nboutput):
        x=(S[i] % (1<<k))^^(S[i]>>k)
        rot = S[i] >> (2 * k - known_up)
        X.append((x >> rot) | ((x << (k - rot)) % (1 << k)))
    return X,S,c


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
    return G.transpose().LLL().transpose()
    
def getInvG(Greduite):
    return Greduite.inverse().n()

def findDS(rot, Greduite, invG):
    tmp = vector([(rot[i+1] - rot[i]) << (2 * k - known_up - known_low) for i in range(nboutput - 1)])
    u = invG * tmp
    tmp = vector([round(u_) for u_ in u])
    return Greduite * tmp

Greduite = getGreduite(nboutput - 1, 2^(2 * k - known_low))
invG = getInvG(Greduite)
cpt = 0
for blabla in range(1000):
    X, S,c = sortiesGenerateur()

    polA = [0];
    powA = [1];
    for i in range (1, nboutput):
        polA.append((polA[i-1] + powA [i-1]) % 2^128)
        powA.append((powA[i-1] * a)% 2^128)

    Sprim = [(S[i] - polA[i] * (c % 2^known_low) - powA[i] * (S[0] % 2^known_low)) % 2^128 for i in range(nboutput)]
    #print((Sprim[1] - Sprim[0]) >> known_low)
    rot = [Sprim[i] >> (k * 2 - known_up) for i in range(nboutput)]
    DS = findDS(rot, Greduite, invG)
    if(DS[0] == ((Sprim[1] - Sprim[0]) >> known_low)):
        cpt += 1
    #print(DS[0])
    #print(12615681514276467327 * 2^64 + 8299778918817149495)
print(nboutput)
print(cpt)