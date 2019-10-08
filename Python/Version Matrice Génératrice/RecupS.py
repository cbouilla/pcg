from parametres import *

###Fonctions utilisées###

def getPolW(W0): #OK 128bits
    polW=[W0]
    for i in range(nbiter - 1):
        polW.append((polW[i] * a) % 2**k)
    return polW

#retrouve Yprim avec vecteur des nb_iter rotations ( < 2**connus_up)
def getY(polW, rot, X): #OK ! 64bits
    Y=[(((polW[i] + polC[i]) % 2**connus_low) ^ (X[i] % 2**connus_low) ) * 2**connus_up + (rot[i] ^ (X[i] // 2**(k//2 - connus_up))) for i in range(nbiter)]
    return Y

def getYprim(Y, polC, polW): #OK ! avec erreurs de retenues ~64bits (polC polW)
    Yprim=[(Y[i] - (polC[i] + polW[i]) // 2**(k//2 - connus_up)) % 2**(connus_up + connus_low) for i in range(nbiter)]
    return Yprim


###Remontee au S correspondant
def FindSprim(X, rot, polW): #rajouter rot dans la version non test ? #OK! ~64bits
    #polW = getPolW(W0)
    Y = getY(polW, rot, X)
    Yprim = getYprim(Y, polC, polW) #OK avec erreurs de retenues!
    tmp = [y * 2**(k//2 - connus_up - connus_low) for y in Yprim]#on rajoute les zéros, recentrage impossible à cause des erreurs de retenues
    u = prodMatVec(invG, tmp)
    print("u")
    print(u)
    Sprim = prodMatVec(Greduite, [round(u_) for u_ in u])
    
    return Sprim

def FindS(X, rot, W0):#OK ! ~64bits
    polW=getPolW(W0)
    Sprim = FindSprim(X, rot,polW)
    Smod = [Sprim[i] * 2**(connus_low) + ((polW[i] + polC[i]) % 2**(k//2 + connus_low)) for i in range(nbiter)]
    S = [((Smod[i] % 2**(k//2)) ^ X[i]) * 2**(k//2) + (Smod[i] % 2**(k//2)) for i in range(nbiter)]
    return S