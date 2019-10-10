from RecupS import *
from math import log

X,S=sortiesGenerateur(k,a,c,nbiter) #S n'est pas utilisé sauf dans la version de vérification des résultats
print("X")
print(X)
rot = [s // 2**(k-connus_up) for s in S]
print("rot")
print(rot)
print(polC)
W0=S[0] % 2**connus_low #version pour tests de "qualité de prédiction"
print("W0")
print(W0)
rot = [s//2**(k - connus_up) for s in S]

polW = getPolW(W0)
print("polW")
print(polW)
Y = getY(polW, rot, X)
print("Y")
print(Y)
Yprim = getYprim(Y, polC, polW)
print("Yprim")
print(Yprim)

Sprim = FindSprim(X,rot,polW)
print("Sprim")
print(Sprim)
S2 = FindS(X, rot, W0)
print("S")
print(S)

