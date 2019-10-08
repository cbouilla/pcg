from RecupS import *
from math import log

cpt=0
for i in range(10000):
    X,S=sortiesGenerateur(k,a,c,nbiter) #S n'est pas utilisé sauf dans la version de vérification des résultats

    W0=S[0] % 2**connus_low #version pour tests de "qualité de prédiction"
    rot = [s//2**(k - connus_up) for s in S]

    S2 = FindS(X, rot, W0)
    cpt += (S == S2)
print(cpt/10000)