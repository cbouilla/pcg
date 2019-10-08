from RecupS import *
from math import log

X,S=sortiesGenerateur(k,a,c,nbiter) #S n'est pas utilisé sauf dans la version de vérification des résultats
nbtour = 2**15

t1 = time.time()
for i in range(nbtour) :
    W0 = r.randint(0, 2**connus_low)
    rot = [r.randint(0, 2**connus_up) for i in range(nbiter)]
    S2 = FindS(X, rot, W0)
    if(S2 == S):
        print("J'AAAIII !!!")
t2=time.time()
print(t2 - t1)