from fonctions import *

####Variables de départ####
k=128
m=2**k
a=2549297995355413924*2**64+4865540595714422341
#mult par defaut de pcg
c=6364136223846793005*2**64+1442695040888963407#increment par defaut de pcg (connu)
connus_up=6
connus_low=18
nbiter=3 #nb d'iterations calculées. On considère qu'il s'agit aussi du nombre de sorties consécutives connues
n=2**((k//2)+connus_low) #modulo utilisé pour les calculs LLL 

#Calcul du polC (commun à tous les itérations et mêmes tous les appels à ce générateur utilisant c comme incrément)
polC=[0]
for j in range(nbiter-1):
    polC.append((polC[j]+(c*(a**j))%2**k)%2**k)

#Matrice génératrice réduite grâce à Sage à rentrer 

Greduite = [[-1241281756092, -5001120657083,  8655886039732],
[ 3827459685972, -2117155768935,  3303731088004],
[ -728312298332,  5479732607037,  6319848582548]]

invG = [[-9.25221813226351e-14,  2.32272588749499e-13,  5.30001389997814e-15],
[-7.81560146462246e-14, -4.52719459143047e-15,  1.09411828555506e-13],
[ 5.71040610630735e-14,  3.06929503514353e-14,  6.39750297008745e-14]]

'''Plus nécéssaire
#grand premiers justes inférieurs à 2**32
#p1=18446744073709551653#4294967291
#p2=18446744073709551629#4294967279

#Calcul des factorisations LU
L1,U1=FactorisationLUMod(Lreduite,p1)
L2,U2=FactorisationLUMod(Lreduite,p2)'''

''' Passage de S à Y , Z et W (idem pour utilisation de polC et polW)
W = S % 2**connus_low
Y = (S % 2**(k//2 + connus_low)) // 2**(k//2 -connus_up)
Z = (S % 2**(k//2 -connus_up)) // 2**connus_low

Smod = Y * 2**(k//2 -connus_up) + Z * 2**connus_low + W '''