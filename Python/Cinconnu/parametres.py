####Variables de départ####

k = 64
M = 2**(k * 2)
a = 2549297995355413924*2**64+4865540595714422341

#mult par defaut de pcg
#c=6364136223846793005*2**64+1442695040888963407#increment par defaut de pcg (connu)

known_up=6
known_low=11
nbiter = 5 #nb d'iterations utilisées pour la résolution. 
nboutput = 100 #nb de sorties consécutives connues
n=2**k #modulo utilisé pour les calculs LLL 

#Matrice génératrice réduite grâce à Sage à rentrer (de taille nbiter - 1) et son inverse

Greduite = [[-186304953996472, -216211368070119,  110964501361298,  131252974561432],
[-126056243766680,   99587582169277,   -5646098666150, -233919070109448],
[   7937589136904, -214303762177807, -268280113597118,  -98716819647784],
[  93078431381544,   -1707551230219,  149382085707466, -134620659538888]]


invG = [[-2.04279952328856e-15, -2.93791683689260e-15,  6.74861642263602e-16,  2.61840245666112e-15],
[-1.98886046520741e-15,  1.10621147658696e-15, -2.12731308263381e-15, -2.30132753131773e-15],
[ 1.44762674070265e-15, -1.54769860122306e-16, -1.55490854567009e-15,  2.82055195172104e-15],
[ 2.19171469167945e-16, -2.21708502251561e-15, -1.23181629826941e-15, -2.45886210697988e-15]]


#Polynôme et puissances de a (utilisés pour le calcul de polW, polC, ...) :
powA = [1]
polA = [0]
for i in range(nboutput - 1):
    polA.append(powA[i] + polA[i])
    powA.append(powA[i] * a)