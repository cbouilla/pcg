import time
import fpylll
import os
from itertools import product
from math import log
import random

k = 64
known_up  = 6
known_low  = 11
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

def dec2bin(d,nb= 2 * k - known_low + 1):
    """Représentation d'un nombre entier en chaine binaire (nb: nombre de bits du mot)"""
    if d == 0:
        return "0".zfill(nb)
    if d<0:
        d += 1<<nb
    b=""
    while d != 0:
        d, r = divmod(d, 2)
        b = "01"[r] + b
    return b.zfill(nb)

####   MATRIX   ####
def getG(n, mod):
    G = fpylll.IntegerMatrix(n, n)
    G[0, 0] = 1
    powA = a % mod
    for i in range(1, n):
        G[i, 0] = powA
        powA = powA * a % mod
        for j in range(1, n):
            G[i, i] = mod
    return G

def getGreduite(n, mod):
    G = getG(n, mod)
    G.transpose()                     # in place
    fpylll.LLL.reduction(G)           # in place
    assert fpylll.LLL.is_reduced(G)
    return G

Greduite1 = getGreduite(nbiter - 1, 2**64)
Greduite2 = getGreduite(nboutput - 1, 2**(128 - known_low))


#### Récupération des données ####
#def recupDonnees():
#   listfiles = os.listdir('results')
#   print(listfiles)
#   listDonnees = []
#   for i in range(listfiles):
#       listDonnees.append[]
#       f = open(listfiles[i], "r")
#       while(1):
#           s = f.readline()
            


def sortiesGenerateur():#OK !
    c = random.randrange(N) | 1
    S = [random.randrange(N)]
    #c=6364136223846793005 * 2^64 + 1442695040888963407#increment par defaut de pcg (connu)
    #S=[8487854484825256858 + 11929896274893053136 * 2^64]
    for i in range (nboutput - 1):
        S.append((S[i] * a + c) % N)
    X = []
    for i in range(nboutput):
        x = (S[i] % K) ^ (S[i] >> 64)
        rot = S[i] >> 122
        X.append((x >> rot) | ((x << (64 - rot)) % K))
    return X, S, c


## Unrotate
def rotateX(X, rot):#OK !
    rX = [];
    for i in range(nbiter):
        rX.append((X[i] >> rot[i]) + ((X[i] << (64 - rot[i])) % K))
    return rX

def unrotateX(X, rot):#OK !
    rot2 = []
    for i in range(nbiter):
        rot2.append(64 - rot[i])
    return rotateX(X, rot2)

def unrotate1(Xi):#OK !
    return (Xi >> 63) | ((Xi << 1) % K)
    
###### Sous-fonctions de FindDS######
def getY(W0, WC, rot, uX):#OK !
    """Y contains bits [58:58+know_low] of the states.

    >>> uX = [1644185799909158604, 17240438759649590578, 11978199211677924477, 5312027464555861226, 16472746463602985594]
    >>> getY(5066, 5697, [35, 1, 55, 20, 50], uX)
    [475558, 434298, 393054, 445126, 20235]
    """
    Y = []
    LOW = 2**known_low
    for i in range(nbiter):
        a = (powA[i] * W0 + polA[i] * WC) % LOW
        y = (a ^ (uX[i] % LOW)) * 64 + (rot[i] ^ (uX[i] >> 58))
        Y.append(y)
    return Y

def getYprim(Y, WC, W0): #OK ! avec erreurs de retenues ~64bits (polC polW)
    """
    >>> getYprim([334551, 344605, 20585, 258928, 328272], 3153, 3443)
    [334551, 209941, 95413, 359500, 210406]
    """
    Yprim=[(Y[i] - ((polA[i] * WC + powA[i]*W0) >> 58)) % 2**(known_up + known_low) for i in range(nbiter)]
    return Yprim

def getDY(Y, WC, W0): #OK ! avec erreurs de retenues ~64bits (polC polW)
    Yprim = getYprim(Y, WC, W0)
    DY = [(Yprim[i+1] - Yprim[i]) % 2**(known_up + known_low) for i in range(nbiter-1)]
    return DY

######FINDDS######
def FindDS64(uX, rot, W0, WC): #rajouter rot dans la version non test ? #OK! ~64bits
    """
    >>> uX = [4176226900736828747, 9952865681489253366, 2688170825735170541, 1480864413810423309, 12311021560119210858]
    >>> DS64, Y0 = FindDS64(uX, [54, 1, 4, 13, 17], 4995, 1897)
    >>> DS64
    (3143953020030867257, 10643623498242246749, 18015690067871096593, 6081360100074780053)
    >>> Y0
    373304
    >>> (DS64[1] - a * DS64[0]) % (2**64)
    0
    >>> (DS64[3] - a**3 * DS64[0]) % (2**64)
    0
    """
    Y = getY(W0, WC, rot, uX)
    DY = getDY(Y, WC, W0) #OK avec erreurs de retenues!
    tmp = tuple([y << (58 - known_low) for y in DY])#on rajoute les zéros, recentrage impossible à cause des erreurs de retenues
    DS64 = fpylll.CVP.closest_vector(Greduite1, tmp)
    return DS64, Y[0]

######FINDROTI######
#DS64ij = ((polA[j] - polA[i])*DSmod0) % 2**k

def FindRoti(DS640, X, i, Y0, W0, WC): #OK !
    LOW = (1 << known_low)
    DS640i = (polA[i] * DS640) % K
    DSmod0i = ((DS640i << known_low) + W0 * powA[i] + WC * polA[i] - WC - W0) % (1 << (k + known_low))
    # Yi = vraiYi ou vraiYi - 1 à cause de la retenue
    Yi1 = (Y0 + (DSmod0i >> 58)) % (1 << (known_low + 6)) #avec ou sans retenue
    Yi2 = Yi1 + 1
    Wi = (W0 * powA[i] + WC * polA[i]) % LOW
    roti = []
    for j in range(k):
        test1 = (((X ^ (Yi1 >> known_up)) % LOW) == Wi) and ((j ^ (X >> 58)) == Yi1 % 64)
        test2 = (((X ^ (Yi2 >> known_up)) % LOW) == Wi) and ((j ^ (X >> 58)) == Yi2 % 64)
        if test1 or test2:
            roti.append(j)
        X = unrotate1(X)
    return roti


def FindRot(DS640, X, Y0, W0, WC): #OK !
    """
    >>> X = [7562813409073772899, 11417283543806325847, 12362888193759185837, 16903337852395683996, 14856662724370164170, 14353358590961859327, 1709017510603189214, 11112388053225929234, 16326783999071815356, 11265942744832116279, 1663761383601515443, 8264994773702130067, 9504586450162996576, 10071306062678369278, 17269754077915402502, 17284552641819218296, 7862210236782409543, 1173693182357848465, 8669989150777025119, 14483674533785677591, 11817711069368805698, 17652236263878993880, 10952708804218529727, 5546964202247030872, 3787337376045307491, 2463598675764221383, 5606708014410888418, 6875937051451039350, 8452423750630374467, 17961190178884696858, 5340192466979356088, 8233231102133794361, 12880144913541639341, 15509087343987396656, 18273038824920986967, 2620114879769493608, 2973201289203915946, 11220433026578214934, 12347822205906378923, 4161809839587283616]
    >>> FindRot(4512847418431975849, X, 335296, 2791, 6689)
    [[21], [43], [27], [49], [31], [36], [17], [7], [16], [12], [4], [41], [56], [57], [17], [21], [15], [51], [14], [62], [40], [3], [50], [50], [12], [21], [26], [48], [39], [20], [3], [6], [60], [50], [53], [28], [24], [33], [9], [50]]
    """
    tabrot =[]
    for i in range(nboutput):
        tabrot.append(FindRoti(DS640, X[i], i, Y0, W0,WC))
        if len(tabrot[i]) == 0:
            raise ValueError("pas de rotations possibles")
    return tabrot

def findDS(rot, W0, WC): #OK!
    """
    >>> rot = (48, 13, 19, 55, 60, 17, 18, 25, 27, 7, 29, 62, 40, 26, 8, 43, 11, 33, 25, 47, 9, 24, 57, 57, 15, 61, 9, 12, 49, 42, 63, 44, 41, 21, 10, 34, 18, 31, 57, 56)
    >>> DS = findDS(rot, 4975, 4063)
    >>> DS
    (-3911495036692940322621796134291199, -7279435059422484119362302816873659, -5980753107089600906205976366383207, 26501861599381476532331924531931709, -31630035624456443253911377653383055, 29213102710206694601105260826063989, -23281560453589451868529477218643575, 17888517654251431755244562532048877, 3228901960492265957430914842941665, 2313333816108045922019596861352613, -38216913742480012196731378496271239, 24068986573173657396568329560531613, -1898490002137598832998023266815407, -12376784405299621377499477592545835, 29039659037871761656409287413276777, -11120000896287372100511132327173555, -17733553412493159685771337165456191, -8741721083865437369382845282560507, 30558047682244154952554650584415577, 1357273502023295372192253817507581, -24843031240572301260463448405579727, -3299110994710694014154259362174155, -1487293042488303561252348415241399, -3759853676895093460585569200088915, 36340912119809025232203444281812129, -36471259272312409369956954296904347, 6259392443418117937466593286640185, 33078586047867229586889558846467933, -39185067306776146290022302686584303, 477254915511356870841557501691029, 14118492461907402245240001077094953, -1960793072427708660280495620799731, -9137193975787581074734956369964927, 14813477132465570555338375023945925, -18712543965827835408441644622214375, 14470108143760423276150009776702397, -789275610922086810087986631368719, 2450579317286955048107917366440437, -12618929602727743299982987170192119)
    >>> (DS[1] - a * DS[0]) % (2**(128 - known_low))
    0
    >>> (DS[5] - a**5 * DS[0]) % (2**(128 - known_low))
    0
    """
    rotprim = []
    for i in range(nboutput):
        rotprim.append((rot[i] - ((powA[i] * W0 + polA[i] * WC) >> 122)) % 64)
    tmp = tuple([(rotprim[i+1] - rotprim[i]) << (122 - known_low) for i in range(nboutput - 1)])
    return fpylll.CVP.closest_vector(Greduite2, tmp)
    
def findDSdebug(rot, W0, WC): #OK!
    """
    >>> rot = (36, 11, 26, 47, 27, 43, 40, 4, 6, 14, 46, 7, 26, 15, 29, 52, 5, 27, 59, 51, 50, 10, 18, 12, 20, 13, 41, 11, 33, 29, 11, 34, 57, 43, 57, 9, 4, 44, 5, 61)
    >>> DS = findDSdebug(rot, 931, 4565)
    >>> DS
    (-2402491240055656244800996055138065, -12494820703177916317385412874753429, 16338193271580590973761587180847319, -15300585971148004077979215350754317, 21684776253948045125360273100268159, -27815550822220705053670430230538693, 11625224150341610375324286543409127, 17844093130120930863081999339086659, -380944688746118817409004598549745, -28522485004078932652170383056544501, 520732373288502378058486581798391, -3634819024663226032969760811575917, 2969698435142367696300086782333599, 28037325802964017993093991345016539, -10910817629932953951556521090665721, -5065536280883514534631998922589469, 22910033460594640618280890256081199, -2053703208298283961366229680604245, -29057244128546792360653529929160937, -3838864775261833707985232041974989, 30542591494432091311532645601492671, 422636935015451582314011290366843, -3834374215185520065523022262513113, -31730206033684723762666891971912061, 1605231086249882265091899797936975, 28745614588815424753948783543535179, 206839261025612279689680544465975, 5000213623392168998536915729513683, -13303445174387063160344871468669217, 16214267541652108868690966934029339, -36078096034273452443520370613514937, 37999702462501891126733863140871715, -13428983410411412751530034747362961, 5641457792955538248534388118420715, -29350387749398409235158943575343785, 21982279210625596709854366303937139, -12448157819225143628043984315675905, -11894269968695156406609582729527109, 8505931709653315845144709262926951)
    >>> (DS[1] - a * DS[0]) % (2**(128 - known_low))
    0
    >>> (DS[5] - a**5 * DS[0]) % (2**(128 - known_low))
    0
    """
    rotprim = []
    for i in range(nboutput):
        rotprim.append((rot[i] - ((powA[i] * W0 + polA[i] * WC) >> 122)) % 64)
    # différence ici !
    rotprim[1] = (rotprim[1] - 1) % 64
    # fin différence
    tmp = tuple([(rotprim[i+1] - rotprim[i]) << (122 - known_low) for i in range(nboutput - 1)])
    return fpylll.CVP.closest_vector(Greduite2, tmp)

if __name__ == '__main__':
    cpt = {}
    cpt[False, False] = 0
    cpt[False, True] = 0
    cpt[True,  False] = 0
    cpt[True,  True] = 0
    
    cptrotfail = 0
    n = 1000
    total_output_size = 0
    for _ in range(n):
      try:
        X, S, c = sortiesGenerateur()    
    
        # génère ce qu'on est pas censé connaître, et que le gros calcul doit retrouver
        W0 = S[0] % LOW
        WC = c % LOW
        rot = []
        for i in range(nboutput):
            rot.append(S[i] >> 122)
    
        uX = unrotateX(X,rot)    
        DS64, Y0 = FindDS64(uX, rot, W0, WC)#OK!

         # refait la partie qui trouve DS64.
#        Y = [((S[i] - (powA[i] * W0 + polA[i] * WC)) >> known_low) % (1 << (128 - known_low)) for i in range(nboutput)]
#        cheat_DS = [(Y[i + 1] - Y[i]) % (1 << (128 - known_low)) for i in range(nboutput - 1)]
#        if cheat_DS[0] % (1 << 64) != DS64[0]:
#            raise ValueError("DS64 incorrect ; cannot proceed further")

        # calcule toutes les rotations futures à partir de DS64
        tabrot = FindRot(DS64[0], X, Y0, W0, WC) #a l'air OK! 
        listrot = list(product(*tabrot))
       
        test_a = False
        test_b = False
    
        output = set()
    
        # vérifie qu'on a bien trouvé un des bon DS complet.
        listDS = [findDS(rot, W0, WC) for rot in listrot]
        output.update(listDS)
        for DS in listDS:
            Sprim = [(S[i] - polA[i] * WC - powA[i] * W0) % N for i in range(nboutput)]
            if DS[0] == ((Sprim[1] - Sprim[0]) >> known_low):
                test_a = True

        listDS = [findDSdebug(rot, W0, WC) for rot in listrot]
        output.update(listDS)
        for DS in listDS:
            Sprim = [(S[i] - polA[i] * WC - powA[i] * W0) % N for i in range(nboutput)]
            if DS[0] == ((Sprim[1] - Sprim[0]) >> known_low):
                test_b = True
        
        cpt[test_a, test_b] += 1
        total_output_size += len(output)

        if not (test_a or test_b):
            raise TypeError("both approaches failed")

      except ValueError as e:
        print(repr(e))

    print(n)
    print("cpt:")
    print(cpt)
    print("total output size = {}".format(total_output_size))

