# this SageMath program produces a good basis for
# the "large" lattice in dimension 64, modulo 2**128;
# The shortest vector has norm 2^127.03...

a = 47026247687942121848144207491837523525

def getG(n, mod):
    """
    Naive basis
    """
    G = matrix(ZZ, n, n)
    for i in range(0, n):
        G[0, i] = pow(a, i, mod)
    for i in range(1, n):
        G[i, i] = mod
    return G

G = getG(64, 2**128)

L = [
 12144252875850345479015002205241987363,
 -15348615718011256064636590606166444173,
 -9084835538682141851035496395127502261,
 17234466561204322281445746373216450637,
 -30709520135048357458223997075862325767,
 13273365585194684751192188342361392061,
 -16657707705878774883995565342082094357,
 10551835762264010302610847026218529611,
 9741843352863159575387669322076694498,
 30401545415842458720824638334153655996,
 -21301422832830108630138234390597720063,
 12943082679906251720732751543854447750,
 -26608736029397995753447945248664058244,
 16504008779817367542456313016700742074,
 -10623665749400971081452814639205027657,
 -42221955640534425567937724106639808954,
 30891705174976345231868365741077473210,
 36583589311333508006721438259000555949,
 -7050528015382194186611161130462441554,
 29023698504320799821125567462694089388,
 -28780651977143078032165325525310692006,
 80641828999898464226011186926508194344,
 -56574857747725737575017746952734912354,
 69356676566431307205578324735925446877,
 -33869683024546744201772555262569839986,
 -26565235055559932529672543662504791838,
 12570516463211853078121725511001272792,
 -3545191593186327315454855798676360051,
 34483716896778783922111688009948767715,
 -20939456075167816323854415536161525816,
 37054527231762903007893323872329597868,
 -65113113872795261158497215763575544978,
 -27970657614360739023057572629108089923,
 -9380443740630786183893558576010683559,
 -31684372116213316987059381973832378607,
 6514548643534369904269689195712545814,
 14847729039548998339264458667105218977,
 20177225057138896843222830352311122871,
 -10200019003678439599482650259589205069,
 10164570837071212387748750871984064078,
 -16485770795775188591072402663394972722,
 -10182226575367101730098853924492492198,
 -21298248172336954279037664114356234573,
 -2151417390168764172560507076988402637,
 20588737496290684579436443152586408051,
 25195201573991617845656823418508896439,
 -9250559638374831564761683687174364975,
 24609591485130417852256436343428786816,
 34329819098992409598701494495874916561,
 21067050821963084044299000164807345041,
 -6187356881994741304068539452724394683,
 82300199668589409058661808102592212051,
 32776452661595718098626668848975457651,
 -20144177049611899844263802898968380329,
 7641717196898404580174495749776072261,
 4540641926702911772918893897827638233,
 -26129828811333947870155679432881913189,
 37836842229341883768177535380751590544,
 28364481592155717163003341291316712226,
 -25569624066399106381499263355295699696,
 -6060757243887189832616629449657156176,
 27944345301214874985591522626368184967,
 87006401896471707412013173253053669642,
 -4241521124916816484733669417798364034
]

def center(u):
    v = u % 2**128
    for i in range(64):
        if v[i] > 2**127:
            v[i] -= 2**128
    return v

# good basis, nearly HKZ reduced. First vector is SVP.
GG = matrix([center(x*G[0]) for x in L])

