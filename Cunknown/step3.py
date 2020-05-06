
k = 64
known_up  = 6
known_low  = 13
a = 2549297995355413924 * 2**64 + 4865540595714422341
a_inv = 566787436162029664* 2**64 + 11001107174925446285
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
    X[ 0] = 0x46b9d66bb1e9fbf0;
    X[ 1] = 0xc0b7382496f0e363;
    X[ 2] = 0x78bf5fb5f4bbc09e;
    X[ 3] = 0x75f352cb784fa167;
    X[ 4] = 0x0c8c9a24d94bbc61;
    X[ 5] = 0xa665faf0ba5eede6;
    X[ 6] = 0x14162574ccdcdc3f;
    X[ 7] = 0x5fae9c165a0f7180;
    X[ 8] = 0x780f6c5c0f9bf744;
    X[ 9] = 0xee89349d3a7cd604;
    X[10] = 0xa2052341850dbc24;
    X[11] = 0x99b08d5d32d5048d;
    X[12] = 0x2a854a252e8284ea;
    X[13] = 0xe2ffa2aab19b966c;
    X[14] = 0xf18760da7640f5da;
    X[15] = 0x89e7da73a09f2ec9;
    X[16] = 0x71044040d182a2b2;
    X[17] = 0x21e4bac7c9d97896;
    X[18] = 0xe077b5e0a770055f;
    X[19] = 0x3c5d2adc91377527;
    X[20] = 0x8d1f3c155a55ca0b;
    X[21] = 0x8db1105671cc6031;
    X[22] = 0x4cd45a248ea6912d;
    X[23] = 0x9af4e51b0af05f0d;
    X[24] = 0xa400b754f2b2b0b3;
    X[25] = 0xc493c31ea136a8f2;
    X[26] = 0x96672aaa1688b031;
    X[27] = 0xd2f7f9527275bfd2;
    X[28] = 0x5187489c1b75e2a7;
    X[29] = 0x627224615c897f12;
    X[30] = 0xa202366a1039e4bc;
    X[31] = 0x48fc42d00ed0c6ed;
    X[32] = 0xc66eb8c29ef48de2;
    X[33] = 0xe0ac6abe677dbdb7;
    X[34] = 0x1dae3a9ed33b97fe;
    X[35] = 0x4359d58482eb24c1;
    X[36] = 0x075b7c7245fddda5;
    X[37] = 0xf9d67dfd4bdb5bdd;
    X[38] = 0x18c778e280d64a55;
    X[39] = 0x7534aeccf7b47c2b;
    X[40] = 0x2937665c341e30c5;
    X[41] = 0x190f02e1acb9233e;
    X[42] = 0x2895523d59cbfdb1;
    X[43] = 0x9a730f18211b67d7;
    X[44] = 0x72277c4b0673a1eb;
    X[45] = 0xda17af3b16f6f009;
    X[46] = 0x7c245c9bdd4a47b7;
    X[47] = 0x20cb56ac9e64fc94;


    original_X = X   

    W_0, W_c, start_DS, rot = 0x00c9, 0x0643, 0x15b3f74a3b47132dcf2389239ebaa9a7, (49, 54, 11, 48, 41, 8, 18, 50, 29, 40, 29, 50, 26, 22, 63, 55, 24, 62, 60, 46, 32, 52, 17, 1, 47, 23, 35, 49, 18, 3, 6, 30, 0, 63, 16, 25, 63, 53, 63, 59)


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

    S_prev = (found_S0 - found_c) * a_inv % N
    seed = (S_prev - found_c - a*found_c) * a_inv % N

    print("c = {:016x} {:016x}".format(found_c >> 65, (found_c >> 1) % K))
    print("seed = {:016x} {:016x}".format(seed >> 64, seed % K))

    XX = sortiesGenerateur(found_S0, found_c, 64)
    for i, x in enumerate(XX):
        print("X[{:2d}] = 0x{:016x};".format(i, x))
