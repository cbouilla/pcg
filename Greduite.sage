##### NE PAS OUBLIER DE RETIRER LA DERNIERE VIRGULE !!!


k = 64
m = 2**(2 * k)
a = 2549297995355413924 * 2**64 + 4865540595714422341

def getA(n):
    A= [[m]]
    for i in range(n-1):
        A[0].append(0)
    powA = a % m
    for i in range(n-1):
        A.append([powA])
        powA = powA * a % m
        for j in range(n-1):
            if(i == j):
                A[i+1].append(-1)
            else:
                A[i+1].append(0)
    return matrix(A)
    
def getG(n):
    G = [[1]]
    for i in range(n-1):
        G[0].append(0)
    powA = a % m
    for i in range(n-1):
        G.append([powA])
        powA = powA * a % m
        for j in range(n-1):
            if(i == j):
                G[i+1].append(m)
            else:
                G[i+1].append(0)
    return matrix(G)
        
def getGreduite(n):
    G = getG(n)
    return G.transpose().LLL().transpose()
    
def getInvG(Greduite):
    return Greduite.inverse().n()

def printFic(n):
    G = getGreduite(n)
    invG = getInvG(G)
    fic = open("G", "a")
    fic.write("\n\n#nbiter = %s\n\n Greduite = {" % n)
    for i in range(n):
        fic.write("\n")
        for j in range(n):
            fic.write("%s, " % G[i][j])
    fic.write("};\n\ninvG = {")
    for i in range(n):
        fic.write("\n")
        for j in range(n):
            fic.write("%s, " % invG[i][j])
    fic.write("};\n\n")
    fic.close()
            
'''A = matrix([[m, 0,0,0],
           [a%m,-1,0,0],
           [(a*a)%m,0,-1,0],
           [(a*a*a)%m,0,0,-1]])
A = matrix([[m, 0,0],
           [a%m,-1,0],
           [(a*a)%m, 0,-1]])


G = matrix([[1, 0, 0, 0,0],
     [a%m, m, 0, 0,0],
     [(a*a)%m, 0, m, 0,0],
     [(a*a*a)%m, 0, 0, m,0],
     [(a*a*a*a)%m, 0, 0, 0, m]])
G = matrix([[1, 0,0,0,0],
     [a%m, m,0,0,1],
     [(a*a)%m,0, m, 0, a+1],
     [(a*a*a)%m, 0, 0,m,a*a+a+1],
     [0,0,0,0,1]])'''
     
#GG = G.transpose().LLL()
#invGG = GG.inverse().n()
printFic(25)

