a = 2549297995355413924 * 2**64 + 4865540595714422341

def mat(n,M,) :

    future_matrice = [[1]+[0]*(n-1)]
    for i in range(1,n):
        lignei=[0]*n
        lignei[0]=a^i % M
        lignei[i] = M
        future_matrice += [lignei]
    future_matrice = matrix(future_matrice).transpose()
    matrice = future_matrice.LLL()
    return(matrice)
    
def normeG(n,M):

    matrice = mat(n,M)
    
    matrice=matrice.BKZ(block_size=20)
    short=vector(matrice[0,:]) 
    normv=norm(short)
    normG = norm(matrice) 
    normGm = norm(matrice^(-1))
    return(normG*normGm, normv)
    
def min_l_babai(n,M):
    cond,normv = normeG(n,M)
    l=1
    truc = (1+cond)*sqrt(n)*(2^(65-l)-1)
    while (truc >= normv):
        l=l+1
        truc = (1+cond)*sqrt(n)*(2^(65-l)-1)
    return(l)
    
def min_m_CVP():
    n=1
    delta = 1
    normv = 0
    while 2*delta >= normv :
        n=n+1
        delta = sqrt(n)*2**(128-6)
        cond,normv = normeG(n,2**(128))
    return(n)

