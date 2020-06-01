import fpylll
from math import sqrt

from g6k import Siever
from g6k.algorithms.pump import pump
from g6k.utils.stats import dummy_tracer

a = 47026247687942121848144207491837523525

def getG(n, mod):
    G = fpylll.IntegerMatrix(n, n)
    for i in range(0, n):
        G[0, i] = pow(a, i, mod)
    for i in range(1, n):
        G[i, i] = mod
    fpylll.LLL.reduction(G) # in place
    return G

def svp_fpylll(G):
    """
    Returns the shortest vector of the lattice spanned by the rows of G.
    """
    x = fpylll.SVP.shortest_vector(G, method = "fast", flags=fpylll.SVP.VERBOSE)
    return x


def svp_g6k(G):
    """
    Returns the shortest vector of the lattice spanned by the rows of G.
    """
    g6k = Siever(G)
    pump(g6k, dummy_tracer, 0, G.nrows, 0)
    return G[0]

def norm(v):
    return sqrt(sum(x**2 for x in v))
    
def show_targets(ell):
    for n in range(30, 70):
        target = 2 * sqrt(n) * (2**(123-ell) - 1)
        print("{}: {}".format(n, target))

def min_m_CVP(ell):
    for n in range(60, 64, 1):
        G = getG(n, 2**(128-ell))
        shortest = svp_g6k(G)
        lambda1 = norm(shortest)
        target = 2 * sqrt(n) * (2**(123 - ell) - 1)
        print("{}: {} vs {}".format(n, target, lambda1))
        if target <= lambda1:
            break
    print("SUCCESS")
    return n

# k = 128
## n = 42 -> 7.10958185672e+37
## n = 48 -> 9.40694696471e+37
## n = 54 -> 1.23309364376e+38
## n = 60 -> 1.48734754252e+38
## n = 

if __name__ == '__main__':
    show_targets(0)
    #min_m_CVP(0)

    n = 64
    G = getG(n, 2**128)
    print(sum([norm(x) for x in G])/n)
    #baseline = norm(G[0])
    g6k = Siever(G)
    #print("baseline = {}".format(baseline))
    for i in range(64):
        pump(g6k, dummy_tracer, 0, n, 0)
        print([norm(x) for x in G[0:8]])
    
    #print(G[0])

    #shortest = G[0]
    #print("pump result = {}".format(norm(shortest)))
    #for i, sqnorm, coeffs in g6k.best_lifts():
    #    v = G.multiply_left(coeffs)
    #    if norm(v) < baseline:
    #        baseline = norm(v)
    #        print("**** new best from best_lifts ***** --> {}".format(baseline))
    #    for coeffs in g6k.itervalues():
    #        v = G.multiply_left(coeffs)
    #        if norm(v) < baseline:
    #            baseline = norm(v)
    #            print("**** new best from itervalues ***** --> {}".format(baseline))
    #        
#