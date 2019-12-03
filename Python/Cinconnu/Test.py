from fonctions import *

X,S,c = sortiesGenerateur()
rot = [S[i] // 2**(k * 2 - known_up) for i in range(nbiter)]
uX = unrotateX(X, rot)
print("uX")
print(uX)

W0 = S[0] % 2**known_low
print(W0)
WC = c % 2**known_low

Y = getY(W0, WC, rot, uX)
print("Y")
print(Y)
Yprim = getYprim(Y, WC, W0)
print("Yprim")
print(Yprim)
DY = getDY(Y, WC, W0)
print("DY")
print(DY)

DS64 = FindDS64(uX, rot, W0,WC)
print("DS64")
print(DS64)
rot2 = FindRot(DS64[0], X, Y[0], W0, WC, 20)
print(rot)
print(rot2)