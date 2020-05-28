#Predictor input --- 3 consecutive calls to lehmer64()
X = [0x1e1495fb7df7b3cf, 0xa267765f13b93bf9, 0xe5b3241c4f81c9a8]

# Prediction algorithm
a = 0xda942042e4dd58b5
r = round(2.64929081169728e-7 * X[0] + 3.51729342107376e-7 * X[1] + 3.89110109147656e-8 * X[2])
s = round(3.12752538137199e-7 * X[0] - 1.00664345453760e-7 * X[1] - 2.16685184476959e-7 * X[2])
t = round(3.54263598631140e-8 * X[0] - 2.05535734808162e-7 * X[1] + 2.73269247090513e-7 * X[2])
u = r * 1556524 + s * 2249380 + t * 1561981
v = r * 8429177212358078682 + s * 4111469003616164778 + t * 3562247178301810180
w = (a**3 * (a*u + v)) >> 64
x = w & 0xffffffffffffffff

# Predictor output --- what the fourth call to lehmer64() would produce
print(hex(x))
