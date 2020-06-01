import re

pattern = re.compile(r' \[(?P<time>[0-9]+\.[0-9])s\]')

match = pattern.search('Done task 53a (W_0=0001 / W_c=0275) [41.6s]')
assert match
assert match['time'] == '41.6'

T = []
with open("breaker1352922.out") as f:
    for line in f:
        match = pattern.search(line) 
        if match:
            T.append(float(match.group('time')))

d = {}
for t in T:
    x = round(10 * t)
    if x not in d:
        d[x] = 1
    else:
        d[x] += 1

for t in d:
    print("{}: {}".format(t, d[t]))
