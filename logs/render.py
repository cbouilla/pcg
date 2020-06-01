from matplotlib import pyplot as plt

d = {41.4: 8405, 
     41.5: 7858, 
     41.6: 574570, 
     41.7: 92711, 
     41.8: 349415, 
     41.9: 11835, 
     42.0: 2763, 
     42.1: 8848, 
     42.2: 644201, 
     42.3: 4384, 
     42.4: 13976, 
     42.5: 133008, 
     42.6: 241950, 
     42.7: 1324, 
     43.4: 491, 
     54.7: 364, 
     54.8: 531, 
     43.6: 50, 
     54.9: 54, 
     43.5: 403, 
     43.3: 5, 
     55.0: 2, 
     43.7: 2, 
     40.0: 1, 
     43.8: 1
}

X = sorted(d.keys())
Y = [d[k] for k in X]

fig, ax = plt.subplots(1, 1, figsize=(16, 10))
ax.bar(X, Y, color='violet', width=0.08)
    
ax.set_yscale('log', basey=2)
#ax.axis(xmin=0, xmax=npasses+1, ymin=ymin, ymax=N*1.15)
#ax.set_yticks(yticks)
ax.set_xticks([i for i in range(39, 57)])
fig.set_tight_layout(True)
#ax.legend(loc='best')
ax.set_ylabel('# Tasks')
ax.set_xlabel('Time per Task (s)')
#ax.yaxis.set_ticks_position('both')

plt.show()
#fig.savefig(outfile)
