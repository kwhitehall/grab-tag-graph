#To do a really quick and dirty memory check, run
#python baselineTiming.py  & top -d 00.2 -b -n 300 -p $! >  usage.txt
#Then you'll have to interrupt (I don't know why this never terminates). Then run
#python getMemUsage.py

import numpy as np
import matplotlib.pyplot as plt

f = open('usage.txt','r')
text = f.read()
text = text.split()

mems = np.array([])
spot = 80

while spot<len(text):
    try:
        mems = np.append(mems,int(text[spot]))
    except:
        break;
    spot += 87
    
#print(mems)
f.close()

x = np.arange(0,mems.shape[0],1)
plt.plot(x,mems/1000)
plt.show()
