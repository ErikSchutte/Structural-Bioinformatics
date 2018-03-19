import random as ran

fw = open('ran.dat','w')
for i in range(1000000):
    fw.write(str(i)+' '+str(ran.random())+'\n')
#    print ran.random()
fw.close()