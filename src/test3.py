"""
Test finite difference using the TEP, for a specific case
"""

import numpy as np

def calcF(atomIndices,cartIndices,disps,ndisps):
    u = np.zeros([natoms,3])
    for d in range(0,ndisps):
        u[atomIndices[d+1]][cartIndices[d+1]] = disps[d]
    #print(u)
    f = np.zeros([natoms,3])
    for i in range(0,natoms):
        for a in range(0,3):
            for j in range(0,natoms):
                for b in range(0,3):
                    f[i][a] -= fc2[i][j][a][b]*u[j][b]
                    for k in range(0,natoms):
                        for c in range(0,3):
                            f[i][a] -= fc3[i][j][k][a][b][c]*u[j][b]*u[k][c]

    return f[atomIndices[0]][cartIndices[0]]

natoms = 8

# Read FC2
fh = open("FC2", 'r')
line = fh.readline()
fc2 = np.zeros([natoms,natoms,3,3])
for i in range(0,natoms):
    for a in range(0,3):
        for j in range(0,natoms):
            for b in range(0,3):
                line = fh.readline().split()
                #i = int(line[0])
                #a = int(line[1])
                #j = int(line[2])
                #b = int(line[3])
                val = float(line[4])
                #print("%d %d %d %d" % (i,j,a,b))
                fc2[i][j][a][b] = val
fh.close()

# Read FC3
fh = open("FC3", 'r')
line = fh.readline()
fc3 = np.zeros([natoms,natoms,natoms,3,3,3])
for i in range(0,natoms):
    for a in range(0,3):
        for j in range(0,natoms):
            for b in range(0,3):
                for k in range(0,natoms):
                    for c in range(0,3):
                        line = fh.readline().split()
                        #i = int(line[0])
                        #a = int(line[1])
                        #j = int(line[2])
                        #b = int(line[3])
                        #k = int(line[4])
                        #c = int(line[5])
                        val = float(line[6])
                        fc3[i][j][k][a][b][c] = val
fh.close()

atomIndices = [0,0,4]
cartIndices = [0,0,1]
#atomIndices = [0,0]
#cartIndices = [0,0]
hvec = [1e-10,1e-9,1e-8,1e-7,1e-6,1e-5,1e-4,1e-3,1e-2,1e-2]
for h in hvec:

    
    """
    disps = h*np.array([1.,1.])
    f1 = calcF(atomIndices,cartIndices,disps,2)
    disps = h*np.array([1.,-1.])
    f2 = calcF(atomIndices,cartIndices,disps,2)
    disps = h*np.array([-1.,1.])
    f3 = calcF(atomIndices,cartIndices,disps,2)
    disps = h*np.array([-1.,-1.])
    f4 = calcF(atomIndices,cartIndices,disps,2)

    d2 = (-1./(0.01**2))*(f1-f2-f3+f3)
    #print(" %f %f %f %f" % (f1,f2,f3,f4))
    print(" h,d2: %e, %e" % (h,d2))
    """
    

    """
    disps = h*np.array([1.])
    f1 = calcF(atomIndices,cartIndices,disps,1)
    disps = h*np.array([-1.])
    f2 = calcF(atomIndices,cartIndices,disps,1)

    d2 = (-1./(0.01**2))*(f1+f2)
    #print(" %f %f %f %f" % (f1,f2,f3,f4))
    print(" d2: %f" % (d2))
    """

    """
    disps = h*np.array([1.])
    f1 = calcF(atomIndices,cartIndices,disps,1)
    disps = h*np.array([-1.])
    f2 = calcF(atomIndices,cartIndices,disps,1)

    d1 = (-1./(2*h))*(f1-f2)
    #print(" %f %f %f %f" % (f1,f2,f3,f4))

    print(" d1: %.20e" % (d1))
    """
