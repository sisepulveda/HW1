from numpy import array, zeros, ix_
from beam_element import beam_element

xy = array([
    [0,0],
    [0,3],
    [0,5],
    [0,6],
    [6,0],
    [6,3],
    [6,5],
    [6,6.5]
    ])

conec = array([
    [0,1],
    [1,2],
    [2,3],
    [4,5],
    [5,6],
    [6,7],
    [1,5],
    [2,6],
    [3,7]
    ],
    dtype=int
    )

beam = 0.2*0.4  #m2
column = 0.3*0.3    #m2
roof = 0.2*0.2  #m2
D = 2400    #kg/m3

properties_0 = {}   # beams
properties_0["E"] = 2400000000.
properties_0["A"] = beam
properties_0["I"] = ((0.4**3)*0.2)/12
properties_0["qx"] = 0.
properties_0["qy"] = -D*beam

properties_1 = {}   #columns
properties_1["E"] = 2400000000.
properties_1["A"] = column
properties_1["I"] = ((0.3**3)*0.3)/12
properties_1["qx"] = 0.
properties_1["qy"] = -D*column  #carga distribuida vertical

properties_2 = {}   #roof
properties_2["E"] = 2400000000.
properties_2["A"] = roof
properties_2["I"] = ((0.2**3)*0.2)/12
properties_2["qx"] = 0.
properties_2["qy"] = -D*roof #carga distribuida vertical


properties = [properties_1, properties_1, properties_1, properties_1, properties_1, properties_1,
              properties_0, properties_0, properties_2]


Nnodes = xy.shape[0]
Nelems = conec.shape[0]

NDOFs_per_node = 3
NDOFs = Nnodes*NDOFs_per_node

K = zeros((NDOFs, NDOFs))
f = zeros((NDOFs, 1))

for e in range(Nelems):
    ni = conec[e,0]
    nj = conec[e,1]

    #print(f"e = {e} ni = {ni} nj = {nj}")

    xy_e = xy[[ni, nj], :]

    ke, fe = beam_element(xy_e, properties[e])

    #print(f"e = {e} ke = {ke}")

    # Node k ---> [ 3*k, 3*k+1, 3+k+2 ]

    d = [3*ni, 3*ni+1, 3*ni+2, 3*nj, 3*nj+1, 3*nj+2]    #global DOFS from local DOFS

    #Direct stiffness method
    for i in range(2*NDOFs_per_node):
        p = d[i]
        for j in range(2*NDOFs_per_node):
            q = d[j]
            K[p, q] += ke[i,j]
        f[p] += fe[i]

#Hand calculation of f

Lb = 6.0
f[4] = -D*beam*Lb/2
f[5] = -D*beam*Lb**2/12
f[7] = -D*beam*Lb/2
f[8] = -D*beam*Lb**2/12
f[10] = -D*beam*Lb/2
f[11] = -D*beam*Lb**2/12
f[16] = -D*beam*Lb/2
f[17] = D*beam*Lb**2/12
f[19] = -D*beam*Lb/2
f[20] = D*beam*Lb**2/12
f[22] = -D*beam*Lb/2
f[23] = D*beam*Lb**2/12

#0,1,2      #1
#3,4,5      #2
#6,7,8      #3
#9,10,11    #4
#12,13,14   #5
#15,16,17   #6
#18,19,20   #7
#21,22,23   #8

#print(K)
#print(f)



# System partitioning and solution

free_DOFs = [3, 4, 5, 6, 7, 8, 9, 10, 11, 12,
             13, 14, 15, 16, 17, 18, 19, 20, 21, 22, 23]
constrained_DOFs = [0, 1, 2, 12]

Kff = K[ix_(free_DOFs, free_DOFs)]
Kfc = K[ix_(free_DOFs, constrained_DOFs)]
Kcf = K[ix_(constrained_DOFs, free_DOFs)]
Kcc = K[ix_(constrained_DOFs, constrained_DOFs)]

print(K)

ff = f[free_DOFs]
fc = f[constrained_DOFs]

# Solve
from scipy.linalg import solve
u = zeros((NDOFs, 1))

u[free_DOFs] = solve(Kff, ff)

#print(u)

#Get reaction forces
R = Kcf @ u[free_DOFs] + Kcc @ u[constrained_DOFs] - fc

print(f"u = {u}")
print(f"R = {R}")