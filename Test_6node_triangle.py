# -*- coding: utf-8 -*-
"""
Created on Sun Oct 21 16:38:14 2018

@author: bjohau
"""

import numpy as np
import triangles_with_TODO as tri

cyclic_ijk = [0,1,2,0,1]      # Cyclic permutation of the nodes i,j,k

# ----- Topology -------------------------------------------------
ex = np.array([0.0,1.0,0.0, 0.5,0.5,0.0])
ey = np.array([0.0,0.0,1.0, 0.0,0.5,0.5])

for i in range(3):
    j = cyclic_ijk[i+1]
    k = cyclic_ijk[i+2]
    ex[i+3] = (ex[i] + ex[j])/2
    ey[i+3] = (ey[i] + ey[j])/2

th = 0.1
ep = [1,th]

E  = 2.1e11
nu = 0.3

D = np.array([
        [ 1.0,  nu,  0.],
        [  nu, 1.0,  0.],
        [  0.,  0., (1.0-nu)/2.0]]) * E/(1.0-nu**2)

eq = [1.0, 3.0]

#Ke, fe = tri.plante(ex,ey,ep,D,eq)

Ke = np.zeros((12,12))
fe = np.zeros((12,1))



rigX = np.zeros((12,1))
rigY = np.zeros((12,1))
rigR = np.zeros((12,1))



for i in range(6):
    rigX[i*2  ,0] = 1.0
    rigY[i*2+1,0] = 1.0
    rigR[i*2  ,0] = ey[i]
    rigR[i*2+1,0] = -ex[i]

Ke, fe = tri.tri6_Kmatrix(ex,ey,D,th,eq)

print('Stiffness matrix:\n', Ke)
print('Consistent forces:\n', fe)

fx = Ke @ rigX
fy = Ke @ rigY
fr = Ke @ rigR

print('Force from rigX translation:\n',fx)
print('Force from rigY translation:\n',fy)
print('Force from rigR rotation:\n',fr)


constEx = np.zeros((12,1))
constEy = np.zeros((12,1))
constGamma1 = np.zeros((12,1))
constGamma2 = np.zeros((12,1))


for i in range(6):
    constEx[i*2  ,0] = ex[i]
    constEy[i*2+1,0] = ey[i]
    constGamma1[i*2  ,0] = ey[i]
    constGamma2[i*2+1,0] = ex[i]
    
zetaInt = np.array([[0.5,0.5,0.0],
                    [0.0,0.5,0.5],
                    [0.5,0.0,0.5]]) 
    
for i in range(3):
    print('zeta', zetaInt[i])
    nx, ny = tri.tri6_shape_function_partials_x_and_y(zetaInt[i],ex,ey)
    print('nx', nx)
    print('ny', ny)

    Be = tri.tri6_Bmatrix(zetaInt[i],ex,ey)
    
    print('Be\n', Be)
    
    Ex = Be @ constEx
    Ey = Be @ constEy
    G1  = Be @ constGamma1
    G2  = Be @ constGamma2
    print('Ex:\n',Ex)
    print('Ey:\n',Ey)
    print('G:\n',G1)
    print('G:\n',G2)

    
