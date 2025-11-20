# -*- coding: utf-8 -*-
"""
Created on Sun Oct 21 16:38:14 2018

@author: bjohau
"""
import numpy as np
import fem_utilities as fem_util
import fem_models
import copy

# Element Type
numElementNodes = 4  # Valid numbers 3, 4, 6, 9

# Number of nodes: Should be odd numbers in order to handle 9 node quad and 6 node triangle
numNodesX = 15
numNodesY = 9

# Cantilever with dimensions H x L x thickness
H         =  2.0
L         = 10.0
thickness =  0.1

#End load, Given as resultant
endLoadXY = np.array([0.0,3.0e6])
#endLoadXY = np.array([3.0e6,0])
#endLoadXY = np.array([4.2e9,0.0]) # Should give unit disp at Poisson = 0


# Distributed load in x and y, load pr unit area
eq = np.array([0.,1.0e3])
eq = np.array([0.,0.])
eq = None  # No distributed load

model = fem_models.CantileverModel(L, H, numElementNodes, numNodesX, numNodesY, endLoadXY, eq)

model.vtu_write_connected_mesh("GeometryOnly_ConnectedMesh.vtu")

# Assemble stiffness matrix
K, R = model.assemble_stiffness_matrix_and_load_vector()

# Solve the system of equations
r, R0 = fem_util.solveq(K, R, model.bcdofs)

# Compute element corner stresses
elementCornerStresses = model.compute_element_corner_stresses(r)

nodMiddle = numNodesY//2 +1  # Mid nod on right edge
xC = r[-(nodMiddle*2)  ,0] # 2 dofs per node, so this is the middle dof on end
yC = r[-(nodMiddle*2)+1,0] # 2 dofs per node, so this is the middle dof on end
print("Displacement center node right end,  x:{:12.3e}   y:{:12.3e}".format(xC, yC))

# Sum uf reaction forces
R0Sum = np.zeros(2,'f')
ndofs = len(R0)
#for i in range(ndofs,2):
for i in range(0,numNodesY*2,2):
    R0Sum[0] += R0[i  ,0]
    R0Sum[1] += R0[i+1,0]
print("Total reaction force: x:{:12.3e} y:{:12.3e})".format(R0Sum[0],R0Sum[1]))

# Draw the displacements and stresses
model.vtu_write_stl_style_mesh("Results.vtu",dispVector=r,elementCornerStresses=elementCornerStresses)

