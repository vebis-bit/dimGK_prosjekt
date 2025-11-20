# -*- coding: utf-8 -*-
"""
Created on Sun Oct 21 08:15:51 2018

@author: bjohau
"""


import numpy as np

def tri3_area(ex, ey):
    """
    Compute the area of a triangle element

    :param list ex: element x coordinates [x1, x2, x3]
    :param list ey: element y coordinates [y1, y2, y3]
    :return area of the triangle
    """

    tmp = np.matrix([[1, ex[0], ey[0]],
                     [1, ex[1], ey[1]],
                     [1, ex[2], ey[2]]])

    A2 = np.linalg.det(tmp)  # Double of triangle area
    A = A2 / 2.0
    return A


def tri3_Bmatrix(ex, ey):
    """
    Compute the strain displacement matrix for a 3 node triangular membrane element.

    :param list ex: element x coordinates [x1, x2, x3]
    :param list ey: element y coordinates [y1, y2, y3]
    :return  [3 x 6] strain displacement matrix
    """

    A = tri3_area(ex, ey)
    A2 = 2.0 * A

    cyclic_ijk = [0, 1, 2, 0, 1]  # Cyclic permutation of the nodes i,j,k

    zi_px = np.zeros(3)  # Partial derivative with respect to x
    zi_py = np.zeros(3)  # Partial derivative with respect to y

    for i in range(3):
        j = cyclic_ijk[i + 1]
        k = cyclic_ijk[i + 2]
        zi_px[i] = (ey[j] - ey[k]) / A2
        zi_py[i] = (ex[k] - ex[j]) / A2

    B = np.array([
        [zi_px[0], 0, zi_px[1], 0, zi_px[2], 0],
        [0, zi_py[0], 0, zi_py[1], 0, zi_py[2]],
        [zi_py[0], zi_px[0], zi_py[1], zi_px[1], zi_py[2], zi_px[2]]])

    return B


def tri3_Kmatrix(ex, ey, D, th, eq=None):
    """
    Compute the stiffness matrix for a two dimensional beam element.

    :param list ex: element x coordinates [x1, x2, x3]
    :param list ey: element y coordinates [y1, y2, y3]
    :param list D : 2D constitutive matrix
    :param list th: element thickness
    :param list eq: distributed loads, local directions [bx, by]
    :return mat Ke: element stiffness matrix [6 x 6]
    :return mat fe: consistent load vector [6 x 1] (if eq!=None)
    """

    A = tri3_area(ex, ey)
    A2 = 2.0 * A

    cyclic_ijk = [0, 1, 2, 0, 1]  # Cyclic permutation of the nodes i,j,k

    zi_px = np.zeros(3)  # Partial derivative with respect to x
    zi_py = np.zeros(3)  # Partial derivative with respect to y

    for i in range(3):
        j = cyclic_ijk[i + 1]
        k = cyclic_ijk[i + 2]
        zi_px[i] = (ey[j] - ey[k]) / A2
        zi_py[i] = (ex[k] - ex[j]) / A2

    B = tri3_Bmatrix(ex, ey)

    Ke = (B.T @ D @ B) * A * th

    if eq is None:
        fe = np.zeros((6,1))
    else:
        fx = A * th * eq[0] / 3.0
        fy = A * th * eq[1] / 3.0
        fe = np.array([[fx], [fy], [fx], [fy], [fx], [fy]])

    return Ke, fe


def tri3_cornerstresses(ex, ey, D, th, elDispVec):
    """
    Compute the corner stresses for all 3 corner nodes

    :param list ex: element x coordinates [x1, x2, x3]
    :param list ey: element y coordinates [y1, y2, y3]
    :param list D : 2D constitutive matrix
    :param list th: element thickness
    :return list of list of corner stresses
    """

    B = tri3_Bmatrix(ex, ey)

    strain = B @ elDispVec
    stress = D @ strain

    cornerStresses = []
    for inod in range(3):
        cornerStresses.append([stress[0], stress[1], stress[2]])

    return cornerStresses
    
def zeta_partials_x_and_y(ex,ey):
    """
    Compute partials of area coordinates with respect to x and y.
    
    :param list ex: element x coordinates [x1, x2, x3]
    :param list ey: element y coordinates [y1, y2, y3]
    """
    
    A = tri3_area(ex, ey)
    A2 = 2.0 * A
       
    cyclic_ijk = [0,1,2,0,1]      # Cyclic permutation of the nodes i,j,k
    
    zeta_px = np.zeros(3)           # Partial derivative with respect to x
    zeta_py = np.zeros(3)           # Partial derivative with respect to y

    for i in range(3):
        j = cyclic_ijk[i+1]
        k = cyclic_ijk[i+2]
        zeta_px[i] = (ey[j] - ey[k]) / A2
        zeta_py[i] = (ex[k] - ex[j]) / A2

    return zeta_px, zeta_py

# Functions for 6 node triangle
    
def tri6_area(ex,ey):
        
    tmp = np.array([[1,ex[0],ey[0]],
                    [1,ex[1],ey[1]],
                    [1,ex[2],ey[2]]])
    
    A = np.linalg.det(tmp) / 2
    
    return A


def tri6_shape_functions(zeta):
    L1, L2, L3 = zeta
    return np.array([
        L1*(2*L1 - 1),       # N1 (hjørne 1)
        L2*(2*L2 - 1),       # N2 (hjørne 2)
        L3*(2*L3 - 1),       # N3 (hjørne 3)
        4*L1*L2,             # N4 (kant 1–2)
        4*L2*L3,             # N5 (kant 2–3)
        4*L3*L1              # N6 (kant 3–1)
    ])


def tri6_shape_function_partials_x_and_y(zeta,ex,ey):
    L1, L2, L3 = zeta
    zeta_px, zeta_py = zeta_partials_x_and_y(ex, ey)  # dL1/dx, dL2/dx, dL3/dx osv.

    dL1dx, dL2dx, dL3dx = zeta_px
    dL1dy, dL2dy, dL3dy = zeta_py

    N6_px = np.zeros(6)
    N6_py = np.zeros(6)

    # Hjørner
    N6_px[0] = (4*L1 - 1) * dL1dx
    N6_py[0] = (4*L1 - 1) * dL1dy

    N6_px[1] = (4*L2 - 1) * dL2dx
    N6_py[1] = (4*L2 - 1) * dL2dy

    N6_px[2] = (4*L3 - 1) * dL3dx
    N6_py[2] = (4*L3 - 1) * dL3dy

    # Kantmidter
    N6_px[3] = 4 * (L1 * dL2dx + L2 * dL1dx)
    N6_py[3] = 4 * (L1 * dL2dy + L2 * dL1dy)

    N6_px[4] = 4 * (L2 * dL3dx + L3 * dL2dx)
    N6_py[4] = 4 * (L2 * dL3dy + L3 * dL2dy)

    N6_px[5] = 4 * (L3 * dL1dx + L1 * dL3dx)
    N6_py[5] = 4 * (L3 * dL1dy + L1 * dL3dy)

    return N6_px, N6_py


def tri6_Bmatrix(zeta,ex,ey):
    
    nx,ny = tri6_shape_function_partials_x_and_y(zeta, ex, ey)

    Bmatrix = np.zeros((3,12))
    for i in range(6):
        c = 2*i
        Bmatrix[0,c  ] = nx[i]
        Bmatrix[1,c+1] = ny[i]
        Bmatrix[2,c  ] = ny[i]
        Bmatrix[2,c+1] = nx[i]


    return Bmatrix


def tri6_Kmatrix(ex,ey,D,th,eq=None):
    """
    Compute the stiffness matrix for a two dimensional beam element.

    :param list ex: element x coordinates [x1, x2, x3, x4, x5, x6]
    :param list ey: element y coordinates [y1, y2, y3, y4, y5, y6]
    :param list D : 2D constitutive matrix
    :param list th: element thickness
    :param list eq: distributed loads, local directions [bx, by]
    :return mat Ke: element stiffness matrix [12 x 12]
    :return mat fe: consistent load vector [12 x 1] (if eq!=None)
    """   

    zetaInt = np.array([[0.5,0.5,0.0],
                        [0.0,0.5,0.5],
                        [0.5,0.0,0.5]])
    
    wInt = np.array([1.0/3.0,1.0/3.0,1.0/3.0])

    A    = tri6_area(ex,ey)

    Ke = np.zeros((12,12))
    fe = np.zeros((12,1))

    body_force = None if eq is None else np.array(eq, dtype=float).reshape(2,1)

    for i_gp, zeta in enumerate(zetaInt):
        weight = wInt[i_gp] * A * th

        B = tri6_Bmatrix(zeta, ex, ey)
        Ke += B.T @ D @ B * weight

        if body_force is not None:
            Nvals = tri6_shape_functions(zeta)
            Nmat = np.zeros((2,12))
            for i_node, Ni in enumerate(Nvals):
                col = 2 * i_node
                Nmat[0, col] = Ni
                Nmat[1, col + 1] = Ni

            fe += Nmat.T @ body_force * weight

    return Ke, fe



def tri6_cornerstresses(ex, ey, D, th, elDispVec):
    """
    Compute the corner stresses for all 3 corner nodes

    :param list ex: element x coordinates [x1, x2, x3]
    :param list ey: element y coordinates [y1, y2, y3]
    :param list D : 2D constitutive matrix
    :param list th: element thickness
    :return list of list of corner stresses
    """
    
    # Corner values
    zetaCorner = np.array([[1.0, 0.0, 0.0],
                           [0.0, 1.0, 0.0],
                           [0.0, 0.0, 1.0]])

    cornerStresses = []
    for zeta in zetaCorner:
        B = tri6_Bmatrix(zeta, ex, ey)
        strain = B @ elDispVec
        stress = D @ strain
        cornerStresses.append(stress.flatten().tolist())

    return cornerStresses

  
