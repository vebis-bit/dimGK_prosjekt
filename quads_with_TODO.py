# -*- coding: utf-8 -*-
"""
Created on Sun Oct 21 08:15:51 2018

@author: bjohau
"""
import numpy as np
import sys

def gauss_points(iRule):
    """
    Returns gauss coordinates and weight given integration number

    Parameters:

        iRule = number of integration points

    Returns:

        gp : row-vector containing gauss coordinates
        gw : row-vector containing gauss weight for integration point

    """
    gauss_position = [[ 0.000000000],
                      [-0.577350269,  0.577350269],
                      [-0.774596669,  0.000000000,  0.774596669],
                      [-0.8611363116, -0.3399810436, 0.3399810436, 0.8611363116],
                      [-0.9061798459, -0.5384693101, 0.0000000000, 0.5384693101, 0.9061798459]]
    gauss_weight   = [[2.000000000],
                      [1.000000000,   1.000000000],
                      [0.555555556,   0.888888889,  0.555555556],
                      [0.3478548451,  0.6521451549, 0.6521451549, 0.3478548451],
                      [0.2369268850,  0.4786286705, 0.5688888889, 0.4786286705, 0.2369268850]]


    if iRule < 1 and iRule > 5:
        sys.exit("Invalid number of integration points.")

    idx = iRule - 1
    return gauss_position[idx], gauss_weight[idx]


def quad4_shapefuncs(xsi, eta):
    """
    Calculates shape functions evaluated at xi, eta
    """
    # ----- Shape functions -----
    N = np.zeros(4)
    for i in range(4):
        xsi_factor = 1
        eta_factor = 1
        if i == 1 or i == 2:
            eta_factor = -1
        if i == 2 or i == 3:
            xsi_factor = -1
        N[i] = 1/4*(1+xsi*xsi_factor)*(1+eta*eta_factor)
    return N

def quad4_shapefuncs_grad_xsi(xsi, eta):
    """
    Calculates derivatives of shape functions wrt. xsi
    """
    # ----- Derivatives of shape functions with respect to xsi -----
    Ndxi = np.zeros(4)
    for i in range(4):
        xsi_factor = 1
        factor = 1
        if i == 0 or i == 1:
            xsi_factor = -1
        if i == 0 or i == 3:
            factor = -1
        Ndxi[i] = factor*1/4*(1 + xsi*xsi_factor)

    return Ndxi


def quad4_shapefuncs_grad_eta(xsi, eta):
    """
    Calculates derivatives of shape functions wrt. eta
    """
    # ----- Derivatives of shape functions with respect to eta -----
    
    Ndeta = np.zeros(4)
    for i in range(4):
        eta_factor = 1
        factor = 1
        if i == 0 or i == 3:
            eta_factor = -1
        if i == 0 or i == 1:
            factor = -1
        Ndeta[i] = factor*1/4*(1 + eta*eta_factor)
    return Ndeta




def quad4_Kmatrix(ex, ey, D, thickness, eq=None):
    """
    Compute the stiffness matrix for a four node membrane element.

    Parameters:

        ex  = [x1 ... x4]           Element coordinates. Row matrix
        ey  = [y1 ... y4]
        D   =           Constitutive matrix
        thickness:      Element thickness
        eq = [bx; by]       bx:     body force in x direction
                            by:     body force in y direction

    Returns:

        Ke : element stiffness matrix (8 x 8)
        fe : equivalent nodal forces (4 x 1)

    """
    t = thickness

    if eq is None:
        f = np.zeros((2,1))  # Create zero matrix for load if load is zero
    else:
        f = np.array([eq]).T  # Convert load to 2x1 matrix

    Ke = np.zeros((8,8))        # Create zero matrix for stiffness matrix
    fe = np.zeros((8,1))        # Create zero matrix for distributed load

    numGaussPoints = 2  # Number of integration points
    gp, gw = gauss_points(numGaussPoints)  # Get integration points and -weight

    for iGauss in range(numGaussPoints):  # Solves for K and fe at all integration points
        for jGauss in range(numGaussPoints):

            xsi = gp[iGauss]
            eta = gp[jGauss]

            Ndxsi = quad4_shapefuncs_grad_xsi(xsi, eta)
            Ndeta = quad4_shapefuncs_grad_eta(xsi, eta)
            N1    = quad4_shapefuncs(xsi, eta)  # Collect shape functions evaluated at xi and eta

            # Matrix H and G defined according to page 52 of Waløens notes
            H = np.array([ex, ey])    # Collect global x- and y coordinates in one matrix
            G = np.array([Ndxsi, Ndeta])  # Collect gradients of shape function evaluated at xi and eta

            
            J = H@G.T          # Jacobian matrix
            invJ = np.linalg.inv(J)  # Inverse of Jacobian
            detJ = np.linalg.det(J)  # Determinant of Jacobian

            dN = invJ @ G  # Derivatives of shape functions with respect to x and y
            dNdx = dN[0]
            dNdy = dN[1]

            # Strain displacement matrix calculated at position xsi, eta

            
            B  = np.zeros((3,8))
            for i in range(0,8,2):
                B[0,i]   = dNdx[i//2]
                B[1,i+1] = dNdy[i//2]
                B[2,i]   = dNdy[i//2]
                B[2,i+1] = dNdx[i//2]

            
            N2 = np.zeros((2,8))
            for i in range(0,8,2):
                N2[0,i]   = N1[i//2]
                N2[1,i+1] = N1[i//2]

            # Evaluates integrand at current integration points and adds to final solution
            Ke += (B.T) @ D @ B * detJ * t * gw[iGauss] * gw[jGauss]
            fe += (N2.T) @ f    * detJ * t * gw[iGauss] * gw[jGauss]

    
    #Ke = np.eye(8) * 1.0e6

    return Ke, fe  # Returns stiffness matrix and nodal force vector


def quad4_cornerstresses(ex, ey, D, th, elDispVec):
    """
    Compute the corner stresses for all 3 corner nodes

    :param list ex: element x coordinates [x1, x2, x3, x4]
    :param list ey: element y coordinates [y1, y2, y3, y4]
    :param list D : 2D constitutive matrix
    :param list th: element thickness
    :return list of list of corner stresses
    """

    # Corner natural coordinates
    xi_eta_corner = np.array([[1.0, 1.0],
                           [-1.0, 1.0],
                           [-1.0, -1.0],
                           [1.0, -1.0]])
    
    cornerStresses = []

    for xsi, eta in xi_eta_corner:

        Ndxsi = quad4_shapefuncs_grad_xsi(xsi, eta)
        Ndeta = quad4_shapefuncs_grad_eta(xsi, eta)
        H = np.array([ex, ey])    # Collect global x- and y coordinates in one matrix
        G = np.array([Ndxsi, Ndeta])  # Collect gradients of shape function evaluated at xi and eta
        J = H@G.T
        invJ = np.linalg.inv(J)  # Inverse of Jacobian
        dN = invJ @ G  # Derivatives of shape functions with respect to x and y
        dNdx, dNdy = dN[0], dN[1]
        B  = np.zeros((3,8))
        for i in range(0,8,2):
            B[0,i]   = dNdx[i//2]
            B[1,i+1] = dNdy[i//2]
            B[2,i]   = dNdy[i//2]
            B[2,i+1] = dNdx[i//2]
        eps = B @ elDispVec  # Strain at corner points
        sigma = D @ eps      # Stress at corner points
        cornerStresses.append(sigma.flatten().tolist())

    return cornerStresses

def shape_funcs_and_grads(xsi, eta):
        """Quadratic 9-node Lagrange shape functions and natural gradients."""
        N1 = 0.25*xsi*eta*(xsi-1)*(eta-1)
        N2 = 0.25*xsi*eta*(xsi+1)*(eta-1)
        N3 = 0.25*xsi*eta*(xsi+1)*(eta+1)
        N4 = 0.25*xsi*eta*(xsi-1)*(eta+1)
        N5 = 0.5*(1 - xsi**2)*eta*(eta-1)
        N6 = 0.5*xsi*(xsi+1)*(1 - eta**2)
        N7 = 0.5*(1 - xsi**2)*eta*(eta+1)
        N8 = 0.5*xsi*(xsi-1)*(1 - eta**2)
        N9 = (1 - xsi**2)*(1 - eta**2)

        dNdxi = np.array([
            0.25*eta*(2*xsi-1)*(eta-1),
            0.25*eta*(2*xsi+1)*(eta-1),
            0.25*eta*(2*xsi+1)*(eta+1),
            0.25*eta*(2*xsi-1)*(eta+1),
           -xsi*eta*(eta-1),
            0.5*(2*xsi+1)*(1 - eta**2),
           -xsi*eta*(eta+1),
            0.5*(2*xsi-1)*(1 - eta**2),
           -2*xsi*(1 - eta**2)
        ])

        dNdeta = np.array([
            0.25*xsi*(xsi-1)*(2*eta-1),
            0.25*xsi*(xsi+1)*(2*eta-1),
            0.25*xsi*(xsi+1)*(2*eta+1),
            0.25*xsi*(xsi-1)*(2*eta+1),
            0.5*(1 - xsi**2)*(2*eta-1),
           -eta*xsi*(xsi+1),
            0.5*(1 - xsi**2)*(2*eta+1),
           -eta*xsi*(xsi-1),
           -2*eta*(1 - xsi**2)
        ])

        N = np.array([N1, N2, N3, N4, N5, N6, N7, N8, N9])
        return N, dNdxi, dNdeta

def quad9_Kmatrix(ex,ey,D,th,eq=None):
    """
    Calculates the stiffness matrix for a 9 node isoparametric element in plane stress

    :param list ex: element x coordinates [x1, x2, x3]
    :param list ey: element y coordinates [y1, y2, y3]
    :param list D : 2D constitutive matrix
    :param list th: element thickness
    :param list eq: distributed loads, local directions [bx, by]
    :return mat Ke: element stiffness matrix [6 x 6]
    :return mat fe: consistent load vector [6 x 1] (if eq!=None)
    """

    t = th
    f = np.zeros((2,1)) if eq is None else np.array(eq, dtype=float).reshape(2,1)

    Ke = np.zeros((18,18))
    fe = np.zeros((18,1))

    

    gp, gw = gauss_points(3)

    for iGauss, xsi in enumerate(gp):
        for jGauss, eta in enumerate(gp):
            wx = gw[iGauss]
            wy = gw[jGauss]

            N, dNdxi, dNdeta = shape_funcs_and_grads(xsi, eta)

            H = np.array([ex, ey])        # 2 x 9
            G = np.array([dNdxi, dNdeta]) # 2 x 9
            J = H @ G.T                   # 2 x 2
            invJ = np.linalg.inv(J)
            detJ = np.linalg.det(J)

            dN = invJ @ G                 # 2 x 9 -> global gradients
            dNdx, dNdy = dN[0], dN[1]

            B = np.zeros((3,18))
            for a in range(9):
                col = 2*a
                B[0, col]   = dNdx[a]
                B[1, col+1] = dNdy[a]
                B[2, col]   = dNdy[a]
                B[2, col+1] = dNdx[a]

            N2 = np.zeros((2,18))
            for a in range(9):
                col = 2*a
                N2[0, col]   = N[a]
                N2[1, col+1] = N[a]

            weight = detJ * t * wx * wy
            Ke += B.T @ D @ B * weight
            fe += N2.T @ f * weight

    return Ke, fe



def quad9_cornerstresses(ex, ey, D, th, elDispVec):
    """
    Compute the corner stresses for all 3 corner nodes

    :param list ex: element x coordinates [x1, x2, x3, x4, x5, x6, x7, x8, x9]
    :param list ey: element y coordinates [y1, y2, y3, y4, y5, y6, y7, y8, y9]
    :param list D : 2D constitutive matrix
    :param list th: element thickness
    :return list of list of corner stresses
    """

    # Corner natural coordinates
    xi_eta_corner = np.array([[1.0, 1.0],
                           [-1.0, 1.0],
                           [-1.0, -1.0],
                           [1.0, -1.0]])

    cornerStresses = []

    for xsi, eta in xi_eta_corner:
        N, dNdxi, dNdeta = shape_funcs_and_grads(xsi, eta)
        H = np.array([ex, ey])    # Collect global x- and y coordinates in one matrix
        G = np.array([dNdxi, dNdeta])  # Collect gradients of shape function evaluated at xi and eta
        J = H@G.T
        invJ = np.linalg.inv(J)  # Inverse of Jacobian
        dN = invJ @ G  # Derivatives of shape functions with respect to x and y
        dNdx, dNdy = dN[0], dN[1]
        B  = np.zeros((3,18))
        for i in range(9):
            c = 2*i
            B[0,c]   = dNdx[i]
            B[1,c+1] = dNdy[i]
            B[2,c]   = dNdy[i]
            B[2,c+1] = dNdx[i]
        eps = B @ elDispVec  # Strain at corner points
        sigma = D @ eps      # Stress at corner points
        cornerStresses.append(sigma.flatten().tolist())

    return cornerStresses
  
