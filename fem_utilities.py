import numpy as np

def assem(edof, K, Ke, f=None, fe=None):
    """
    Assemble element matrices Ke ( and fe ) into the global
    stiffness matrix K ( and the global force vector f )
    according to the topology matrix edof.

    Parameters:
        edof        dof topology list (0-based, i.e. first dof is 0)
        K           the global stiffness matrix
        Ke          element stiffness matrix
        f           the global force vector
        fe          element force vector

    Output parameters:

        K           the updated global stiffness matrix
        f           the updated global force vector
    """

    idx = np.array(edof)
    K[np.ix_(idx, idx)] += Ke

    if (f is not None) and (fe is not None):
         f[np.ix_(idx)] += fe


def solveq(K, f, bcList, bcValList=None):
    """
    Solve static FE-equations considering boundary conditions.

    Parameters:
        K           global stiffness matrix, dim(K)= nd x nd
        f           global load vector, dim(f)= nd x 1

        bcPrescr    1-dim integer list containing prescribed dofs.
        bcVal       1-dim float list containing prescribed values.
                    If not given all prescribed dofs are assumed 0.
    Returns:
        a           solution including boundary values
        Q           reaction force vector
                    dim(a)=dim(Q)= nd x 1, nd : number of dof's
    """

    bcPrescr = np.array(bcList,dtype=int)  # Prescribed dofs (0-based)

    nDofs = K.shape[0]
    nPdofs = bcPrescr.shape[0]

    if bcValList is None:
        bcVal = np.zeros([nPdofs], dtype=float)
    else:
        bcVal    = np.array(bcValList)

    bcIsFree = np.ones(nDofs, dtype=bool) # Setting all dofs as free initially
    allDofs = np.arange(nDofs)

    bcIsFree[np.ix_(bcPrescr)] = False  # Since bcPrescr is 0-based
    freeDofs = allDofs[bcIsFree]

    fsys = f[freeDofs] - K[np.ix_((freeDofs), (bcPrescr - 1))] * np.asmatrix(bcVal).reshape(nPdofs, 1)
    asys = np.linalg.solve(K[np.ix_((freeDofs), (freeDofs))], fsys)

    a = np.zeros([nDofs, 1])
    a[np.ix_(bcPrescr)] = np.asmatrix(bcVal).reshape(nPdofs, 1)
    a[np.ix_(freeDofs)] = asys

    Q = K @ np.asmatrix(a) - f

    return (np.asmatrix(a), Q)


