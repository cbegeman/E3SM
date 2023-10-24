import numpy as np
import matplotlib.pyplot as plt

def _mpas_tridiag(vertViscTopOfEdge, layerThickEdgeMean, dt, normalVelocity):

    # tridiagonal matrix algorithm currently in mpas-ocean
    # vertViscTopOfEdge
    # layerThickEdgeMean
    # dt
    # normalVelocity of size layerThickEdgeMean
    # crop each of the vectors so they range from minLevelEdgeBot to maxLevelEdgeTop
    print('enter tridiag')
    Nsurf = 0  # top index
    N = len(normalVelocity) - 1  # bottom index

    C = np.zeros((len(normalVelocity)))
    bTemp = np.zeros((len(normalVelocity)))
    rTemp = np.zeros((len(normalVelocity)))

    C[Nsurf] = ( -2. * dt * vertViscTopOfEdge[Nsurf+1] /
               (layerThickEdgeMean[Nsurf] + layerThickEdgeMean[Nsurf+1]) /
               layerThickEdgeMean[Nsurf] )
    bTemp[Nsurf] = 1. - C[Nsurf]
    rTemp[Nsurf] = normalVelocity[Nsurf]

    # first pass: set the coefficients
    for k in range(Nsurf+1, N):
        A = ( -2. * dt * vertViscTopOfEdge[k] /
            (layerThickEdgeMean[k-1] + layerThickEdgeMean[k]) /
             layerThickEdgeMean[k] )
        m = A / bTemp[k-1]
        C[k] = ( -2. * dt * vertViscTopOfEdge[k+1] /
               (layerThickEdgeMean[k] + layerThickEdgeMean[k+1]) /
               layerThickEdgeMean[k] )
        bTemp[k] = 1. - A - C[k] - m * C[k-1] # This is the same as the denominator in E14 of Schopf / h[k]
        # bTemp[k] = 1. - A - A[k-1](1 - C[k-1]/bTemp[k-1]) # this isn't entirely correct because m = f(A[k])
        # E[k] = C[k] / bTemp[k]
        # bTemp[k] = 1. - A - A[k-1](1 - E[k-1])
        # bTemp[k] = 1. - A - alpha[k-1]
        # E[k] = C[k] / (1. - A - C[k] - m * C[k-1]) # this is the version currently in MPAS
        rTemp[k] = normalVelocity[k] - m * rTemp[k-1]

    print('C')
    print(C)
    print('rTemp')
    print(rTemp)
    print('bTemp')
    print(bTemp)

    A = ( -2. * dt * vertViscTopOfEdge[N] /
        (layerThickEdgeMean[N-1] + layerThickEdgeMean[N]) /
        layerThickEdgeMean[N] )
    m = A / bTemp[N-1]

    # We do not apply bottom drag, unlike mpas
    normalVelocity[N] = ( (normalVelocity[N] - m * rTemp[N-1]) /
                          (1. - A - m * C[N-1]) )
    # second pass: back substitution
    for k in range(N-1, Nsurf-1, -1):
        normalVelocity[k] = (rTemp[k] - C[k] * normalVelocity[k+1]) / bTemp[k]

    return normalVelocity

def _mod_tridiag(vertViscTopOfEdge, layerThickEdgeMean, dt, normalVelocity):

    # tridiagonal matrix algorithm currently in mpas-ocean
    # vertViscTopOfEdge
    # layerThickEdgeMean
    # dt
    # normalVelocity of size layerThickEdgeMean
    # crop each of the vectors so they range from minLevelEdgeBot to maxLevelEdgeTop
    print('enter tridiag')
    Nsurf = 0  # top index
    N = len(normalVelocity) - 1  # bottom index

    C = np.zeros((len(normalVelocity)))
    bTemp = np.zeros((len(normalVelocity)))
    rTemp = np.zeros((len(normalVelocity)))
    alpha = np.zeros((len(normalVelocity)))

    C[Nsurf] = ( 2. * dt * vertViscTopOfEdge[Nsurf+1] /
               (layerThickEdgeMean[Nsurf] + layerThickEdgeMean[Nsurf+1]) /
               layerThickEdgeMean[Nsurf] )
    bTemp[Nsurf] = 1. + C[Nsurf] # - alpha[k-1]
    rTemp[Nsurf] = normalVelocity[Nsurf]
    alpha[Nsurf] = -C[Nsurf]*(1 + C[Nsurf]/bTemp[Nsurf])
    # first pass: set the coefficients
    for k in range(Nsurf+1, N):
        C[k] = ( 2. * dt * vertViscTopOfEdge[k+1] /
               (layerThickEdgeMean[k] + layerThickEdgeMean[k+1]) /
               layerThickEdgeMean[k] )
        bTemp[k] = 1. + C[k] - alpha[k-1]
        alpha[k] = -C[k]*(1 + C[k]/bTemp[k])
        rTemp[k] = normalVelocity[k] + C[k] * rTemp[k-1]/bTemp[k-1]

    print('C')
    print(C)
    print('rTemp')
    print(rTemp)
    print('alpha')
    print(alpha)
    print('btemp')
    print(bTemp)

    # We do not apply bottom drag, unlike mpas
    bTemp[N] = (1. - alpha[N-1]) # + C[k], C[k]=0
    normalVelocity[N] = ( (normalVelocity[N] + (C[N-1] / bTemp[N-1]) * rTemp[N-1]) /
                          bTemp[N] )
    # second pass: back substitution
    for k in range(N-1, Nsurf-1, -1):
        normalVelocity[k] = (rTemp[k] + C[k] * normalVelocity[k+1]) / bTemp[k]

    return normalVelocity

def _schopf_tridiag(vertViscTopOfEdge, layerThickEdgeMean, dt, normalVelocity):

    # tridiagonal matrix algorithm currently in mpas-ocean
    # vertViscTopOfEdge
    # layerThickEdgeMean
    # dt
    # normalVelocity of size layerThickEdgeMean
    # crop each of the vectors so they range from minLevelEdgeBot to maxLevelEdgeTop
    print('enter tridiag')
    Nsurf = 0  # top index
    N = len(normalVelocity) - 1  # bottom index
    h = layerThickEdgeMean

    alpha = np.zeros((len(normalVelocity)))
    # 0   : B A      R
    # k-1 : C B A    R
    # k   :   C B A  R
    # N   :     C B  R
    A = np.zeros((len(normalVelocity)))
    B = np.zeros((len(normalVelocity)))
    C = np.zeros((len(normalVelocity)))
    D = np.zeros((len(normalVelocity)))
    DD = np.zeros((len(normalVelocity)))
    E = np.zeros((len(normalVelocity)))
    EE = np.zeros((len(normalVelocity)))
    F = np.zeros((len(normalVelocity)))
    FF = np.zeros((len(normalVelocity)))
    U = np.zeros((len(normalVelocity)))

    # Using POP notation now:

    # This time, we don't divide through by layerThickEdgeMean
    k = Nsurf
    dz = (layerThickEdgeMean[k] + layerThickEdgeMean[k+1]) / 2.
    A[k] = dt * vertViscTopOfEdge[k+1] / dz
    alpha[k] = 0.0
    D[k] = h[k] * normalVelocity[k]
    DD[k] = h[k] + A[k] # POP: hfac_u(k) = dz(k)/c2dtu = dz/(2 * dt)
    E[k] = A[k] / DD[k] # alpha should really be at k-1 here
    EE[k] = A[k] / DD[k]
    B[Nsurf] = h[Nsurf] * E[Nsurf]
    FF[k] = h[k] * normalVelocity[k] / DD[k]
    F[k] = D[Nsurf] / DD[k] # alpha should really be at k-1 here

    # first pass: set the coefficients
    for k in range(Nsurf+1, N+1):
        C[k] = A[k-1]
        dz = (layerThickEdgeMean[k-1] + layerThickEdgeMean[k]) / 2.
        A[k] = 2. * dt * vertViscTopOfEdge[k] / dz
        D[k] = h[k] * normalVelocity[k]
        DD[k] = h[k] + A[k] + B[k-1]
        alpha[k] = A[k]* (h[k]+alpha[k-1]) / (h[k] + A[k] + alpha[k-1])
        EE[k] = A[k] / DD[k]
        E[k] = A[k] / (h[k] + A[k] + alpha[k-1])
        B[k] = A[k] * (h[k] + B[k-1]) / DD[k]
        FF[k] = ( h[k] * normalVelocity[k] + 
                 A[k-1] * FF[k-1] ) / DD[k]
        F[k] = ( D[k] + A[k-1]*F[k-1] ) / (h[k] + A[k] + alpha[k-1])

    #print('alpha')
    #print(alpha)
    #print('A')
    #print(A)
    #print('B')
    #print(B)
    #print('D')
    #print(D)
    print('EE')
    print(EE)
    print('E')
    print(E)
    print('FF')
    print(FF)
    print('F')
    print(F)

    # second pass: back substitution

    # Treat the bottom index separately
    # We do not apply bottom drag, unlike mpas
    normalVelocity[N] = F[N]
    U[N] = FF[N]

    for k in range(N-1, Nsurf-1, -1):
        normalVelocity[k] = F[k] + E[k] * F[k+1]
        U[k] = FF[k] + EE[k] * FF[k+1]
    print('U')
    print(U)
    return normalVelocity

def main():
    dt = 1.
    nz = 20
    h = 2. * np.ones((nz))
    nu = 2. * np.ones((nz))
    z = np.linspace(0, nz, num=nz)
    u = np.linspace(1, 0, num=nz)
    u_old = u.copy()
    
    print('MPAS')
    print('Input:')
    print(u_old)
    u_new = _mpas_tridiag(nu, h, dt, u_old)
    u_mpas = u_new.copy()
    print('Output:')
    print(u_mpas)
    print('Schopf')
    print('Input:')
    print(u)
    #u_mod = _schopf_tridiag(nu, h, dt, u)
    u_mod = _mod_tridiag(nu, h, dt, u)
    print('Output:')
    print(u_mod)
    fig = plt.figure()
    plt.plot(u_new - u_old, z, 'k:', label='old')
    plt.plot(u_mod - u_old, z, 'k:', label='mod')
    #plt.plot(u_old, z, 'k:', label='u_old')
    #plt.plot(u_new, z, 'k-', label='u_new')
    plt.legend()
    plt.savefig('tridiag.png')

if __name__ == "__main__":
    main()
