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
        bTemp[k] = 1. - A - C[k] - m * C[k-1]
        rTemp[k] = normalVelocity[k] - m * rTemp[k-1]

    print('C')
    print(C)
    print('rTemp')
    print(rTemp)
    #print('bTemp')
    #print(bTemp)

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

    alpha = np.zeros((len(normalVelocity)))
    A = np.zeros((len(normalVelocity)))
    B = np.zeros((len(normalVelocity)))
    C = np.zeros((len(normalVelocity)))
    D = np.zeros((len(normalVelocity)))
    E = np.zeros((len(normalVelocity)))
    F = np.zeros((len(normalVelocity)))
    h = np.zeros((len(normalVelocity)))

    # Using POP notation now:

    # This time, we don't divide through by layerThickEdgeMean
    k = Nsurf
    dz = (layerThickEdgeMean[Nsurf] + layerThickEdgeMean[Nsurf+1]) / 2.
    h[Nsurf] = layerThickEdgeMean[Nsurf]
    A[Nsurf] = dt * vertViscTopOfEdge[Nsurf+1] / dz
    #D[Nsurf] = h[Nsurf] + A[Nsurf] # hfac_u(k) = dz(k)/c2dtu = dz/(2 * dt)
    alpha[Nsurf] = 0.0
    D[Nsurf] = h[Nsurf] * normalVelocity[Nsurf]
    E[Nsurf] = A[Nsurf] / (h[Nsurf] + A[Nsurf]) # alpha should really be at k-1 here
    # B[Nsurf] = h[Nsurf] * E[Nsurf]
    # B[Nsurf] = h[Nsurf] * A[Nsurf]
    #F[Nsurf] = h[Nsurf] * normalVelocity[Nsurf] / h[Nsurf] + D[Nsurf]
    F[Nsurf] = D[Nsurf] / (h[Nsurf] + A[Nsurf]) # alpha should really be at k-1 here

    # first pass: set the coefficients
    for k in range(Nsurf+1, N+1):
        h[k] = layerThickEdgeMean[k]
        C[k] = A[k-1]
        A[k] = ( 2. * dt * vertViscTopOfEdge[k] /
            (layerThickEdgeMean[k-1] + layerThickEdgeMean[k]))
        D[k] = h[k] * normalVelocity[k]
        #D[k] = h[k] + A[k] + B[k-1]
        alpha[k] = A[k]* (h[k]+alpha[k-1]) / (h[k] + A[k] + alpha[k-1])
        # E[k] = A[k] / D[k]
        E[k] = A[k] / (h[k] + A[k] + alpha[k-1])
        # B[k] = E[k] * (h[k] + B[k-1])
        #F[k] = ( h[k] * normalVelocity[k] + 
        #         C[k] * F[k-1] ) / D[k]
        F[k] = ( D[k] + A[k-1]*F[k-1] ) / (h[k] + A[k] + alpha[k-1])

    print('alpha')
    print(alpha)
    print('A')
    print(A)
    print('B')
    print(B)
    print('D')
    print(D)
    print('E')
    print(E)
    print('F')
    print(F)

    # second pass: back substitution

    # Treat the bottom index separately
    # We do not apply bottom drag, unlike mpas
    Ftemp = ( A[N]*F[N] ) / (alpha[N])
    # this is F[k+1] but without terms below boundary 
    normalVelocity[N] = F[N] + E[k] * Ftemp

    for k in range(N-1, Nsurf-1, -1):
        normalVelocity[k] = F[k] + E[k] * F[k+1]

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
    u_mod = _schopf_tridiag(nu, h, dt, u)
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
