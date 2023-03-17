import numpy as np

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
    print('bTemp')
    print(bTemp)
    print('rTemp')
    print(rTemp)

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

def main():
    print('Entering main')
    dt = 10.
    nz = 100
    h = np.ones((nz))
    nu = np.ones((nz))
    u_old = np.linspace(1, 0, num=nz)
    print(u_old)
    
    u_new = _mpas_tridiag(nu, h, dt, u_old)
    print(u_new)
    # return print(u_new)

if __name__ == "__main__":
    main()
