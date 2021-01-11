#Simple DVR code to caclulate the P(x) for a 1d harmonic potential
#Can compare to NVT PIMD simulation
#DVR hamiltonian follows D. T. Colbert and W. H. Miller JCP 1992

#I calculate full DVR hamiltonian as an example, but to calculate P(x) only need the diagonal terms
#i.e. P(x_i) = e^{-\beta H_ii} / Q, where Q = \sum_i e^{-\beta H_ii} and H is the DVR Hamiltonian matrix

import numpy as np
import utils
from scipy.linalg import expm

def get_Vel(x):

    k1 = 4e-5
    k2 = 3.2e-5
    R1 = -1.75
    R2 = 1.75
    eps1 = 0.0
    eps2 = 2.28e-5

    c = 5e-5
    alpha = 0.4
    R12 = 0.0

    Vel = np.zeros([2,2])
    Vel[0,0] = 0.5 * k1 * (x-R1)**2 + eps1
    Vel[1,1] = 0.5 * k2 * (x-R2)**2 + eps2
    Vel[0,1] = c * np.exp( -alpha * (x-R12)**2 )
    Vel[1,0] = Vel[0,1]

    return Vel


Nstates = 2
m       = 3600.0
temp    = 8.0 / 3.1577455e5 
beta    = 1.0 / temp

N    = 51 #number of grid points, have the extra 1 to include the min of the potential as a grid point
xmin = -6.0
xmax = 6.0
delx = (xmax - xmin) / (N-1) #grid spacing


#Calculate DVR Hamiltonian
Hdvr = np.zeros([2*N,2*N])
xpos = np.zeros(N)
for state1 in range(Nstates):
    for i in range(N):

        indx1 = i + state1*N

        xi = xmin + i*delx
        xpos[i] = xi

        Vel = get_Vel(xi)

        for state2 in range(state1,Nstates):
            for j in range(i,N):

                indx2 = j + state2*N

                xj = xmin + j*delx
    
                fctr = (-1)**(i-j) / ( 2.0*m*delx**2 )

                if( indx1 == indx2 ): #Diagonal terms of dvr matrix
                    Hdvr[indx1,indx2] = fctr * np.pi**2/3.0 #kinetic term
                    Hdvr[indx1,indx2] = Vel[ state1, state1 ] #state-dependent potential term
                else:
                    if( state1 == state2 ): #Kinetic contribution from nuclei, so off-diagonal in nuclei, but still diagonal for electronic state
                        Hdvr[indx1,indx2] = fctr * 2.0 / (i-j)**2
                    elif( i == j ): #Off diagonal coupling of states, so off diagonal in state, but still diagonal for nuclei
                        Hdvr[indx1,indx2] = Vel[ state1, state2 ]

Hdvr = Hdvr + np.transpose( np.triu( Hdvr, 1 ) )


densmat = expm( -beta*Hdvr )
Q = np.trace( densmat ) * delx
densmat = densmat / Q

#State-dependent probabilities
px_1 = np.zeros([N,2])
px_1[:,0] = np.copy(xpos)
px_1[:,1] = np.diag(densmat)[:N]

px_2 = np.zeros([N,2])
px_2[:,0] = np.copy(xpos)
px_2[:,1] = np.diag(densmat)[N:]

#Joint probability
px = np.zeros([N,2])
px[:,0] = np.copy(xpos)
px[:,1] = px_1[:,1] + px_2[:,1]

utils.printarray( px_1, 'prob_x_1.dat', True )
utils.printarray( px_2, 'prob_x_2.dat', True )
utils.printarray( px, 'prob_x.dat', True )

