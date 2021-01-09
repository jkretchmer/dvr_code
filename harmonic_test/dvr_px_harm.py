#Simple DVR code to caclulate the P(x) for a 1d harmonic potential
#Can compare to NVT PIMD simulation
#DVR hamiltonian follows D. T. Colbert and W. H. Miller JCP 1992

#I calculate full DVR hamiltonian as an example, but to calculate P(x) only need the diagonal terms
#i.e. P(x_i) = e^{-\beta H_ii} / Q, where Q = \sum_i e^{-\beta H_ii} and H is the DVR Hamiltonian matrix

import numpy as np
import utils


k  = 2.0
R0 = 1.0
m  = 1.0
beta = 1.0

N    = 101 #number of grid points, have the extra 1 to include the min of the potential as a grid point
xmin = -2.0 #min value of grid along x
xmax = 4.0 #max value of grid along x
delx = (xmax - xmin) / (N-1) #grid spacing


#Calculate DVR Hamiltonian
Hdvr = np.zeros([N,N])
prob_x = np.zeros([N,2])
for i in range(N):
    xi = xmin + i*delx
    prob_x[i,0] = xi
    for j in range(i,N):
        xj = xmin + j*delx

        fctr = (-1)**(i-j) / ( 2.0*m*delx**2 )

        if( i==j ): 
            Hdvr[i,j] = 0.5 * k * (xi-R0)**2 #potential term
            Hdvr[i,j] += fctr * np.pi**2/3.0 #kinetic term
        else:
            Hdvr[i,j] = fctr * 2.0 / (i-j)**2 #kinetic term

Hdvr = Hdvr + np.transpose( np.triu( Hdvr, 1 ) )

#Calculate probability distribution along x
prob_x[:,1] = np.exp( -beta * np.diag(Hdvr) )
Q = np.sum( prob_x[:,1] ) * delx #need delx here due to numerical integration of partition function to properly normalize
prob_x[:,1] = prob_x[:,1] / Q

utils.printarray( prob_x, 'prob_x.dat', True )
