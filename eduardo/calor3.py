# - - coding: utf8 --
import numpy
import scipy
import scipy.linalg
import matplotlib
import matplotlib.pyplot as plt

def TDMAsolver(a, b, c, d):
    '''
    TDMA solver, a b c d can be NumPy array type or Python list type.
    refer to http://en.wikipedia.org/wiki/Tridiagonal_matrix_algorithm
    and to http://www.cfd-online.com/Wiki/Tridiagonal_matrix_algorithm_-_TDMA_(Thomas_algorithm)
    '''
    nf = len(d) # number of equations
    ac, bc, cc, dc = map(numpy.array, (a, b, c, d)) # copy arrays
    for it in range(1, nf):
        mc = ac[it-1]/bc[it-1]
        bc[it] = bc[it] - mc*cc[it-1]
        dc[it] = dc[it] - mc*dc[it-1]

    xc = bc
    xc[-1] = dc[-1]/bc[-1]

    for il in range(nf-2, -1, -1):
        xc[il] = (dc[il]-cc[il]*xc[il+1])/bc[il]

    return xc

def calor(x_end, t_end, N_x, N_t):
    dx = x_end/(N_x - 1)
    dt = t_end/(N_t - 1)

    mu = 1.0

    theta = 0.5
    thetaC = 1.0 - theta

    x = numpy.linspace(0.0, x_end, N_x)
    t = numpy.linspace(0.0, t_end, N_t)

    def f(x, t):
        return numpy.zeros(x.shape)

        AA[j, j] = 1.0 + 2.0*thetaC*mu*dt/numpy.square(dx)
        AA[j, j+1] = -thetaC*mu*dt/numpy.square(dx)

    d1 = numpy.ones(N_x-1)*(-thetaC*mu*dt/numpy.square(dx))
    d2 = numpy.ones(N_x)*(1.0 + 2.0*thetaC*mu*dt/numpy.square(dx))
    d3 = numpy.ones(N_x-1)*(-thetaC*mu*dt/numpy.square(dx))
    d1[-1] = 0.0
    d3[0] = 0.0
    d2[0] = 1.0
    d2[-1] = 1.0

    u = numpy.ones(N_x)

    #fig, ax = plt.subplots()

    f1 = theta*mu*dt/numpy.square(dx)
    for i in range(1, N_t-1):
        D = numpy.zeros(N_x)
        #for j in range(1,N_x-1):
        #    D[j] = f1*(u[j+1] - 2*u[j] + u[j-1]) + theta*dt*f(x[j], t[i]) + thetaC*dt*f(x[j], t[i+1]) + u[j]
        D[1:-1] = f1*(u[2:] - 2.0*u[1:-1] + u[0:-2]) + u[1:-1] + theta*dt*f(x[1:-1], t[i]) + thetaC*dt*f(x[1:-1], t[i+1])
        D[0] = 0.0
        D[-1] = 1.0

        #u = numpy.linalg.solve(AA, D)
        #u = scipy.linalg.lu_solve(AAA, D)
        u = TDMAsolver(d1, d2, d3, D)

        #if i%50 == 0:
        #    ax.plot(x, u)

    #ax.plot(x, u)
    #plt.show()
    return [x, u]

#ax.set(xlabel='time (s)', ylabel='voltage (mV)',
       #title='About as simple as it gets, folks')
#ax.grid()

#fig.savefig("test.png")
#plt.show()
