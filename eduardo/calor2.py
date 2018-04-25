# - - coding: utf8 --
import numpy
import scipy
import scipy.linalg
import matplotlib
import matplotlib.pyplot as plt

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

    AA = numpy.zeros((N_x, N_x))
    AA[0, 0] = 1.0
    AA[-1, -1] = 1.0

    for j in range(1, N_x-1):
        AA[j, j-1] = -thetaC*mu*dt/numpy.square(dx)
        AA[j, j] = 1.0 + 2.0*thetaC*mu*dt/numpy.square(dx)
        AA[j, j+1] = -thetaC*mu*dt/numpy.square(dx)

    AAA = scipy.linalg.lu_factor(AA)

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

        u = numpy.linalg.solve(AA, D)
        u = scipy.linalg.lu_solve(AAA, D)

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
