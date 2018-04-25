
function calor(x_end, t_end, N_x, N_t)
    dx = x_end/(N_x - 1)
    dt = t_end/(N_t - 1)
    mu = 1.0
    theta = 0.5
    thetaC = 1.0 - theta
    x = linspace(0.0, x_end, N_x)
    t = linspace(0.0, t_end, N_t)

    f(x, t) = 0.0

    #AA = zeros(N_x, N_x)
    #AA[1, 1] = 1.0
    #AA[end, end] = 1.0
    #for j in 2:N_x-1
    #    AA[j, j-1] = -thetaC*mu*dt/dx^2
    #    AA[j, j] = 1.0 + 2.0*thetaC*mu*dt/dx^2
    #    AA[j, j+1] = -thetaC*mu*dt/dx^2
    #end
    d1 = ones(N_x-1)*(-thetaC*mu*dt/dx^2)
    d2 = ones(N_x)*(1.0 + 2.0*thetaC*mu*dt/dx^2)
    d3 = ones(N_x-1)*(-thetaC*mu*dt/dx^2)
    d1[end] = d3[1] = 0.0
    d2[1] = d2[end] = 1.0

    AA = Tridiagonal(d1, d2, d3)
    AAA = factorize(AA)

    u = ones(N_x)
    f1 = theta*mu*dt/dx^2
    for i in 1:N_t-1
        D = zeros(N_x)
        for j in 2:N_x-1
            D[j] = f1*(u[j+1] - 2*u[j] + u[j-1]) + theta*dt*f(x[j], t[i]) + thetaC*dt*f(x[j], t[i+1]) + u[j]
        end
        D[1] = 0.0
        D[end] = 1.0

        u = AAA \ D
    end
    return x, u
end

tic()
x, u = calor(1.0, 0.1, 201, 2001)
toc()
#println("$u")
