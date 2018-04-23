//Implementamos a edp com condições de Newmann

function [y] = c0(t)
    y = t^2
endfunction

function [y] = cf(t)
    y = t^3
endfunction

function [y] = g(x,t)
    y = x + t 
endfunction

//Definimos as constantes necessárias
function [y] = edpNewmann(t0,tf,Nt,x0,xf,Nx,eta,gama,K,cNx0,cNxf,f)
    
dt=(tf-t0)/Nt                                                 //Incremento no tempo
dx=(xf-x0)/Nx                                                 //Incremento no espaço
x=linspace(x0,xf,Nx)                                          //Malha no espaço
t=linspace(t0,tf,Nt)                                          //Malha no tempo

c1 = (eta*dt)/(dx^2) - (gama*dt)/(2*dx)
c2 = 1 - (2*eta*dt)/(dx^2)
c3 = (eta*dt)/(dx^2) + (gama*dt)/(2*dx)

y_0=zeros(Nx,1)
for i=1:Nx
    y_0(i)=x(i)*(1-x(i))                                        //U(x,0) = x*(1-x)
end
plot(x',y_0)

//Inicializa o vetor das diferenças finitas
D=zeros(size(y_0,1)-2,K); D0=zeros(size(y_0,1)-2,K)

//Inicialização com RK4
y=zeros(size(y_0,1),K)
y(:,1)=y_0
for j=2:K
    
    for i=2:Nx-1
        k1= c1*y_0(i-1) + c2*y_0(i) + c3*y_0(i+1) + f(x(i),t(j))
        k2= c1*(y_0(i-1) + .5*dt*k1) + c2*(y_0(i) + + .5*dt*k1) + c3*(y_0(i+1) + + .5*dt*k1) + f(x(i),t(j) + .5*dt)
        k3= c1*(y_0(i-1) + .5*dt*k2) + c2*(y_0(i) + + .5*dt*k2) + c3*(y_0(i+1) + + .5*dt*k2) + f(x(i),t(j) + .5*dt)
        k4= c1*(y_0(i-1) + dt*k3) + c2*(y_0(i) + + dt*k3) + c3*(y_0(i+1) + + dt*k3) + f(x(i),t(j) + dt)
        y(i,j)=y_0(i)+dt/6*(k1+2*(k2+k3)+k4)
        
        D0(i-1,j-1) = k1
    end

    y(1,j) = y(2,j) - dx*cNx0(t(j))
    y(Nx,j) = y(Nx-1,j) + dx*cNxf(t(j))
    
    y_0=y(:,j)
    
    plot(x',y_0)
end

//Calcula o vetor das diferenças finitas
for j=1:K
    for i=2:Nx-1
        D(i-1,j) = c1*y(i-1,K) + c2*y(i,K) + c3*y(i+1,K) + f(x(i),t(K))
    end
end
fi=1
for i=2:K
    fi=fi*(i-1); fj=fi; fji=1
    for j=i:K
        fj = fj*(j-1)
        D(:,j) = D(:,j) + (-1)^(i-1)*fj/(fji*fi)*D0(:,K-i+1)
        fji = fji*(j+1-i)
    end 
end

//Calcula os coeficientes para k passos.
A=zeros(K,K); coef=zeros(K); b=zeros(K)
for i=1:K
    for j=1:K
        A(i,j)=(j-1)^(i-1)
    end
end
for i=1:K
    b(i)=(-1)^(i-1)*(1/i)
end
coef = A\b

//Por fim, o método!
for l=K+1:Nt
    y_0 = y(:,K)
    for j=1:K
        for i=2:Nx-1
            y_0(i) = y_0(i) + dt*coef(j)*D(i-1,j)
        end 
    end
    y_0(1) = y_0(2) - dt*cNx0(t(j))
    y_0(Nx) = y_0(Nx-1) + dt*cNxf(t(j))
    
    //atualiza y_0.
    for j=1:K-1
        y(:,j) = y(:,j+1)
    end
    y(:,K) = y_0
    
    //atualiza D.
    D0 = D
    for i=2:Nx-1
        D(i-1,1) = c1*y_0(i-1) + c2*y_0(i) + c3*y_0(i+1) + f(x(i),t(l))
    end
    for j=2:K
        D(:,j) = D(:,j-1) - D0(:,j-1)
    end
    
    plot(x',y_0)
end

endfunction
