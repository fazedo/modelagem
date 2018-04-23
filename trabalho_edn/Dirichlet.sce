//Primeiramente, implementamos os métodos de resolução de edo's.





function [y]=RK4(y0,x0,h,f,Nmax)
//O rk4 foi impletado pelo professor Guidi
//nas aulas de MMCCII. Fiz algumas alterações
//para que a codição inicial saia na matriz resultante. 

x_0=x0
y_0=y0
y=zeros(size(y0,1),Nmax)
y(:,1)=y_0

for i=2:Nmax
    k1=f(x_0,y_0)
    k2=f(x_0+0.5*h,y_0+0.5*h*k1)
    k3=f(x_0+0.5*h,y_0+0.5*h*k2)
    k4=f(x_0+h,y_0+h*k3)
    yi=y_0+h/6*(k1+2*(k2+k3)+k4)

    y(:,i)=yi
 
    x_0=x_0+h 
    y_0=yi
 
end

endfunction

function [y]=Adams_Bashfort(y0,x0,xf,f,Nmax)
//y0 é a matriz cujo as linhas correspondem a discretização no espaço
//e as colunas o ponto no tempo em que foi calculada, neste caso,
//as colunas vão ser os valores já inicializados pelo método RK.
//x0 e xf corresponderá ao tempo inicial e final.
//f é a função ut.
//Nmax é o número de pontos na malha do tempo.
    
x_0=x0
y_0=y0
K=size(y_0,2)
y=zeros(size(y_0,1),1)
h=(xf-x0)/Nmax
x=linspace(x0,xf,Nmax)
for i=1:K
    plot(y_0(:,i))
end 

//Calcula os coeficientes para k passos.
//A relação matricial obtive em "A Matrix System for Computing the Coefficients
//of the Adams Bashforth-Moulton Predictor-Corrector formulae" de Baba Seidu.
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


//Para calcular as diferenças finitas criaremos o vetor D nos quais a i-ésima
//entrada é Delta^ifn, no caso a i-ésima diferença finita de fn.
//Inicializando o vetor e o vetor auxiliar das diferenças finitas.
D=zeros(size(y_0,1),K); D0=zeros(size(y_0,1),K)
for i=1:K
    D(:,i) = f(x(K),y_0(:,K))
end
fi=1
for i=2:K
    //Para evitar o calculo dos fatorias que aparecem ao expandir os termos
    //das diferenças finitas dentro de cada laço, usamos as variáveis 
    //fi, fj, fji correspondendo, respectivamente, a (i-1)!, (j-1)! e
    //(j-i)!.
    fi=fi*(i-1); fj=fi; fji=1
    for j=i:K
        fj = fj*(j-1)
        D(:,j) = D(:,j) + (-1)^(i-1)*fj/(fji*fi)*f(x(K-i+1),y_0(:,K-i+1))
        fji = fji*(j+1-i)
    end 
end

//Por fim, o método!
for i=K+1:Nmax
    y = y_0(:,K)
    for j=1:K
        y = y + h*coef(j)*D(:,j) 
    end
    
    //atualiza y_0.
    for j=1:K-1
        y_0(:,j) = y_0(:,j+1)
    end
    y_0(:,K) = y
    
    //atualiza D.
    D0 = D
    D(:,1) = f(x(i),y)
    for j=2:K
        D(:,j) = D(:,j-1) - D0(:,j-1)
    end
    
    plot(y)
end
    
endfunction





//Implementamos a edo com condições de Dirichlet

//Definimos as demais constantes necessárias
eta = 0.5
gama = 1

tf = .1                                                       //Tempo final
t0 = 0                                                        //Tempo inicial
Nt = 801                                                     //Pontos no tempo
dt=(tf-t0)/Nt                                                 //Incremento no tempo

x0 = 0                                                        //Posição inicial
xf = 1                                                        //Posição final
Nx = 41                                                      //Pontos no espaço
dx=(xf-x0)/Nx                                                 //Incremento no espaço
x=linspace(x0,xf,Nx)                                          //Malha no espaço

c1 = (eta*dt)/(dx^2) - (gama*dt)/(2*dx)
c2 = 1 - (2*eta*dt)/(dx^2)
c3 = (eta*dt)/(dx^2) + (gama*dt)/(2*dx)

u=zeros(Nx,1)
for i=1:Nx
    u(i)=x(i)*(1-x(i))                                        //U(x,0) = x*(1-x)
end

function y = edo1(t,u)
    
    //Definimos as condições de Dirichlet.
    y(1) = 2*t                                                //U(0,t) = t^2
    y(Nx) = 3*t^2                                             //U(0,t) = t^3
    
    //Interior da malha.
    for i=2:Nx-1
        y(i) = c1*u(i-1) + c2*u(i) + c3*u(i+1) + x(i)+(t^2-0.05)       //f(x,t) = x+t
    end
    
endfunction

//Resolvendo a edo com condições de Dirichlet
[yi] = RK4(u,t0,dt,edo1,6)
[y] = Adams_Bashfort(yi,t0,tf,edo1,Nt)



