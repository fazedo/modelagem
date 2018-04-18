//Código para sistemas

function z=c_inicial_1(x)
    z=0.1//u
endfunction

function z=forcante_1(x,t)
    z=0
endfunction

function z=c_inicial_2(x)
    z=1//v
endfunction

function z=forcante_2(x,t)
    z=0
endfunction

function [y] = sistema(L,tf,N_x,N_t,theta)
    dt=tf/N_t
    dx=L/(N_x-1)

    x=[0:dx:L]   //Malha
 
    for k=1:N_x
        w(k)=c_inicial_1(x(k))
        w(N_x+k)=c_inicial_2(x(k))
    end
    
    alpha=zeros(2,2)
    alpha(1,1)= -1
    alpha(1,2)= .5    
    alpha(2,1)= .4
    alpha(2,2)= -1
    
    M= spzeros(2*N_x,2*N_x)  

    //diagonal principal de cada submatriz
    for k=2:N_x-1
        M(k,k)=-alpha(1,1)+2/dx^2 
        M(k,N_x+k)=-alpha(1,2)
    
        M(N_x+k,N_x+k)=-alpha(2,2)+2/dx^2
        M(N_x+k,k)=-alpha(2,1)
    end


    
        for k = 2:N_x-1 //preenche linhas de 2 a N_x-1
            M(k,k-1) = -1/dx^2    //diagonal inferior     
            M(k,k+1) = -1/dx^2    //diagonal superior     
            
            M(N_x+k,N_x+k-1) = -1/dx^2    //diagonal inferior     
            M(N_x+k,N_x+k+1) = -1/dx^2    //diagonal superior     
        end
        
        M=M*dt*(1-theta)

        for k=2:N_x-1
            M(k,k)=M(k,k)+1
            M(N_x+k,N_x+k)=M(N_x+k,N_x+k)+1
        end

        M(1,1)=1        //u à esquerda
        M(N_x,N_x)=1    //u à direita
        M(N_x+1,N_x+1)=1   //v à esquerda 
        M(2*N_x,2*N_x)=1  //v à direita
        
        b=zeros(2*N_x,1)
        b(1)= 0   //   u à esquerda
        b(N_x)= 1   // u à direita
        b(N_x+1)= 0   //   v à esquerda
        b(2*N_x)= 2   //   v à direita


    for j=1:N_t        
          
        for k = 2:N_x-1 //preenche linhas de 2 a N_x-1
            b(k)    =w(k)     +dt*theta * ((w(k-1)-2*w(k)+w(k+1))/dx^2 + alpha(1,1)*w(k)+alpha(2,1)*w(N_x+k))
            b(N_x+k)=w(N_x+k) +dt*theta * ((w(N_x+k-1)-2*w(N_x+k)+w(N_x+k+1))/dx^2 + alpha(2,1)*w(k)+alpha(2,2)*w(N_x+k))    
        end 
          
    
    //      w=lusolve(M,b)
            w=M\b
    end
y=M

endfunction
