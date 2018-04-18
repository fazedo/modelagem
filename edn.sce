function y = fun(x)
    y = exp(-x)
endfunction

u = linspace(0,1,101)
v = fun(u)
