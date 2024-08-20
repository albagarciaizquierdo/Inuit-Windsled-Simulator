function [X Error EXITO]=my_fzero(fun,X0,tol,Max_it,dh)

NT     = length(X0);
EXITO  = 0;
Error0 = 1000000;

for j=1:1:Max_it
    
    
    F0    = feval(fun,0,X0);
    Error = sum(abs(F0));
    
    if (j>1 && Error>Error0)
        X     = zeros(length(X0),1);
        EXITO =0;
        break
    end
    
    if Error<tol
        X     = X0;
        EXITO = 1;
        break
    else
  
       %Calculo del Jacobiano   
       J  = fun_jac_num(fun,0,X0,dh);
%        "Nuevo valor de x"
       X0 = X0 - J\F0;
%        pause
        
    
    end
   
    
end

if EXITO ==0
    X       = 0;
    Error   = 0;
end




end

