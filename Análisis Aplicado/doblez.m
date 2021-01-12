function [ps] = doblez(B,g,Delta)
% T�cnica del doblez para el problema cuadr�tico de regi�n de confianza
% Suponemos que B ya es sim�trica positiva definida
% Min (1/2)*p'*B*p + g'*p + f(x)
% sa     ||p||<=Delta

pN= -B\g; %direcci�n de Newton
pC= -((g'*g)/(g'*B*g))*g;

if norm(pN)<= Delta
    ps=pN;
else
    if norm(pC)>= Delta
        ps=(Delta/norm(pC))*(pC);
    else
        %Resolver la ecuaci�n de segundo grado
        % ||pC+ t*(pN-pC)||^2 - (Delta)^2 = 0
        % La cual equivale a:
        % (pC+t*u)'(pC-t*u)-(Delta)^2=0 con u=pN-pC
        % a*(t^2)+b*t+c=0
        
        u=pN-pC; 
        a=u'*u; b=2*pC'*u; c=pC'*pC-Delta^2;
        t=roots([a b c]); %Obtenemos las ra�ces del polinomio
        %Recuerda que obtienes dos ra�ces, nos interesa la +
        if t(1)>0
            ts=t(1);
        else
            ts=t(2);
            
        end
        
        ps=pC+ts*u;
    end
end

end
        
        