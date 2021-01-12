function [ps] = doblez(B,g,Delta)
% Técnica del doblez para el problema cuadrático de región de confianza
% Suponemos que B ya es simétrica positiva definida
% Min (1/2)*p'*B*p + g'*p + f(x)
% sa     ||p||<=Delta

pN= -B\g; %dirección de Newton
pC= -((g'*g)/(g'*B*g))*g;

if norm(pN)<= Delta
    ps=pN;
else
    if norm(pC)>= Delta
        ps=(Delta/norm(pC))*(pC);
    else
        %Resolver la ecuación de segundo grado
        % ||pC+ t*(pN-pC)||^2 - (Delta)^2 = 0
        % La cual equivale a:
        % (pC+t*u)'(pC-t*u)-(Delta)^2=0 con u=pN-pC
        % a*(t^2)+b*t+c=0
        
        u=pN-pC; 
        a=u'*u; b=2*pC'*u; c=pC'*pC-Delta^2;
        t=roots([a b c]); %Obtenemos las raíces del polinomio
        %Recuerda que obtienes dos raíces, nos interesa la +
        if t(1)>0
            ts=t(1);
        else
            ts=t(2);
            
        end
        
        ps=pC+ts*u;
    end
end

end
        
        