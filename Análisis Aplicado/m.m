function fx = m(p,B,g,x)
%Funci�n correspondiente al modelo cuadr�tico
%In
% p - vector de longitud 2
% B - matriz Hessiana evaluada en un punto
% g - gradiente de la funci�n evaluada en un punto
%Out
% fx.- n�mero real con el valor de la funci�n.

fx=(1/2)*p'*B*p+g'*p+feval('rosenbrock',x);

