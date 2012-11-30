% Esta funcion genera una parabola de amplitud A.
% M: Numero de Renglones.
% N: Numero de Columnas.
% A: Amplitud de la parabola.
% Autor: Orlando Medina.
% fecha: 30 Noviembre 2012.

function [Parabola] = makeParabola(M,N,A)
    [X Y]    = meshgrid(1:N,1:M); 
    b        = (X-M/2).^2 + (Y-N/2).^2;
    Parabola = (-b / (max(max(b))/1) +1 )*A;
end