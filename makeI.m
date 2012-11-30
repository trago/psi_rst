% Esta funcion genera una secuencia de interferogramas de PSI clasico.

% Date: 30 Nov 2012.
% Autor: Orlando Medina.

% ----- Inputs ------
% DC    : Termino de iluminacion de fondo.
% b     : Termino de contraste.
% phase : Fase de entrada de MxN.
% step  : Valor del paso entre interferogramas.
% k     : Numero de pasos.
% nv    : Nivel de varianza para el ruido.

% ----- Return ------
% I : Interferogramas de generados de MxNxk.



function [I] = makeI(DC,b,phase,step,k,nv)
    [M N] = size(phase);
    I = zeros(M,N);
    steps = (0:1:k-1)*step;

    for n=1:k
        noise = randn(M,N) * nv;
        I(:,:,n) =  DC + b .* cos( phase + steps(n) ) + noise;
    end
end