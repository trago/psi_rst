%% Generando Interferogramas Sinteticos.
close all;
clear all;
M       = 512; % Number of rows of each interferogram.
N       = 512; % Number of columns of each interferogram.
k       = 5;   % Number of frames.
A       = 25;  % Amplitud para la fase tipo Peaks.

step    = pi/2; % Valor del paso.
nv      = .80; % Varianza del Ruido.

DC      = makeParabola(M,N,15);
rampa   = makeRampa(0.051,0.051,M,N);
phase   = makePeaks(N,M,A)+rampa;
b       = 1;
step_noise = 0.0;

[I,steps] = makeI(DC,b,phase,step,step_noise,k,nv);
%steps = atan2(sin(steps),cos(steps));


%% Inicializando parametros para Minimos Cuadrados Regularizado.

iters = 100;
lambdaDC = 0; % Parametro de regulacizacion para el DC.
lambdaf  = 0; % Parametro de regulacizacion para campo complejo.

Phi   = ones(M,N);
Psi   = ones(M,N);
a     = ones(M,N);
f     = complex(Phi,Psi);
s = 1;
Sk    = sin(s*(0:1:k-1));
Ck    = cos(s*(0:1:k-1));

%% Aplicando metodo Minimos Cuadrados Regularizado
for n=1:iters
    [a f] = MinCuaReg(I,f,a,Sk,Ck,lambdaf,lambdaDC);
    disp(n);
%     drawnow expose update
%     imshow(angle(f),[]),title(['Fase Min Cua Reg en Iteracion: ' num2str(n)]);
end


fase = angle(exp(1j*phase));
SP   = angle(f);
errorFase = mean2(abs(fase-SP));
disp('Pasos');
disp(steps);
disp('Error de fase');
disp(errorFase);
figure,imshow(I(:,:,1),[]),title('Interferograma de Entrada');
figure,imshow(a,[]),title('DC Calculado');
figure,imshow(fase,[]),title('Fase Esperada');
figure,imshow(SP,[]),title('Fase Min Cua Reg');

