%% Generando Interferogramas Sinteticos.
close all;
clear all;
M       = 512; % Number of rows of each interferogram.
N       = 512; % Number of columns of each interferogram.
k       = 5;   % Number of frames.
A       = 25;  % Amplitud para la fase tipo Peaks.

step    = pi/3; % Valor del paso.
nv      = 0.8; % Varianza del Ruido.

DC      = makeParabola(M,N,10);
rampa   = makeRampa(0.051,0.051,M,N);
phase   = makePeaks(N,M,A)+rampa;
b       = 1;
step_noise = 0.0;

[I,steps] = makeI(DC,b,phase,step,step_noise,k,nv);
steps = atan2(sin(steps),cos(steps));


%% Inicializando parametros para Minimos Cuadrados Regularizado.

iters = 50;
lambdaDC = 1; % Parametro de regulacizacion para el DC.
lambdaf  = 10; % Parametro de regulacizacion para campo complejo.

Phi   = zeros(M,N);
Psi   = zeros(M,N);
a     = ones(M,N);
f     = complex(Phi,Psi);

Sk    = sin(steps);
Ck    = cos(steps);

%% Aplicando metodo Minimos Cuadrados Regularizado
for n=1:iters
    [a f] = MinCuaReg(I,f,a,Sk,Ck,lambdaf,lambdaDC);
%     drawnow expose update
%     imshow(angle(f),[]),title(['Fase Min Cua Reg en Iteracion: ' num2str(n)]);
end
fase = exp(1j*phase);
figure,imshow(I(:,:,1),[]),title('Interferograma de Entrada');
figure,imshow(a,[]),title('DC Calculado');
figure,imshow(angle(fase),[]),title('Fase Esperada');
figure,imshow(angle(f),[]),title('Fase Min Cua Reg');

