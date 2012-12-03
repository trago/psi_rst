%% Generando Interferogramas Sinteticos.
close all;
clear all;
M       = 512; % Number of rows of each interferogram.
N       = 512; % Number of columns of each interferogram.
k       = 5;   % Number of frames.
A       = 25;  % Amplitud para la fase tipo Peaks.

step    = pi/2; % Valor del paso.
nv      = 0.0; % Varianza del Ruido.

DC      = makeParabola(M,N,5);
rampa   = makeRampa(0.051,0.051,M,N);
phase   = makePeaks(N,M,A)+rampa;
b       = 1;
step_noise = 0.5;

[I,steps]       = makeI(DC,b,phase,step,step_noise,k,nv);

figure,imshow(I(:,:,1),[]),title('Interferograma de Entrada');

%% Inicializando parametros del metodo RST.

Muestreo = 6; % Numero de pixeles a satar para el muestreo.
iters1   = 25*2; % Numero de iteraciones para el metodo completo.
iters2   = 50; % Numero de iteraciones para el calculo de los pasos.
lambda   = 0; % Parametro de regulacizacion.
%% Inicializando parametros del metodo AIA.

iters = 20;
v     = 1;
Sk    = sin( v* (0:1:k-1) );
Ck    = cos( v* (0:1:k-1) );
Show  = 1; % 1 si se decea mostrar resultados parciales.

%% Aplicando metodos RST y AIA

% Aplicando algoritmo RST.
tic
[pasosRST f_RST] = RST(I,Sk,Ck,lambda,Muestreo,iters1,iters2,Show);
tRST = toc;
% Aplicando algoritmo AIA.
tic
[pasosAIA f_AIA] = AIA(I,Sk,Ck,iters,Show);
tAIA = toc;

%% Eliminando Piston de fase.
pasosRST = pasosRST-pasosRST(1);
Sk = sin(pasosRST);
Ck = cos(pasosRST);
[a1 f_RST] = MinCuaCpp(I,Sk,Ck);

pasosAIA = pasosAIA-pasosAIA(1);
Sk = sin(pasosAIA);
Ck = cos(pasosAIA);
[a1 f_AIA] = MinCuaCpp(I,Sk,Ck);

%% Mostrando Resultados.

SP_RST = angle(f_RST);
SP_AIA = angle(f_AIA);

figure,imshow(SP_RST,[]),title('fase Estimada RST');
figure,imshow(SP_AIA,[]),title('fase Estimada AIA');
figure,imshow(angle(exp(-1i*phase)),[]),title('fase Esperada');

disp('Estimados AIA');
disp(pasosAIA);

disp('Estimados RST');
disp(pasosRST);

disp('Esperados');
disp(steps);

disp('Error AIA');
disp(abs(steps - pasosAIA));

disp('Error RST');
disp(abs(steps - pasosRST));

figure;
imshow(I(:,:,1),[]),title('Interferograma de Entrada');