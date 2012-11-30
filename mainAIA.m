%% Generando Interferogramas Sinteticos.
close all;
clear all;
M       = 512; % Number of rows of each interferogram.
N       = 512; % Number of columns of each interferogram.
k       = 4;   % Number of frames.
A       = 25;  % Amplitud para la fase tipo Peaks.

step    = pi/2; % Valor del paso.
nv      = 0.0; % Varianza del Ruido.

DC      = makeParabola(M,N,2);
rampa   = makeRampa(0.051,0.051,M,N);
phase   = makePeaks(N,M,A)+rampa;
b       = 1;

I       = makeI(DC,b,phase,step,k,nv);

%% Inicializando parametros del metodo RST.

iters = 50;
v     = pi/2;
Sk    = sin( v* (0:1:k-1) );
Ck    = cos( v* (0:1:k-1) );

%% Aplicando algoritmo AIA.

tic
[pasos f] = AIA(I,Sk,Ck,iters,1);
tAIA = toc;

%% Mostrando Resultados.

SP = angle(f);
figure,imshow(SP,[]),title('fase Estimada');
figure,imshow(angle(exp(-1i*phase)),[]),title('fase Esperada');

disp('Tiepo de Procesamiento');
disp(tAIA);

disp('Estimados');
disp(pasos-pasos(1));

disp('Esperados');
disp(step*(0:k-1));
