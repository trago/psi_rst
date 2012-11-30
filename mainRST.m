%% Generando Interferogramas Sinteticos.
close all;
clear all;
M       = 512; % Number of rows of each interferogram.
N       = 512; % Number of columns of each interferogram.
k       = 4;   % Number of frames.
A       = 25;  % Amplitud para la fase tipo Peaks.

step    = pi/4; % Valor del paso.
nv      = 0.0; % Varianza del Ruido.

DC      = makeParabola(M,N,2);
rampa   = makeRampa(0.051,0.051,M,N);
phase   = makePeaks(N,M,A)+rampa;
b       = 1;

I       = makeI(DC,b,phase,step,k,nv);
figure,imshow(I(:,:,1),[]);

%%
%clear all
%load vela.mat
%load placa.mat
% load ExperimentalResultsSorted.mat
% I(:,:,1) = I2(:,:,1);
% I(:,:,2) = I2(:,:,2);
% I(:,:,3) = I2(:,:,3);
% I(:,:,4) = I2(:,:,4);
% I(:,:,5) = I2(:,:,5);
% I(:,:,6) = I2(:,:,6);
% I(:,:,7) = I2(:,:,7);
% I = I(:,1:end,:);
% [M N k]=size(I);
% 
% figure,imshow(I(:,:,1),[]);
% h = fspecial('average',11);
% %I = removeDC(I,105);
% 
% %figure,imshow(I(:,:,1),[]);
% I = imfilter(I,h,'replicate');


%% Inicializando parametros del metodo RST.

Muestreo = 8; % Numero de pixeles a satar para el muestreo.

Sk = sin((0:1:k-1)*(1));
Ck = cos((0:1:k-1)*(1));

iters1   = 20; % Numero de iteraciones para el metodo completo.
iters2   = 50; % Numero de iteraciones para el calculo de los pasos.
lambda   = 0; % Parametro de regulacizacion.

%% Aplicando algoritmo RST.
tic
[pasos f] = RST(I,Sk,Ck,lambda,Muestreo,iters1,iters2,1);
tRST = toc;

%% Mostrando Resultados.

SP = angle(f);
figure,imshow(SP,[]),title('fase Estimada');
figure,imshow(angle(exp(-1i*phase)),[]),title('fase Esperada');

disp('Tiepo de Procesamiento');
disp(tRST);

disp('Estimados');
disp(pasos-pasos(1));

disp('Esperados');
disp(step*(0:k-1));

    
    


