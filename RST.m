function [steps f] = RST(I,Sk,Ck,lambda,Muestreo,iters1,iters2,Show)

    q = Muestreo;
    Itmp = I(1:q:end,1:q:end,:);% Submuestreando los Interferogramas.

    [M N k] = size(Itmp);
    S   = zeros(M,N,k);
    C   = zeros(M,N,k);
    a   = zeros(M,N);
    
%     fi = I(:,:,1)-I(:,:,3);
%     psi= I(:,:,4)-I(:,:,2);
%     f  = complex(fi,psi);
    figure,
    for n=1:iters1

       [a1 f] = MinCuaCpp(I,Sk,Ck);
       ftmp = f(1:q:end,1:q:end,:);% Submuestreando fase calculada.

        for x=1:iters2
            [a S C] = gs_aCkSk(Itmp,ftmp,a,S,C,lambda);
        end

        steps = get_ModaPasos(S,C);
        Sk = sin(steps);
        Ck = cos(steps);
        
        % Mostrando Resultados Parciales si Show vale 1.
        if (Show == 1)
            drawnow expose
            imshow(angle(f),[]),title(['Fase RST en Iteracion: ' num2str(n)]);
            disp([n steps]);
        end
    end

end