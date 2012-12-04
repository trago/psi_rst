function [steps f S C a] = RST(I,Sk,Ck,lambdaDC,lambdaSC,Muestreo,iters1,iters2,Show)

    q = Muestreo;

    Itmp = I(1:q:end,1:q:end,:);% Submuestreando los Interferogramas.

    [M N k] = size(Itmp);
    S   = zeros(M,N,k);
    C   = zeros(M,N,k);
    a   = ones(M,N);
    
%     fi = I(:,:,1)-I(:,:,3);
%     psi= I(:,:,4)-I(:,:,2);

    [M N k] = size(I);
    a1      = ones(M,N);
    f       = complex(a1,a1);
    
    for x=1:200
            [a1 f] = MinCuaReg(I,f,a1,Sk,Ck,7.5,1);
    end
    
    figure,
    for n=1:iters1
        
        for x=1:100
            [a1 f] = MinCuaReg(I,f,a1,Sk,Ck,7.5,1);
        end
%         a = a1;
%         ftmp = f;
%        [a1 f] = MinCuaCpp(I,Sk,Ck);
        ftmp = f(1:q:end,1:q:end);% Submuestreando fase calculada.
        a   = a1(1:q:end,1:q:end);

        for x=1:iters2
            [a S C] = gs_aCkSk(Itmp,ftmp,a,S,C,lambdaDC,lambdaSC);
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