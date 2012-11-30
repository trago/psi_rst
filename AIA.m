function [steps f] = AIA(I,Sk,Ck,iters,Show)

    figure,
    for n=1:iters
        [DC f] = MinCuaCpp(I,Sk,Ck);
        [Sk Ck]  = getStepsAIACpp(I,f);

        steps = AntiAliasing(atan2(Sk,Ck));

        % Mostrando Resultados Parciales si Show vale 1.
        if(Show == 1)
            faseMC = angle(f);
            drawnow expose
            imshow(faseMC,[]),title(['fase AIA en Iteracion: ' num2str(n)]);
            disp([n steps]);
        end

    end
end