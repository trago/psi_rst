function ModaPasos = get_ModaPasos(S,C)

    [M N k] = size(S);
    pasos = atan2(-S,C);
    
%     for x=1:M
%         for y=1:N
%             pasos(x,y,:) = AntiAliasing(pasos(x,y,:));
%         end
%     end
    for n=1:k
        mi = min(min(pasos(:,:,n)));
        ma = max(max(pasos(:,:,n)));
        tmp = pasos(:,:,n);

        [values xPlace] = hist(tmp(:),mi:.0001:ma);
        [mModa mPlace]  = max(values);
        ModaPasos(n)    = xPlace(mPlace);
    end
    ModaPasos = AntiAliasing(ModaPasos);


end