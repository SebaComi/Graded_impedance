function out = c1_qua_1D(alpha_medio,alpha_out,chi)
    n = @(r) r*0+1; % questo lo definisco solo per sport perché in questo caso non considero una lente

    c0 = @(c1) alpha_out* (c1+1).^2;
    alpha =  @(c1,x)  c0(c1) ./ (c1+x).^2;
    
    inte = @(c1) integral(@(x)  alpha(c1,x) ,0,1) - alpha_medio;
    
    c1 = fsolve(inte,-5);

    if chi == 0
        out = c0(c1);
    elseif chi == 1
        out = c1;
    elseif chi == 2
        out = @(r) alpha(c1,r) .* n(r);
    elseif chi == 3
        out = @(r) alpha(c1,r);
    end
end

