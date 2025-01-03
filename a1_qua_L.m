function out = a1_qua_L(alpha_medio,alpha_out,chi)
    n = @(r) sqrt(2 - r.^2);
    % N = @(r1) integral(@(r) n(r)./r ,0,r1);
    N = @(r) n(r) - sqrt(2)*atanh(n(r)/sqrt(2));
    
    a0 = @(a1,r) alpha_out*(a1+N(r)).^2;
    alpha =  @(a1,r) a0(a1,1) ./ (a1+N(r)).^2;
    inte = @(a1) integral(@(r) alpha(a1,r) .* n(r) .* (2*pi*r)  ,0,1) - alpha_medio *integral(@(r) n(r) .*(2*pi*r) ,0,1);
    
    a1 = fsolve(inte,-0.5);

    if chi == 0
        out = a0(a1,1);
    elseif chi == 1
        out = a1;
    elseif chi == 2
        out = @(r) alpha(a1,r) .* n(r);
    elseif chi == 3
        out = @(r) alpha(a1,r);
    end
    % Ricorda che a1 in realtà è a1*R, quindi poi devi riscalarlo
end

% alpha_medio = 0.3;
% a1 = a1_exp(alpha_medio); a0 = exp(-a1);

% Controllo integrale
% integral(@(r) 2*pi*a0* exp(a1*r).*r ,0,1) / pi - alpha_medio

% figure
% rr = linspace(0,1);
% plot(rr,a0 * exp(a1*rr))

