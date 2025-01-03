function out = a1_exp_L(alpha_medio,alpha_out,chi)
    n = @(r) sqrt(2 - r.^2);

    a0 = @(a1) alpha_out * exp(-a1);
    alpha = @(a1,r) a0(a1)*exp(a1*r);

    inte = @(a1) integral(@(r) alpha(a1,r) .*n(r) .*(2*pi*r)  ,0,1) -  alpha_medio *integral(@(r) n(r) .*(2*pi*r) ,0,1);
    % inte = @(a1) (a1 - 1 + exp(-a1))./a1.^2 - alpha_medio/2;
    
    a1 = fsolve(inte,3);
    % Ricorda che a1 in realtà è a1*R, quindi poi devi riscalarlo
    if chi == 0
        out = a0(a1);
    elseif chi == 1
        out = a1;
    elseif chi == 2
        out = @(r) alpha(a1,r) .* n(r);
    elseif chi == 3
        out = @(r) alpha(a1,r);
    end
end

% alpha_medio = 0.3;
% a1 = a1_exp(alpha_medio); a0 = exp(-a1);

% Controllo integrale
% integral(@(r) 2*pi*a0* exp(a1*r).*r ,0,1) / pi - alpha_medio

% figure
% rr = linspace(0,1);
% plot(rr,a0 * exp(a1*rr))

% alpha_medio = 0.2;
% a1 = a1_qua_BH(alpha_medio,1,1/3,1);
% a0 = a1_qua_BH(alpha_medio,0,1/3,1);
% rho = a1_qua_BH(alpha_medio,2,1/3,1);
% rr = linspace(1/3,1,100);
% plot(rr,rho(rr))
% hold on
% a1_exp = a1_exp_BH(alpha_medio,1,1/3);
% a0_exp = exp(-a1_exp*1);
% rho_exp = @(r) a0_exp*exp(a1_exp*r) .* 1^2./r.^2;
% rr = linspace(1/3,1,100);
% plot(rr,rho_exp(rr))

