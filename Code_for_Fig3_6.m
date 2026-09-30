clc;
clear;
T=50;
plot_region(T);
function plot_region(T)  
lam_vec=1.02.^(-1000:107); % range of λ
dt_vec=1.02.^(-1000:40); % range of  ΔT 
[Lam, DT] = meshgrid(lam_vec, dt_vec);
 
alpha = nan(size(Lam));
beta  = nan(size(Lam));
cond  = false(size(Lam));

for i = 1:numel(Lam)
    lam = Lam(i);
    dt  = DT(i);
    z   = lam * dt;
    Nt=(T/dt);
    % compute β 
    denom_beta = 1 + dt * (lam - 1);   % = 1 + z - dt
    if denom_beta <= 0
        continue;    
    end
    beta_val = 1 / denom_beta;
    beta(i)  = beta_val;
    
   w=1+z;
    % --- compute alpha ---
    termSum =0;
    for j = 0:3
        coeff = dt^j / factorial(j);
        diffTerm = abs( exp(-z) - factorial(j) * w.^(-j-1) );
        termSum = termSum + coeff * diffTerm;
    end
 
    denom_alpha = 1 - dt / (1+z);    
    if denom_alpha <= 0
        continue;
    end
    alpha_val= termSum + dt^4 * (1+z)^(-4) / denom_alpha;
    alpha(i) =alpha_val;
    
    if abs(beta_val - 1) < 1e-10       
        RHS = 2 / (Nt - 1);
    else
        num = (1 - beta_val) * (1 - beta_val^Nt);
        den = 1 - beta_val^Nt - Nt * (1 - beta_val) * beta_val^(Nt-1);
        if abs(den) < 1e-10
            RHS = inf;    
        else
            RHS = num / den;
        end
    end

    if alpha(i) <= RHS
        cond(i) = true;
    end
end

figure(1);
pcolor(Lam, DT, double(cond));
shading flat;
set(gca,'fontsize',16)
xlabel('$\lambda$','Interpreter','latex','FontSize',22)
ylabel('$\Delta T$','Interpreter','latex','FontSize',22);
title(['$T=',num2str(T),'$'],'Interpreter','latex','FontSize',22)
colormap([1 1 1; 0.75 0.88 0.98]);    
grid on;
set(gca,'xtick',[0,1,2]);
hold on;
contour(Lam, DT, double(cond), [0.5 0.5], 'k', 'LineWidth', 1.5);
hold off;
xlim([0,6]);
ylim([0,2]);
end