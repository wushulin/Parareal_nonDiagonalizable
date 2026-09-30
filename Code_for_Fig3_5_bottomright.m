clc;
clear;

lam_vec = linspace(0, 6, 2000);          % range of λ 
dt_vec  = linspace(0, 2, 2000);          % range of ΔT 
[Lam, DT] = meshgrid(lam_vec, dt_vec);

 
cond = false(size(Lam));
 
for i = 1:numel(Lam)
    lam = Lam(i);
    dt  = DT(i);
    z   = lam * dt;
    
    % coumpute β 
    denom_beta = 1 + dt * (lam - 1);   % = 1 + z - dt
    if denom_beta <= 0
        continue;    
    end
    beta_val = 1 / denom_beta;
    
    % compute α
    s = 0;
    for j = 0:3
        term = (dt^j / factorial(j)) * abs( exp(-z) - factorial(j) * (1+z)^(-j-1) );
        s = s + term;
    end
    
    denom_alpha = 1 - dt / (1+z);    
    if denom_alpha <= 0
        continue;
    end
    alpha_val = s + dt^4 * (1+z)^(-4) / denom_alpha;
 
    if alpha_val + beta_val < 1
        cond(i) = true;
    end
end
 
cond(dt_vec == 0, :) = false;    

figure(1);
pcolor(Lam, DT, double(cond));
set(gca,'fontsize',16)
xlabel('$\lambda$','Interpreter','latex','FontSize',22)
ylabel('$\Delta T$','Interpreter','latex','FontSize',22); 
colormap([1 1 1; 0.75 0.88 0.98]);   
axis tight;
grid on;
yticks([8e-4,0.2:0.2:2])
yticklabels({'0', '0.2', '0.4', '0.6', '0.8','1','1.2','1.4','1.6','1.8','2'})
shading flat; 
colormap([1 1 1; 0.75 0.88 0.98]);   
axis tight;
grid on;
 
hold on;
contour(Lam, DT, double(cond), [0.5 0.5], 'k', 'LineWidth', 1.5);
hold off;
 