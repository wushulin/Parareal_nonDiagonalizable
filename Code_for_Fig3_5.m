clc;
clear;
Nt=8;
  
lam_vec = linspace(0, 6, 2004);     % range of λ 
dt_vec  = linspace(0, 2, 2000);      % range of ΔT 
cond_my=zeros(length(lam_vec),length(dt_vec));
 
for i1= 1:length(lam_vec)
    for i2=1:length(dt_vec)
        lam = lam_vec(i1);
        dt=dt_vec(i2);
        z   = lam * dt;
        beta_val = 1 / (1+z-dt);
 
 
        % --- compute  alpha 
        termSum =0;
        for j = 0:3
            coeff = dt^j / factorial(j);
            diffTerm = abs( exp(-z) - factorial(j) * (1+z).^(-j-1) );
            termSum = termSum + coeff * diffTerm;
        end
        alpha_val= termSum + dt^4 * (1+z)^(-4) / (1 - dt / (1+z));
 
        num = (1 - beta_val) * (1 - beta_val^Nt);
        den = 1 - beta_val^Nt - Nt * (1 - beta_val) * beta_val^(Nt-1);
        RHS = num / den; 
     
        if alpha_val<=RHS
            cond_my(i1,i2) =true;
        else
            cond_my(i1,i2) =false;
        end
        if dt>=0.85 && lam<=0.55
            cond_my(i1,i2)=false;
        end
    end
end
 
figure(1);
pcolor(lam_vec, dt_vec, double(cond_my'));
shading flat;
set(gca,'fontsize',16)
xlabel('$\lambda$','Interpreter','latex','FontSize',22)
ylabel('$\Delta T$','Interpreter','latex','FontSize',22);
title(['$N_t=',num2str(Nt),'$'],'Interpreter','latex','FontSize',22)
colormap([1 1 1;0.75 0.88 0.98]);    
grid on;
set(gca,'xtick',[0,1,2]); 
hold on;
contour(lam_vec, dt_vec, double(cond_my'), [0.5 0.5], 'k', 'LineWidth', 1.5);
hold off;

 