clc;
clear;
r=1.5;
%dtt=0.9.^(2:100);
dtt=[1/r:-0.01:0.25,0.92.^(17:110)];
dtt=[1/4,1/128,1/4096];
LL=0.92.^(-175:110);
%LL=r.^(-20:23);
alpha_scalar=zeros(length(dtt),length(LL));
beta_scalar=zeros(length(dtt),length(LL));
alpha_sys=zeros(length(dtt),length(LL));
beta_sys=zeros(length(dtt),length(LL));
for j1=1:length(dtt)
    dt=dtt(j1);
    DeltaT=dt;
    for j2=1:length(LL)
        lambda=LL(j2);
        z=dt*lambda; 
        w = 1 + z;    
        termSum = zeros(size(z));
        for j = 0:3
            coeff = DeltaT^j / factorial(j);
            diffTerm = abs( exp(-z) - factorial(j) * w.^(-j-1) );
            termSum = termSum + coeff * diffTerm;
        end
      
        secondTerm = (DeltaT^4) * w.^(-4) ./ (1 - DeltaT * w.^(-1));
        
        alpha_sys(j1,j2)= termSum + secondTerm;
        alpha_scalar(j1,j2)= abs(exp(-z)-1/(1+z));

       % --- compute beta ---
        beta_sys(j1,j2)= 1 / (1+z-DeltaT);
        beta_scalar(j1,j2)= 1 / (1+z);
    end
end
% figure(1);
% mesh(LL,dtt,beta_sys);shg
% hold on;
% mesh(LL,dtt,beta_scalar);shg
% hold off;
% set(gca, 'XScale', 'log');
% set(gca, 'YScale', 'log');
% set(gca,'fontsize',16);
% xlabel('$\lambda$','FontSize',21,'Interpreter','latex');
% ylabel('$\Delta T$','FontSize',21,'Interpreter','latex');
% zlabel('$\beta$','FontSize',21,'Interpreter','latex');
% xlim([min(LL),max(LL)]);
% ylim([min(dtt),max(dtt)]);
% set(gca,'xtick',10.^(-4:4));
% set(gca,'ytick',10.^(-4:4));
% figure(2);
% mesh(LL,dtt,alpha_sys);shg
% hold on;
% mesh(LL,dtt,alpha_scalar);shg
% hold off;
% set(gca, 'XScale', 'log');
% set(gca, 'YScale', 'log');
% set(gca,'fontsize',16);
% xlim([min(LL),max(LL)]);
% ylim([min(dtt),max(dtt)]);
% xlabel('$\lambda$','FontSize',21,'Interpreter','latex');
% ylabel('$\Delta T$','FontSize',21,'Interpreter','latex');
% zlabel('$\alpha$','FontSize',21,'Interpreter','latex');
% set(gca,'xtick',10.^(-4:4));
% set(gca,'ytick',10.^(-4:4));

figure(3);
plot(LL,beta_sys(1,:),'r',LL,beta_sys(2,:),'b',LL,beta_sys(3,:),'k','LineWidth',1.5);shg
hold on;
plot(LL,beta_scalar(1,:),'r--',LL,0.995*beta_scalar(2,:),'b--',LL,beta_scalar(3,:),'k--','LineWidth',1.5);shg
hold off;
set(gca, 'XScale', 'log');
%set(gca, 'YScale', 'log');
set(gca,'fontsize',16);
xlabel('$\lambda$','FontSize',21,'Interpreter','latex'); 

set(gca,'xtick',10.^(-5:5));
xlim([1e-4,1e+5]);

figure(4);
plot(LL,alpha_sys(1,:),'r',LL,alpha_sys(2,:),'b',LL,alpha_sys(3,:),'k','LineWidth',1.5);shg
hold on;
plot(LL,alpha_scalar(1,:),'r--',LL,0.983*alpha_scalar(2,:),'b--',LL,alpha_scalar(3,:),'k--','LineWidth',1.5);shg
hold off;
set(gca, 'XScale', 'log');
%set(gca, 'YScale', 'log');
set(gca,'fontsize',16);
xlabel('$\lambda$','FontSize',21,'Interpreter','latex'); 
xlim([1e-4,1e+5]);
set(gca,'xtick',10.^(-5:5));
%set(gca,'ytick',10.^(-4:4));