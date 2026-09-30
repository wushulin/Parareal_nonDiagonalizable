clc
clear;
global T Dt J  r Kmax Tol c flag_F
flag_F=2;
if flag_F==1
    c=2/2;
else
    c=1-1/sqrt(2);
end
%rng(0);
Kmax=100;
Err=zeros(1,Kmax);
Tol=1e-12;
Nt=2000;
Dt=0.1;
T=Nt*Dt;
J=8;
dt=Dt/J;
Nt=T/Dt;
r=5;
m0=100/r;
mm=m0*ones(r,1);
m=sum(mm);
Im=eye(m);
b=5;
bb=ones(r,1)*b;
%aa=bb+[0.1;1;10;100;4000];
%aa=bb+[0.1;4000;4000;4000;4000];
 aa=bb+[0.1;4000;4000;4000;4000];
aa=bb+[0.1;1;10;100;4000];
A=getA(mm, aa, bb, r);
U_Ini=sin((1:m)');
U_ref=zeros(m,Nt+1);
U_ref(:,1)=U_Ini;
inv_Ac=Im/(Im+Dt*A);
inv_Af=Im/(Im+c*dt*A);
for n=1:Nt
    Z0=U_ref(:,n);
    for j=1:J
        if flag_F==1
            Z1=inv_Af*((Im-(1-c)*dt*A)*Z0+dt*(c*getF(m,(n-1)*Dt+j*dt)+(1-c)*getF(m,(n-1)*Dt+(j-1)*dt)));
            Z0=Z1;
        else
           Z10=inv_Af*(Z0+c*dt*getF(m,(n-1)*Dt+(j-1)*dt+c*dt));
           hZ=((1-c)/c)*Z10+(1-(1-c)/c)*Z0;
           Z1=inv_Af*(hZ+c*dt*getF(m,(n-1)*Dt+j*dt));
           Z0=Z1;
        end
    end
    U_ref(:,n+1)=Z0;
end
Uk=random('Uniform',0,1,m,Nt+1);
Uk(:,1)=U_Ini;
Uk1=Uk;
k=1;
Err(k)=max(max(abs(Uk-U_ref)))/max(max(abs(U_ref)));
for k=2:Kmax
    Gk=zeros(m,Nt+1);
    Fk=zeros(m,Nt+1);
    for n=1:Nt
        Z0=Uk(:,n);
        for j=1:J
            if flag_F==1
                Z1=inv_Af*((Im-(1-c)*dt*A)*Z0+dt*(c*getF(m,(n-1)*Dt+j*dt)+(1-c)*getF(m,(n-1)*Dt+(j-1)*dt)));
                Z0=Z1;
            else
               Z10=inv_Af*(Z0+c*dt*getF(m,(n-1)*Dt+(j-1)*dt+c*dt));
               hZ=((1-c)/c)*Z10+(1-(1-c)/c)*Z0;
               Z1=inv_Af*(hZ+c*dt*getF(m,(n-1)*Dt+j*dt));
               Z0=Z1;
            end
        end
        Fk(:,n+1)=Z0;
        Gk(:,n+1)=inv_Ac*(Uk(:,n)+Dt*getF(m,n*Dt));
    end
    for n=1:Nt
        temp=inv_Ac*(Uk1(:,n)+Dt*getF(m,n*Dt));
        Uk1(:,n+1)=temp+Fk(:,n+1)-Gk(:,n+1);
    end
    Uk=Uk1;
    Err(k)=max(max(abs(Uk-U_ref)))/max(max(abs(U_ref)));
    fprintf('Error at %d-th iteration is %2.15d\n',k, Err(k));
    if Err(k)<=Tol
        break;
    end
end
if c==1/2
    semilogy(1:k,Err(1:k),'r-.o','MarkerSize',10);shg
    hold on;
elseif c==1-1/sqrt(2)
    semilogy(1:k,Err(1:k),'b-.+','MarkerSize',10);shg
else
    semilogy(1:k,Err(1:k),'k-.*','MarkerSize',10);shg
    hold off;
end
set(gca,'FontSize',15);
xlim([1,k]);
ylim([Tol,2*max(Err)]);
 xlabel('Iteration Index $k$','Interpreter','latex','FontSize',22);
 ylabel('Error','Interpreter','latex','FontSize',22);
 set(gca,'ytick',10.^(-12:2:10));
set(gca,'xtick',[1,2:3:k]);
leg=legend('use TR for $\mathcal{F}$','use SDIRK2 for $\mathcal{F}$','use backward Euler for $\mathcal{F}$');
set(leg,'fontsize',15,'interpreter','latex');
title(['$b=',num2str(b),'$'],'Interpreter','latex','FontSize',22);
grid on;
function F=getF(m,t)
F=zeros(m,1);
F(1)=0.5*(1+sin(2*pi*t)^2);
end
function A = getA(mm, aa, bb, r)
    n = sum(mm);
    nnz_total = sum(2*mm - 1) + (r-1);   
    I = zeros(nnz_total, 1);
    J = zeros(nnz_total, 1);
    V = zeros(nnz_total, 1);
    
    idx = 1;  
    row_offset = 0;
    
    for i = 1:r
        mi = mm(i);
        ai = aa(i);
        bi = bb(i);
        for j = 1:mi
            I(idx) = row_offset + j;
            J(idx) = row_offset + j;
            V(idx) = ai;
            idx = idx + 1;
        end
        for j = 1:mi-1
            I(idx) = row_offset + j + 1;
            J(idx) = row_offset + j;
            V(idx) = -bi;
            idx = idx + 1;
        end
        if i < r
            I(idx) = row_offset + mi + 1;
            J(idx) = row_offset + mi;
            V(idx) = -bb(i+1);
            idx = idx + 1;
        end
        
        row_offset = row_offset + mi;
    end

    I = I(1:idx-1);
    J = J(1:idx-1);
    V = V(1:idx-1);
    A = sparse(I, J, V, n, n);
end