clc
clear;
global T Dt J  Kmax Tol 
Kmax=30;
Err=zeros(1,Kmax);
Tol=1e-12;
Nt=30;
b0=1;
flag=2;
if flag==1
    Dt=0.5*b0;
    a=1/b0;
    b=1;
else
    Dt=0.5;
    a=1;
    b=b0;
end

T=Nt*Dt;
J=100;
dt=Dt/J;
Nt=T/Dt;
m=100;
Im=eye(m);

A=toeplitz([a;zeros(m-1,1)],[a,b,zeros(1,m-2)]);
U_Ini=random('Uniform',0,1,m,1);
U_ref=zeros(m,Nt+1);
U_ref(:,1)=U_Ini;
inv_Ac=Im/(Im+Dt*A);
inv_Af=Im/(Im+dt*A);
for n=1:Nt
    Z0=U_ref(:,n);
    for j=1:J
        getF=g((n-1)*Dt+j*dt)*ones(m,1);
        Z1=inv_Af*(Z0+dt*getF);
        Z0=Z1;
    end
    U_ref(:,n+1)=Z0;
end
Uk=random('Uniform',0,1,m,Nt+1);
Uk(:,1)=U_Ini;
Uk1=Uk;
for k=1:Kmax
    Gk=zeros(m,Nt+1);
    Fk=zeros(m,Nt+1);
    for n=1:Nt
        Z0=Uk(:,n);
        for j=1:J
            getF=g((n-1)*Dt+j*dt)*ones(m,1);
            Z1=inv_Af*(Z0+dt*getF);
            Z0=Z1;
        end
        Fk(:,n+1)=Z0;
        getF=g(n*Dt)*ones(m,1);
        Gk(:,n+1)=inv_Ac*(Uk(:,n)+Dt*getF);
    end
    for n=1:Nt
        getF=g(n*Dt)*ones(m,1);
        temp=inv_Ac*(Uk1(:,n)+Dt*getF);
        Uk1(:,n+1)=temp+Fk(:,n+1)-Gk(:,n+1);
    end
    Uk=Uk1;
    Err(k)=max(max(abs(Uk-U_ref)))/max(max(abs(U_ref)));
    fprintf('Error at %d-th iteration is %2.15d\n',k, Err(k));
    if Err(k)<=Tol
        break;
    end
end
if flag==2
    if b0==2
            semilogy(1:k,Err(1:k),'r-.o','MarkerSize',10);shg
            hold on;
        elseif b0==1.5
            semilogy(1:k,Err(1:k),'b-.o','MarkerSize',10);shg
        else
            semilogy(1:k,Err(1:k),'k-.o','MarkerSize',10);shg
            hold off
    end
else
    if b0==2
        semilogy(1:k,Err(1:k),'r-.*','MarkerSize',10);shg
        hold on;
    elseif b0==1.5
        semilogy(1:k,Err(1:k),'b-.*','MarkerSize',10);shg
    else
        semilogy(1:k,Err(1:k),'k-.*','MarkerSize',10);shg
        hold off
    end
end

 set(gca,'fontsize',16);
 xlabel('Iteration Index $k$','Interpreter','latex','FontSize',22);
 ylabel('Error','Interpreter','latex','FontSize',22);
 xlim([1,k]);
 set(gca,'ytick',10.^(-12:2:10));
 ylim([Tol,5*1.2315e+03]);
set(gca,'xtick',[1,2:2:30]);
leg=legend('$b=2$','$b=1.5$','$b=1$');
set(leg,'fontsize',16,'interpreter','latex');
function val=g(t)
global T
val=cos(5.5*acos(2*t-T)/T);
end
