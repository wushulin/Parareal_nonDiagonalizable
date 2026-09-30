clc;
clear;
global m DT lambda Jc Jf Nt A T theta_c theta_f
theta_c=1;
theta_f=1;
Tol=1e-14;
Kmax=50;
m=1;
dem=m;
DT=0.1;
Jc=1;
Jf=40;
dt_c=DT/Jc;
dt_f=DT/Jf;
lambda=0.05;
Nt=200;T=Nt*DT;

 
A = lambda * eye(m) + diag(ones(m-1,1),1);
r_c=stabFun_G(DT*A);
r_f=stabFun_F(dt_f*A);
U0=rand(m,1);
U_ref=zeros(m,Nt+1);
U_ref(:,1)=U0;
t=0:DT:T;
for n=1:Nt
    u_n=U_ref(:,n);
    for j=1:Jf
        u_n=F(u_n,t(n)+j*dt_f,dt_f);
    end
    U_ref(:,n+1)=u_n;
end
Umax_ref=max(max(abs(U_ref))); 

Error=zeros(1,Kmax);
Uk=random('Uniform',-1,1,dem,Nt+1);
Uk(:,1)=U0;

Uk1=Uk;
k=1;
Error(k)=max(max(abs(Uk-U_ref)))/Umax_ref;
fprintf('Initial Error is: %2.15f\n',Error(k));
for k=2:Kmax
    Gk=zeros(dem,Nt+1);
    Fk=zeros(dem,Nt+1);
    for n=1:Nt
        u_n=Uk(:,n);
        for j=1:Jc
            u_n=G(u_n, t(n)+j*dt_c,dt_c);
        end
        Gk(:,n+1)=u_n;
        u_n=Uk(:,n);
        for j=1:Jf
            u_n=F(u_n, t(n)+j*dt_f,dt_f);
        end
        Fk(:,n+1)=u_n;
    end
    for n=1:Nt
        temp=Uk1(:,n);
        for j=1:Jc
            temp=G(temp, t(n)+j*dt_c,dt_c);
        end
        Uk1(:,n+1)=temp+Fk(:,n+1)-Gk(:,n+1);
    end
    Uk=Uk1;
    Error(k)=max(max(abs(Uk-U_ref)))/Umax_ref;
    fprintf('Error at %d-th iteration is: %2.15f\n',k,Error(k));
    if Error(k)<=Tol
        break;
    end
end
semilogy(1:k,Error(1:k),'r-.o','Markersize',12);shg
set(gca,'fontsize',15)
function val=G(un,tn,h)
global m A
Ix=eye(m);
val=(Ix+h*A)\(un+h*sf(tn)); % u'+Au=f: u_{n+1}-u_n+dt*A*u_{n+1}=dt*f(t_n)
end

function val=F(un,tn,h)
global m A
Ix=eye(m);
val=(Ix+h*A)\(un+h*sf(tn));  
end
function val = sf(t)
global m T
val=zeros(m,1);
%val(m,1)= 1/(1+t^10); 
end
function val=stabFun_G(z)
global theta_c  m
Im=speye(m);
val=(Im+theta_c*z)\(Im-(1-theta_c)*z);
end
function val=stabFun_F(z)
global theta_f  m
Im=speye(m);
val=(Im+theta_f*z)\(Im-(1-theta_f)*z);
end