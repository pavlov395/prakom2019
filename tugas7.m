function tugas7
%Tugas 7
clear all;clc
%Data
k=0.1;
h=0.0007;
D=1;
Ts=400;
Tu=50;
L=20;
N=101;
eps=0.001;
tol=0.00001;
% Boundaary Condition
% x=0 T=Ts
% x=L dT/dX |x=L = -h/k * (T|x=L - Tu)
%Initial Condition
% Trial dT/dX |x=0
x(1)=0;
T(1)=Ts;

%Persamaan Diferensial
%d^2T/dx^2 - 4h/kD ( T - Tu ) =0
% pemisalan
%dT/dx = Z ,sehingga dZ/dz=d^2T/dx^2  PD menjadi :
%dZ/dx = 4h/kD ( T - Tu )          (1)
%dT/dx = Z                         (2)

% PD (1) dan (2) diselesaikan secara simultan
% Runge Kutta
Zold=0;
er=100;
stepX=(L-0)/(N-1);
while er>tol
    for i=1:N
        Z(1)=runge(Zold);
        x(i+1)=x(i)+stepX;
        %Mencari k1,k2,k3,k4 untuk mencari rata2 slope
        k1h=frk(x(i)        ,Z(i)            ,T(i)            );
        k2h=frk(x(i)+stepX/2,Z(i)+stepX/2*k1h(1),T(i)+stepX/2*k1h(2));
        k3h=frk(x(i)+stepX/2,Z(i)+stepX/2*k2h(1),T(i)+stepX/2*k2h(2));
        k4h=frk(x(i)+stepX  ,Z(i)+stepX  *k3h(1),T(i)+stepX  *k3h(2));
        %Mencari nilai Z dan T step berikutnya
        Z(i+1)=Z(i)+1/6*(k1h(1)+2*k2h(1)+2*k3h(1)+k4h(1))*stepX;
        T(i+1)=T(i)+1/6*(k1h(2)+2*k2h(2)+2*k3h(2)+k4h(2))*stepX;
        b1(i)=Z(i)-h/k*(T(i)-Tu);
        xfinal(i)=x(i);
        Tfinal(i)=T(i);
    end
    for i=1:N
        Z(1)=runge(Zold+eps);
        x(i+1)=x(i)+stepX;
        %Mencari k1,k2,k3,k4 untuk mencari rata2 slope
        k1h=frk(x(i)        ,Z(i)            ,T(i)            );
        k2h=frk(x(i)+stepX/2,Z(i)+stepX/2*k1h(1),T(i)+stepX/2*k1h(2));
        k3h=frk(x(i)+stepX/2,Z(i)+stepX/2*k2h(1),T(i)+stepX/2*k2h(2));
        k4h=frk(x(i)+stepX  ,Z(i)+stepX  *k3h(1),T(i)+stepX  *k3h(2));
        %Mencari nilai Z dan T step berikutnya
        Z(i+1)=Z(i)+1/6*(k1h(1)+2*k2h(1)+2*k3h(1)+k4h(1))*stepX;
        T(i+1)=T(i)+1/6*(k1h(2)+2*k2h(2)+2*k3h(2)+k4h(2))*stepX;
        bplus(i)=Z(i)-h/k*(T(i)-Tu);
    end
    c=bplus(N);
    a=b1(N);
    df=(bplus(N)-b1(N))/eps;
    Znew=Zold-b1(N)/df;
    er=abs((Zold-Znew)/Znew)*100;
    Zold=Znew;
end

fprintf('Distribusi Suhu terhadap Jarak(x) Sepanjang Rod\n')
fprintf('Jarak(cm) \t\t Suhu Rod (K)\n')
fprintf('%5.2f \t\t\t %5.4f\n',[xfinal(1:10:101); Tfinal(1:10:101)])
plot(xfinal,Tfinal)

%Mencari heat lost
for j=1:N
    ql(j)=Tfinal(j)-Tu;
end
for l=2:N-1
    if (-1)^l>0
        ql(l)=4*ql(l);
    else
        ql(l)=2*ql(l);
    end
end
dz=(L-0)/(N-1);
q=(pi*D*h*dz/3*sum(ql))+(pi/4*D^2*h*(Tfinal(N)-Tu));
qideal=((pi*D*L)+(pi/4*D^2))*h*(Ts-Tu);
%Efisiensi Fin
ef=q/qideal*100;

%menulis hasil kalor hilang dan efisiensi fin
fprintf('\n')
fprintf('Kalor hilang pada Fin %5.4f cal/det \n',[q])
fprintf('Kalor hilang pada Fin (kondisi ideal) %5.4f cal/det \n',[qideal])
fprintf('Efisiensi Fin %4.2f persen \n',[ef])
fprintf('\n')

%daftar fungsi
    function y=frk(x,Z,T)
    dZdX=+4*h/(k*D) * ( T - Tu );
    dTdX=Z;
    y=[dZdX; dTdX];
    end
    function g=runge(A)
    g=A;
    end
    
end

