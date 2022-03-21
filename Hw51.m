clear all; close all; clc


%% Plot relative
%Constants picked arbitrarily
n=40;
d=3;
error=1e-4;
N=n^d;
b=rand(N,1);

figure(1)
[resj,sol2]=Jacobi(n,d,error,b);
xJ=linspace(0,1,length(resj));
semilogy(xJ,resj)
hold on

[resC,sol]=Conjugate(n,d,error,b);
xC=linspace(0,1,length(resC));
semilogy(xC,resC)
hold on
legend('Jacobi','Conjugate')
title("Relative residual")
hold off

norm(sol-sol2)

%% Plot absolute
n=30;
d=3;
error=10e-14;
N=n^d;
b=rand(N,1);
figure(1)
[resj,sol2]=Jacobi(n,d,error,b);
xJ=linspace(0,length(resj),length(resj));
semilogy(xJ,resj)
hold on

[resC,sol]=Conjugate(n,d,error,b);
xC=linspace(0,length(resC),length(resC));
semilogy(xC,resC)
hold on
legend('Jacobi','Conjugate')
title("Iteration curve")
hold off
%% 
%Idenfying constants 
val1=resj(50);

L=resj<val1*0.1;
index=find(L,1,'first')
val2=L(index);
CJ=val1/val2*1/(index-50)

val1=resC(50);

L=resC<val1*0.1;
index=find(L,1,'first')
val2=L(index);
CC=val1/val2*1/(index-50)



%% Comparing various speeds with n
clear all;close all;clc
%Also used for determining which method converges faster
%Jacobi
%constants are set arbitrarily
d=2;
epsilon=0.00001;
ItersJ=[];
spacingN=linspace(5,50,50-4);
for n=5:1:50
    J=Jacobi(n,d,epsilon);
    ItersJ=[ItersJ length(J)];
    
    
end

% Conjugate
ItersC=[];
for n=5:1:50
    [C,b,sol]=Conjugate(n,d,epsilon);
    ItersC=[ItersC length(C)];
    
    
end

plot(spacingN,ItersC)
hold on 
plot(spacingN,ItersJ)



%% Comparing various speeds with d, with fixed N 
clear all;close all;clc
%Also used for determining which method converges faster
%Jacobi
%constants are set arbitrarily
N=10000;
epsilon=0.0001;
ItersJ=[];
spacingd=linspace(3,7,7-3+1);
for d=3:7
    n=round(nthroot(N,d));
    J=Jacobi(n,d,epsilon);
    ItersJ=[ItersJ length(J)];
    
    
end

% Conjugate
ItersC=[];
for d=3:7
    n=round(nthroot(N,d));
    [C,b,sol]=Conjugate(n,d,epsilon);
    ItersC=[ItersC length(C)];
    
    
end

plot(spacingd,ItersC)
hold on 
plot(spacingd,ItersJ)


%% b part
clear all; close all; clc

%Varying d
n=5;
error=1e-10;
for d=1:6
    N=n^d;
    b=rand(N,1);
    [resC,sol]=Conjugate(n,d,error,b);
    A=lap(n,d);
    L=["Time for backslash"];
    disp(L)
    tic
    sol1=A\b;
    toc
    norm(sol1-sol)
end
    


%% b part
clear all; close all; clc

%Varying n
d=3;
error=1e-10;
for n=2:5:27
    n
    N=n^d
    b=rand(N,1);
    [resC,sol]=Conjugate(n,d,error,b);
    A=lap(n,d);
    L=["Time for backslash"];
    disp(L)
    tic
    sol1=A\b;
    toc
    norm(sol1-sol)
end






