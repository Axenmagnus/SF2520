clear all; close all; clc

%constants
a=1;
tau=2.2;
D=4.5;
T=6;
N=100;
SpaceingT=linspace(0,T,N/2);
DeltaT=T/(N*2);
SpaceingX=linspace(0,D,N);
DeltaX=(D/(N));
Lambda= DeltaT/DeltaX;
% We make sure the CFL number satisfies, the stability condition
% I.e DeltaT\leq DeltaX


for i=1:N*2+1;
    u(1,i)=sin(2*pi*i*DeltaT/(tau));
end

for i=1:N
    u(i,1)=0;
end



for i=1:N -1
    
    for k=2:N*2
        u(i+1,k)=1/2*(u(i,k-1)+u(i,k+1))-a*Lambda/2*(u(i,k-1)-u(i,k+1));
        
    end
    u(i,N*2)=2*u(i,N*2-1)-u(i,N*2-2);
end
mesh(u)
vals=u(:,end)
plot(SpaceingX,vals)
