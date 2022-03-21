clear all; close all; clc

%constants
tau = 2.2;
D = 4.5;
a = 1;
tend = 6;

N = 100;
dx = D/(N-1);
% dt <= dx
DeltaT = dx*0.5;
Nt = tend/DeltaT;
Lambda = DeltaT/dx;
SpaceingTX=linspace(0,D,2*N+1);
% We make sure the CFL number satisfies, the stability condition
% I.e DeltaT\leq DeltaX



for i=0:N;
    u(1,i+1)=sin(2*pi*i*DeltaT/(tau));
    l(1,i+1)=sin(2*pi*i*DeltaT/(tau));
    c(1,i+1)=sin(2*pi*i*DeltaT/(tau));
end

for i=0:2*N
    u(i+1,1)=0;
    l(i+1,1)=0;
    c(i+1,1)=0;
end

for n=1:Nt
    
    for j=2:N
        u(n+1,j)= 1/2 *( u(n,j+1)+ u(n,j-1)) - a*Lambda/2*(u(n,j+1)-u(n,j-1));
        l(n+1,j)= l(n,j)-a*Lambda*(l(n,j)-l(n,j-1));
        c(n+1,j)= c(n,j)-a*Lambda/2*(c(n,j+1)-c(n,j-1))+a*a*Lambda*Lambda/2*(c(n,j+1)-2*c(n,j)+c(n,j-1));
    end
    %n
    %j
    n+1;
    
    
end
for i=1:2*N
    u(i,N+1)=2*u(i,N)-u(i,N-1);
    l(i,N+1)=2*l(i,N)-l(i,N-1);
    c(i,N+1)=2*c(i,N)-c(i,N-1);
end

vals=l(:,end)
%plot(SpaceingTX,vals)
mesh(l)

% for j=1:2
%     
%     for k=1:2
%         
%     end
%     
%     fh(k) = figure(j);
%     ax(k) = axes('parent',fh(j));
%     
%     hold(ax(k),'on')
%     hold off
%     lh(k) = figure(j);
%     ax(k) = axes('parent',fh(j));
%     
%     hold(ax(k),'on')
%     
% end



