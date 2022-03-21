%Labb 1 applied Numa
clear all; close all; clc;

N=10000; %steps
T=40; %end time
h=T/N;
alpha=0.1;
a=0.5*[1,1,sqrt(2)];



c=a*alpha;
m=[1,0,0];
m2=[1,0,0];

A = nan(3,N);
absA=nan(0);
i=0;
%Test %steps
%T=60; %end time
%h=2.14;
%N=T/h
%N=round(N);
%A = nan(3,N);

for x = 0:h:T
    i=i+1;
    k1=deriv(a,m,c);
    k2=deriv(a,m+h*k1,c);
    k3=deriv(a,m+h*k1/4+h*k2/4,c);
    
    m=m+h/6*(k1+k2+4*k3);
    A(:,i)=m;
    
    
end
%plot(0:h:T,A)
%title("m vector with time")
%xlabel('time') 
%ylabel('Value of components of m vector') %

%saveas(gcf,'spiral2.png')
absA=vecnorm(A);
%M=A/norm(A);
%A=A/vecnorm(A);

%plot3(1/2,1/2,0.5*sqrt(2))

k=1;
p=vecnorm(A);
while k<length(A)-1
    A(:,k)=A(:,k)/p(:,k);
    k=k+1;
    
   
    
end


plot3([A(1:1,:) 0],[A(2:2,:) 0],[A(3:3,:) 0]);
hold on
plot3([0 1/2],[0 1/2],[0 sqrt(2)/2] ,'Color','r');
title("m/|m| with t=40")
%ylabel('Value of components of m vector') 
legend( 'm with time','Initial vector a')
axis equal
saveas(gcf,'spiral3d1.png')


%Starting next bullet point
N=[40,80,160,320,640];
N2=[80,160,320,640,1280];
T=40*ones(1,6);
h1=T/vecnorm(N);
h2=T/vecnorm(N2);
k=1
while k<length(N)+1
    h1(:,k)=T(:,k)/N(:,k);
    h2(:,k)=T(:,k)/N2(:,k);
    k=k+1;
    
   
    
end

i=0;
r=0;
Diff = nan(1,length(h1));
for p=1:length(h1)
    
    for x = 0:h:T
        i=i+1;
        k1=deriv(a,m,c);
        k2=deriv(a,m+h1(1,p)*k1,c);
        k3=deriv(a,m+h1(1,p)*k1/4+h1(1,p)*k2/4,c);

        m=m+h/6*(k1+k2+4*k3);

        %c=c+1;
        k4=deriv(a,m2,c);
        k5=deriv(a,m2+h2(1,p)*k1,c);
        k6=deriv(a,m2+h2(1,p)*k1/4+h2(1,p)*k2/4,c);

        m2=m2+h/6*(k4+k5+4*k6);


    end
    Diff(1,p)=vecnorm(m2-m);
end



%loglog(Diff,h1)
%hold on
% loglog(h1.*h1.*h1.*h1.*h1,h1)
% legend('Difference', 'h^5')
%saveas(gcf,'Loglog.png')
%Solved analytically (by hand)
AMatrix = 1/2*[-0.1*3/2 -sqrt(2)+0.1/2 1+sqrt(2)/2*0.1;
    sqrt(2)+0.05 -0.1*3/2 -1+0.1*sqrt(2)/2;
    -1+sqrt(2)/2*0.1 1+sqrt(2)/2*0.1 -0.1];

eigenvalues=eig(AMatrix);

eiga=eigenvalues(1:1);
h=0.1;
z=h*eiga;
criteria = norm(1+z+z.^2/2+z.^3/6);

while criteria<1
    
    h=h+0.001;
    z=h*eiga;
    criteria = norm(1+z+z.^2/2+z.^3/6);

    
    
    
    
    
end














function f=deriv(a,m,c)
    
    f= cross(a,m) + cross(c,cross(a,m));
    
end



