clear all;
close all;
clc

T=1;
N=[125,250,500,1000,2000];

r1=0.04;
r2=10*10*10*10;
r3=3*10.^7;

deriv1 = @(x1,x2,x3) -r1*x1+r2*x2*x3;
deriv2 = @(x1,x2,x3) r1*x1-r2*x2*x3-r2*x2*x2;
deriv3 = @(x1,x2,x3) r3*x2*x2;

%deriv= @(x1,x2,x3) [-r1x1


l=1;
endx1=[];
endx2=[];
endx3=[];
while l<length(N)+1
    
    x1=[1];
    x2=[0];
    x3=[0];
    n=1;
    
    while n<(N(l))
        h=T/N(l);


        k1=deriv1(x1(:,n),x2(:,n),x3(:,n));
        k2=deriv1(x1(:,n)+h*k1,x2(:,n),x3(:,n));
        k3=deriv1(x1(:,n)+h*k1/4+h*k2/4,x2(:,n),x3(:,n));

        x1temp=x1(:,n) + h/6 * (k1+k2+4*k3);

        k1=deriv2(x1(:,n),x2(:,n),x3(:,n));
        k2=deriv2(x1(:,n),x2(:,n)+h*k1,x3(:,n));
        k3=deriv2(x1(:,n),x2(:,n)+h*k1/4+h*k2/4,x3(:,n));

        x2temp=x2(:,n) + h/6 * (k1+k2+4*k3);

        k1=deriv3(x1(:,n),x2(:,n),x3(:,n));
        k2=deriv3(x1(:,n),x2(:,n),x3(:,n)+h*k1);
        k3=deriv3(x1(:,n),x2(:,n),x3(:,n)+h*k1/4+h*k2/4);

        x3temp=x3(:,n) + h/6 * (k1+k2+4*k3);

        x1(:,n+1)=x1temp;
        x2(:,n+1)=x2temp;
        x3(:,n+1)=x3temp;
        n=n+1;
    end
    
    endx1(l)=x1temp;
    endx2(l)=x2temp;
    endx3(l)=x3temp;
    l=l+1;
    x1temp
    x2temp
    x3temp

end
%endx1
%endx2
%endx3
%plot3(x1,x2,x3)
t=linspace(0,1,2000)
loglog(x1,t)
hold on
loglog(x2,t)
hold on 
loglog(x3,t)

