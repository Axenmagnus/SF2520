clear all;
close all;
clc;
mu=1/82.45;
r0=[-mu,0];
r1=[1-mu,0];

N=[0,1
   -1,0];
r=[1.2, 0];
rprime=[0, -1];
%instansiating stuff
c=deriv(r,rprime,mu,r0,r1,N);

position= nan(2,100000);
position(:,1)=r;
derivative=nan(2,100000);
derivative(:,1)=rprime;
secondDerivative=[];

i=2;
h=0.01
while i<5
    secondDerivative(:,i)=deriv(position(:,i-1)',derivative(:,i-1)',mu,r0,r1,N);
    
    position(:,i)=position(:,i-1)+h*derivative(:,i-1);
    derivative(:,i)=derivative(:,i-1)+h*secondDerivative(:,i);
    i=i+1
    
end



function ODE=deriv(r,rprime,mu,r0,r1,N)

    ODE=-(1-mu)*(r-r0)/(norm(abs(r-r0)).^3)-mu*(r-r1)/(norm(abs(r-r1)).^3) + 2 * transpose(mtimes(N,rprime')) + r

end





