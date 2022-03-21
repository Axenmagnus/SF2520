clear all; close all; clc;

xbound=5;
ybound=2;

T0y=40;
T5y=400;
dtdy0=0;
dtdy2=0;

f=100;
HH=[0.1,0.05,0.025];

h=0.1;
N=(xbound/h)-1;
M=(ybound/h)+1;



T1 = diag(4*ones(1,N*M)) + diag(-1*ones(1,(N*M)-1),1) + diag(-1*ones(1,(N*M)-1),-1);


f=f*(ones(1,N*M))';
%BC=1/(h*h)*zeros(N*M);

n=1;

while n<N*M
    
    bc(n)=T0y;
    n=n+N-1;
    bc(n)=T5y;
    n=n+1;
    
    
end
bc=(bc/(h*h))';
F=f+bc;

I = -diag(ones(1,N*M-N));
Ttemp1=zeros(M*N,M*N);
Ttemp1(N+1:M*N,1:M*N-N)=I;
Ttemp2=zeros(M*N,M*N);
Ttemp2(1:M*N-N,N+1:M*N)=I;


T=Ttemp1+Ttemp2+T1;
n=1;
while n<N+1
    T(n,N+n)=-2;
    T(M*N-N+n,N*M-2*N+n)=-2;
    n=n+1;
    
end

n=N;

while n<N*M
    T(n,n+1)=0;
    T(n+1,n)=0;
    n=n+N;
end

C=T/(h*h)\F;
C=reshape(C,N,M);
Start=ones(M,1)*T0y';
End=ones(M,1)*T5y;
C=[Start';C;End'];

%mesh(C)
%zlabel('T(x,y)')

fcn = @(x,y) 6000* exp(-5*(x-1)*(x-1)-10*(y-1.5)*(y-1.5));

j=1;
l=1;
t=1;
NewF=zeros(N*M,1);
while l<M
    
    j=1;
   while j<N+1
       
       NewF(t)=fcn(j*h,l*h);
       
       t=t+1;
       j=j+1;
   end
   l=l+1;
end

NewF=NewF+bc;




n=0;

while n<N
    j=0;
    while j<M
        u(M*n+j+1)=40+72*xbound/M*j;
        j=j+1;
    end
    n=n+1;
    
    
end

cc=u
msss=reshape(u,M,N)
mesh(msss)


%% crank

AppendArray=[];
tau=10;
C=1;
dt=C*h^2;
t=0;
I=eye(N*M);
u=u';
A=-T*(1/h^2);
A=A;

un=u
answer=[];
%answer=[un,answer];
while t<tau
    %u=((I-1/2*dt*-A))\((I+1/2*dt*-A).*u'+1/2*dt*(2*NewF));
    
    %u=(sparse((I-0.5*dt.*A)))\( sparse(((I+1/2*dt.*A)*u) + (0.5*dt.*(2*NewF)) ));
    un1=(sparse((I-0.5*dt.*A)))\( sparse(((I+((1/2)*dt.*A))*un) + (0.5*dt.*(2*NewF)) ));
    answer=[answer,un1];
    un=un1;
    
    %AppendArray=[AppendArray,u];
    t=t+dt;
    
end
answer=[u,answer];
T=answer(:,100);

CC=reshape(T,N,M);
mesh(CC);
%Start=ones(M,1)*T0y';
%End=ones(M,1)*T5y;
%CC=[Start';CC;End'];

%mesh(CC)










