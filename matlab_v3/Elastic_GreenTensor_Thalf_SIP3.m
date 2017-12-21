function G=Elastic_GreenTensor_Thalf_SIP3(omega,kp,ks,y,x)
%% x2=0
%% Traction free wave Green Tensor
% the integral path is SIP
% 向量形式积分
%y is point source
% 2016 12 17 
%modify 2017 11 28
n=size(x);
G=zeros(4,n(2));
L=3*(kp+ks);
mu = omega.^2./(ks.^2);
lambda=0.01;

Gc1=integral(@(t) A1(t,x(:,:),y(:,:),kp,ks,mu) ,-L+lambda*1i,L-lambda*1i,'Waypoints',[-lambda+lambda*1i,lambda-lambda*1i], 'ArrayValued',true);
Gc2=integral(@(t) A2(t,x(:,:),y(:,:),kp,ks,mu) ,-L+lambda*1i,L-lambda*1i,'Waypoints',[-lambda+lambda*1i,lambda-lambda*1i], 'ArrayValued',true);
Gc3=integral(@(t) A3(t,x(:,:),y(:,:),kp,ks,mu) ,-L+lambda*1i,L-lambda*1i,'Waypoints',[-lambda+lambda*1i,lambda-lambda*1i], 'ArrayValued',true);
Gc4=integral(@(t) A4(t,x(:,:),y(:,:),kp,ks,mu) ,-L+lambda*1i,L-lambda*1i,'Waypoints',[-lambda+lambda*1i,lambda-lambda*1i], 'ArrayValued',true);

G(1,:)=Gc1;
G(2,:)=Gc2;
G(3,:)=Gc3;
G(4,:)=Gc4;
end


function f=A1(t,x,y,kp,ks,mu)
mus=(ks^2-t^2)^(1/2);
mup=(kp^2-t^2)^(1/2);
beta=ks^2-2*t^2;
delta=beta^2+4*t^2*mus*mup;
B=mus*beta;
C=2*mus*t^2;
Es=exp(1i*mus*(y(2,:))+1i*t*(x(1,:)-y(1,:)));
Ep=exp(1i*mup*(y(2,:))+1i*t*(x(1,:)-y(1,:)));
f=1i*(B*Es+C*Ep)/(2*pi*delta*mu);
end

function f=A2(t,x,y,kp,ks,mu)
mus=(ks^2-t^2)^(1/2);
mup=(kp^2-t^2)^(1/2);
beta=ks^2-2*t^2;
delta=beta^2+4*t^2*mus*mup;
B=2*t*mus*mup;
C=-t*beta;
Es=exp(1i*mus*(y(2,:))+1i*t*(x(1,:)-y(1,:)));
Ep=exp(1i*mup*(y(2,:))+1i*t*(x(1,:)-y(1,:)));
f=1i*(B*Es+C*Ep)/(2*pi*delta*mu);
end

function f=A3(t,x,y,kp,ks,mu)
mus=(ks^2-t^2)^(1/2);
mup=(kp^2-t^2)^(1/2);
beta=ks^2-2*t^2;
delta=beta^2+4*t^2*mus*mup;
C=-2*t*mus*mup;
B=t*beta;
Es=exp(1i*mus*(y(2,:))+1i*t*(x(1,:)-y(1,:)));
Ep=exp(1i*mup*(y(2,:))+1i*t*(x(1,:)-y(1,:)));
f=1i*(B*Es+C*Ep)/(2*pi*delta*mu);
end

function f=A4(t,x,y,kp,ks,mu)
mus=(ks^2-t^2)^(1/2);
mup=(kp^2-t^2)^(1/2);
beta=ks^2-2*t^2;
delta=beta^2+4*t^2*mus*mup;
B=2*t^2*mup;
C=mup*beta;
Es=exp(1i*mus*(y(2,:))+1i*t*(x(1,:)-y(1,:)));
Ep=exp(1i*mup*(y(2,:))+1i*t*(x(1,:)-y(1,:)));
f=1i*(B*Es+C*Ep)/(2*pi*delta*mu);
end