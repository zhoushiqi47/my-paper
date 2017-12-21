function T=TractionNGreenTensor_2D_e1_e2(omega,kp,ks,y,x)
% 向量形式积分
%y is point source
% 2017 12 20
% 积分路径是复平面上的
% y2=0 source on the gamma0

n=size(x);
T=zeros(8,n(2));

L=3*(kp+ks);
lambda=0.1;


% 设置间接变量

a= x(1,:)-y(1,:);
b= x(2,:);


%integral part
Tc1=integral(@(t) A1(t,a,b,kp,ks) ,-L+lambda*1i,L-lambda*1i,'Waypoints',[-lambda+lambda*1i,lambda-lambda*1i], 'ArrayValued',true);
Tc2=integral(@(t) A2(t,a,b,kp,ks) ,-L+lambda*1i,L-lambda*1i,'Waypoints',[-lambda+lambda*1i,lambda-lambda*1i], 'ArrayValued',true);
Tc3=integral(@(t) A3(t,a,b,kp,ks) ,-L+lambda*1i,L-lambda*1i,'Waypoints',[-lambda+lambda*1i,lambda-lambda*1i], 'ArrayValued',true);
Tc4=integral(@(t) A4(t,a,b,kp,ks) ,-L+lambda*1i,L-lambda*1i,'Waypoints',[-lambda+lambda*1i,lambda-lambda*1i], 'ArrayValued',true);

Tc6=integral(@(t) B2(t,a,b,kp,ks) ,-L+lambda*1i,L-lambda*1i,'Waypoints',[-lambda+lambda*1i,lambda-lambda*1i], 'ArrayValued',true);
Tc8=integral(@(t) B4(t,a,b,kp,ks) ,-L+lambda*1i,L-lambda*1i,'Waypoints',[-lambda+lambda*1i,lambda-lambda*1i], 'ArrayValued',true);

T(1,:)=Tc1;
T(2,:)=Tc2;
T(3,:)=Tc3;
T(4,:)=Tc4;
%traction tensor is symmetry
T(5,:)=T(2,:);
T(6,:)=Tc6;
T(7,:)=T(4,:);
T(8,:)=Tc8;
end

function f=A1(t,a,b,kp,ks)
mus=(ks^2-t^2)^(1/2);
mup=(kp^2-t^2)^(1/2);
beta=ks^2-2*t^2;
delta=beta^2+4*t^2*mus*mup;
alpha=ks^2-2*kp^2+2*t^2;
B=2*t*mus*beta;
C=2*t*mus*alpha;
Es=exp(1i*mus*b+1i*t*a);
Ep=exp(1i*mup*b+1i*t*a);
f=(B*Es+C*Ep)/(-2*pi*delta);
end

function f=A2(t,a,b,kp,ks)
mus=(ks^2-t^2)^(1/2);
mup=(kp^2-t^2)^(1/2);
beta=ks^2-2*t^2;
delta=beta^2+4*t^2*mus*mup;
alpha=ks^2-2*kp^2+2*t^2;
B=beta^2;
C=4*t^2*mus*mup;
Es=exp(1i*mus*b+1i*t*a);
Ep=exp(1i*mup*b+1i*t*a);
f=(B*Es+C*Ep)/(-2*pi*delta);
end

function f=A3(t,a,b,kp,ks)
mus=(ks^2-t^2)^(1/2);
mup=(kp^2-t^2)^(1/2);
beta=ks^2-2*t^2;
delta=beta^2+4*t^2*mus*mup;
alpha=ks^2-2*kp^2+2*t^2;
B=-4*t.^2.*mus*mup;
C=alpha*beta;
Es=exp(1i*mus*b+1i*t*a);
Ep=exp(1i*mup*b+1i*t*a);
f=(B*Es+C*Ep)/(-2*pi*delta);end

function f=A4(t,a,b,kp,ks)
mus=(ks^2-t^2)^(1/2);
mup=(kp^2-t^2)^(1/2);
beta=ks^2-2*t^2;
delta=beta^2+4*t^2*mus*mup;
alpha=ks^2-2*kp^2+2*t^2;
B=-2*t*mup*beta;
C=2*t*mup*beta;
Es=exp(1i*mus*b+1i*t*a);
Ep=exp(1i*mup*b+1i*t*a);
f=(B*Es+C*Ep)/(-2*pi*delta);
end

function f=B2(t,a,b,kp,ks)
mus=(ks^2-t^2)^(1/2);
mup=(kp^2-t^2)^(1/2);
beta=ks^2-2*t^2;
delta=beta^2+4*t^2*mus*mup;
alpha=ks^2-2*kp^2+2*t^2;
B=-2*t*mus*beta;
C=2*t*mus*beta;
Es=exp(1i*mus*b+1i*t*a);
Ep=exp(1i*mup*b+1i*t*a);
f=(B*Es+C*Ep)/(-2*pi*delta);
end

function f=B4(t,a,b,kp,ks)
mus=(ks^2-t^2)^(1/2);
mup=(kp^2-t^2)^(1/2);
beta=ks^2-2*t^2;
delta=beta^2+4*t^2*mus*mup;
alpha=ks^2-2*kp^2+2*t^2;
B=4*t^2*mus*mup ;
C=beta^2;
Es=exp(1i*mus*b+1i*t*a);
Ep=exp(1i*mup*b+1i*t*a);
f=(B*Es+C*Ep)/(-2*pi*delta);
end