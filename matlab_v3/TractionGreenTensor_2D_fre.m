function T=TractionGreenTensor_2D_fre(omega,kp,ks,y,x)
% 向量形式积分
%y is point source
% 2017 12 18
% 积分路径是复平面上的

w2 = 1/omega^2 ;
n=size(x);
T=zeros(8,n(2));

L=5*(kp+ks);
lambda=1;


ss=sign(x(2,:)-y(2,:));
% 设置间接变量

d1= x(1,:)-y(1,:);
d2=abs(x(2,:)-y(2,:));

mus=@(t) (ks^2-t.^2).^(1/2);
mup=@(t) (kp^2-t.^2).^(1/2);
beta=@(t) ks^2-2.*t.^2;
delta=@(t) beta(t).^2 + 4*t.^2.*mus(t).*mup(t);
alpha=@(t) ks^2-2*kp^2+2.*t.^2;
mu=omega^2/ks^2;
forec=@(t) mu*w2./(4*pi);

Ess=@(t) exp(1i*mus(t)*d2+1i*t*d1);
Epp=@(t) exp(1i*mup(t)*d2+1i*t*d1);


A1=@(t) forec(t).*(-2.*t.*mus(t).*Ess(t)+(-alpha(t).*t./mup(t)).*Epp(t));
A2=@(t) forec(t).*(  -beta(t).*Ess(t)-2*t.^2.*Epp(t));
A3=@(t) forec(t).*( 2*t.^2.*Ess(t)-alpha(t).*Epp(t));
A4=@(t) forec(t).*( t.*beta(t)./mus(t).*Ess(t)-2*t.*mup(t).*Epp(t));

B2=@(t) forec(t).*(2*t.*mus(t).*Ess(t)-t.*beta(t)./mup(t).*Epp(t));
B4=@(t) forec(t).*( -2*t.^2.*Ess(t)-beta(t).*Epp(t));
%integral part
Tc1=integral(@(t) A1(t) ,-L+lambda*1i,L-lambda*1i,'Waypoints',[-lambda+lambda*1i,lambda-lambda*1i], 'ArrayValued',true);
Tc2=integral(@(t) A2(t) ,-L+lambda*1i,L-lambda*1i,'Waypoints',[-lambda+lambda*1i,lambda-lambda*1i], 'ArrayValued',true);
Tc3=integral(@(t) A3(t) ,-L+lambda*1i,L-lambda*1i,'Waypoints',[-lambda+lambda*1i,lambda-lambda*1i], 'ArrayValued',true);
Tc4=integral(@(t) A4(t) ,-L+lambda*1i,L-lambda*1i,'Waypoints',[-lambda+lambda*1i,lambda-lambda*1i], 'ArrayValued',true);

Tc6=integral(@(t) B2(t) ,-L+lambda*1i,L-lambda*1i,'Waypoints',[-lambda+lambda*1i,lambda-lambda*1i], 'ArrayValued',true);
Tc8=integral(@(t) B4(t) ,-L+lambda*1i,L-lambda*1i,'Waypoints',[-lambda+lambda*1i,lambda-lambda*1i], 'ArrayValued',true);

T(1,:)=Tc1;
T(2,:)=Tc2.*ss;
T(3,:)=Tc3.*ss;
T(4,:)=Tc4;
%traction tensor is symmetry
T(5,:)=T(2,:);
T(6,:)=Tc6;
T(7,:)=T(4,:);
T(8,:)=Tc8.*ss;



    

