function T=TractionNGreenTensor_2D(omega,kp,ks,y,x)
% 向量形式积分
%y is point source
% 2017 12 18
% 积分路径是复平面上的

w2 = 1/omega^2 ;
n=size(x);
T=zeros(8,n(2));

L=5*(kp+ks);
lambda=1;

yy=y;
yy(2,:)=-y(2,:);
%%sigma.*e1
T1=TractionGreenTensor_2D_e1(omega,kp,ks,y,x);
Ts1=TractionGreenTensor_2D_e1(omega,kp,ks,yy,x);
%%sigma.*e2
T2=TractionGreenTensor_2D(omega,kp,ks,y,x);
Ts2=TractionGreenTensor_2D(omega,kp,ks,yy,x);
%%set value as zero when x=y
for i=1:n(2)
    if isnan([1 1 1 1]*T1(:,i))==1
        T1(:,i)=[0;0;0;0];
    end
end

for i=1:n(2)
    if isnan([1 1 1 1]*T2(:,i))==1
        T2(:,i)=[0;0;0;0];
    end
end

% 设置间接变量

mus=@(t) (ks^2-t.^2).^(1/2);
mup=@(t) (kp^2-t.^2).^(1/2);
beta=@(t) ks^2-2.*t.^2;
delta=@(t) beta(t).^2 + 4*t.^2.*mus(t).*mup(t);
alpha=@(t) ks^2-2*kp^2+2.*t.^2;
mu=omega^2/ks^2;
forec=@(t) -mu*w2./(2*pi*delta(t));

Ess=@(t) exp(1i*mus(t)*(x(2,:)+y(2,:))+1i*t*(x(1,:)-y(1,:)));
Epp=@(t) exp(1i*mup(t)*(x(2,:)+y(2,:))+1i*t*(x(1,:)-y(1,:)));
Esp=@(t) exp(1i*(mus(t)*x(2,:)+mup(t)*y(2,:))+1i*t*(x(1,:)-y(1,:)));
Eps=@(t) exp(1i*(mup(t)*x(2,:)+mus(t)*y(2,:))+1i*t*(x(1,:)-y(1,:)));

A1=@(t) forec(t).*((2.*t.*mus(t).*beta(t).^2).*Ess(t)+(4.*t.^3.*mus(t).*alpha(t)).*Epp(t)+        (4.*t.^3.*mus(t).*beta(t)).*Esp(t)+ (2.*t.*mus(t).*beta(t).*alpha(t)).*Eps(t));
A2=@(t) forec(t).*(              (beta(t).^3).*Ess(t)+      (8*mus(t).*mup(t).*t.^4).*Epp(t)+              (2*t.^2.*beta(t).^2).*Esp(t)+  (4*t.^2*mus(t).*mup(t).*beta(t)).*Eps(t));
A3=@(t) forec(t).*( (-8*mus(t).*mup(t).*t.^4).*Ess(t)+     (beta(t).^2.*alpha(t)).*Epp(t)+ (-4*t.^2*mus(t).*mup(t).*beta(t)).*Esp(t)+       (2*t.^2.*alpha(t).*beta(t)).*Eps(t));
A4=@(t) forec(t).*( (-4*t.^3*mup(t).*beta(t)).*Ess(t)+  (2*t.*mup(t).*beta(t).^2).*Epp(t)+        (-2*t.*mup(t).*beta(t).^2).*Esp(t)+          (4*t.^3*mup(t).*beta(t)).*Eps(t));

B2=@(t) forec(t).*(-(2.*t.*mus(t).*beta(t).^2).*Ess(t)+(4.*t.^3.*mus(t).*beta(t)).*Epp(t)-        (4.*t.^3.*mus(t).*beta(t)).*Esp(t)+ (2.*t.*mus(t).*beta(t).*beta(t)).*Eps(t));
B4=@(t) forec(t).*( (8*mus(t).*mup(t).*t.^4).*Ess(t)+     (beta(t).^3).*Epp(t)+ (4*t.^2*mus(t).*mup(t).*beta(t)).*Esp(t)+       (2*t.^2.*beta(t).*beta(t)).*Eps(t));
%integral part
Tc1=integral(@(t) A1(t) ,-L+lambda*1i,L-lambda*1i,'Waypoints',[-lambda+lambda*1i,lambda-lambda*1i], 'ArrayValued',true);
Tc2=integral(@(t) A2(t) ,-L+lambda*1i,L-lambda*1i,'Waypoints',[-lambda+lambda*1i,lambda-lambda*1i], 'ArrayValued',true);
Tc3=integral(@(t) A3(t) ,-L+lambda*1i,L-lambda*1i,'Waypoints',[-lambda+lambda*1i,lambda-lambda*1i], 'ArrayValued',true);
Tc4=integral(@(t) A4(t) ,-L+lambda*1i,L-lambda*1i,'Waypoints',[-lambda+lambda*1i,lambda-lambda*1i], 'ArrayValued',true);

Tc6=integral(@(t) B2(t) ,-L+lambda*1i,L-lambda*1i,'Waypoints',[-lambda+lambda*1i,lambda-lambda*1i], 'ArrayValued',true);
Tc8=integral(@(t) B4(t) ,-L+lambda*1i,L-lambda*1i,'Waypoints',[-lambda+lambda*1i,lambda-lambda*1i], 'ArrayValued',true);

T(1,:)=T1(1,:)-Ts1(1,:)+Tc1;
T(2,:)=T1(2,:)-Ts1(2,:)+Tc2;
T(3,:)=T1(3,:)-Ts1(3,:)+Tc3;
T(4,:)=T1(4,:)-Ts1(4,:)+Tc4;
%traction tensor is symmetry
T(5,:)=T(2,:);
T(6,:)=T2(2,:)-Ts2(2,:)+Tc6;
T(7,:)=T(4,:);
T(8,:)=T2(4,:)-Ts2(4,:)+Tc8;



    

