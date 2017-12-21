clear;


n_src = 1 ;
n_recv = 10;

flag = 1;
theta_r = (0:n_recv-1)*2*pi/n_recv;
theta_s = (0:n_src-1)*2*pi/n_src ;

n = 256; %% discretization point of integral equation
pp = 1;

lamda=1/2;

mu=1/4;
omega=pi/2;
ks=omega/sqrt(mu);
kp=omega/sqrt(lamda+2*mu);

% for R=[10 8 6 4 3.5 3.2 3.1 2.99 2.49];
%source   = zeros(2,n_src);
receiver = zeros(2,n_recv);
R = 20;
a = 0;
% source(1,:)   = R*cos(theta_s); source(2,:)   = R*sin(theta_s);

%receiver(1,:) = R*cos(theta_r); receiver(2,:) = R*sin(theta_r)-10;
receiver(1,:) = linspace(0,R,n_recv);    receiver(2,:) =  a*ones(1,n_recv);
%source=[0;10];
source=[0;0];
receiver;
Q=[1 ,0;0 ,1];
bctype = 1;
bctype2 = 2;
gamma=0.25;
%% solve via Nystrom's method

%% U =   NystromImpedance(npts, n_src, n_recv, bctype, freq, source, receiver);
%Data = NyElasticWave_src_half_multi(omega, lamda, mu, n, bctype1,bctype2, n_src, n_recv, source, receiver);
%[Data] = NyElasticWave_src_half_neumann_lite(omega, lamda, mu, n, bctype, n_src, n_recv, source, receiver);
 Data =NyElasticWave_src_half_transmission(omega, lamda, mu, n, bctype, n_src, n_recv, source, receiver,gamma);
Data2= Elastic_GreenTensor_Thalf_SIP2(omega,kp,ks,source(:,1)*ones(1,n_recv),receiver);

%Data3= Elastic_GreenTensor_Thalf_SIP2(omega,kp,ks,source(:,2)*ones(1,n_recv),receiver);

[Data(:,:,1,1) -Data2(1,:).']
