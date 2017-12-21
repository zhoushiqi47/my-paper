%% ---------------- Parameters Setting---------------------------------------------%%

clear;
n_src = 301;
n_recv = 301;

%lamda=1;
%mu=1;
lamda=1/2;
mu=1/4;
%omega=[2,2.5,3,3.5,4,4.5,5,5.5,6,6.5,7,7.5,8]*pi;
omega=4*pi;
ks=omega/sqrt(mu);
kp=omega/sqrt(lamda+2*mu);
Nomegas = length(omega);
npts = 256 ; %% discretization point of integral equation

source   = zeros(2,n_src);
receiver = zeros(2,n_recv);

Q=[1 ,0;0 ,1];
a =50;

%% half space reflection
source(1,:)   = linspace(-a,a,n_src);     source(2,:)   =  zeros(1,n_src);
receiver(1,:) = linspace(-a,a,n_recv);    receiver(2,:) =  zeros(1,n_recv);
bctype1 = 1;  %% 1 circle;2 Penut;3 p_leaf 4 rectangle      
bctype2 = 1;


tic
%Data1 = zeros(n_recv,n_src,2,2,Nomegas);
 %% the obstacle is under the half space
 % 合成散射数据
 %Data = NyElasticWave_src_half_multi(omega , lamda, mu, npts, bctype1,bctype2, n_src, n_recv, source, receiver);
 %save scatteredfield_256_401_a50_r1_h10_kp8_ks16_Penut Data  
 %save scatteredfield_256_401_a50_r1_h10_multi_2_8_circle_5_times_5_peanut Data 
 load scatteredfield_256_401_a50_r1_h10_multi_2_8_bi_circle Data 
 %load scatteredfield_256_401_a100_r1_h10_kp2_ks4_Penut Data;
 %Data1(:,:,:,:,1)=Data;
 %load scatteredfield_256_401_a100_r1_h10_kp4_ks8_Penut Data;
 %Data1(:,:,:,:,2)=Data;
 Data1=Data(:,:,:,:,5);
 %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
 %% Add additive Gaussian random noise.
%     data(:,:,k) = U(:,:,k);
%     maxU = max(max(abs( U(:,:,k))));
%     noiselevel = 0.2*maxU;
 %  U(:,:,k) =  U(:,:,k) + noiselevel/sqrt(2)*maxU*((2*rand(size( U(:,:,k)))-1)+i*(2*rand(size( U(:,:,k)))-1));
level=0;
data=zeros(size(Data1));
for kk=1:Nomegas
    for j=1:2
        for l=1:2
            data(:,:,j,l,kk)=Data1(:,:,j,l,kk);
            maxU=max(max(abs(Data1(:,:,j,l,kk))));
            noiselevel = level*maxU;
            %data(:,:,j,l,kk)=data(:,:,j,l,kk)+noiselevel/sqrt(2)*((2*rand(n_recv,n_src)-1)+1i*(2*rand(n_recv,n_src)-1));
            data(:,:,j,l,kk)=data(:,:,j,l,kk)+1/sqrt(2)*(normrnd(0,noiselevel,n_recv,n_src)+1i*normrnd(0,noiselevel,n_recv,n_src));
        end
    end
end





%% Sampling domain;
Nx = 401;
Nz = 201;
x = linspace( -4, 4, Nx);
z = linspace( 8, 12, Nz);

%% 逆时偏移成像
Irtm = RTM_halfspace_elastic2(n_src, n_recv,source, receiver, x, z , data, omega , kp ,ks);
%load circle_r1_401_-10_256_201 Irtm


[xx,zz]=meshgrid(x,z);
 
figure, imagesc(x,z,imag(Irtm)');colorbar; 
colormap(jet)

%% Plot the true scatter
n = npts;
node = 0:2*n-1;
t = pi*node(:)/n;
if bctype1==1
    [xt,zt]=circlebcl(t,1);   
else if bctype1==2
    [xt,zt]=Penutl(t,1);   
    else if bctype1==3
            [xt,zt]=p_leafl(t,1);   
        else
            [xt,zt]=myrectanglel(t,1);   
        end
    end
end
hold on;
plot(xt,zt,'r');

if bctype2==1
    [xt,zt]=circlebcr(t,1);   
else if bctype2==2
    [xt,zt]=Penutr(t,1);   
    else if bctype2==3
            [xt,zt]=p_leafr(t,1);   
        else if bctype2==4
            [xt,zt]=myrectangler(t,1);   
            else
                [xt,zt]=kiter(t,1);
            end
        end
    end
end
hold on;
plot(xt,zt,'r');

%figure, imagesc(x,z,real(Irtm)');colorbar; 
%colormap(jet)

%hold on;
%plot(xt,zt,'r');

