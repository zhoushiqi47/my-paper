function Data = NyElasticWave_src_half_impedance_lite(omegas, lambda, mu, n, bctype, n_src, n_recv, source, receiver,gamma)
%% the exact solution is given by the radiation solution
%% 2017 12 21 by zhou
%% gamma is impedance constanct





%% discrete point
node = 0:2*n-1;
t = pi*node(:)/n;
if bctype==1
    [x1,x2]=circlebc(t,1);   
    [dx1,dx2]=circlebc(t,2);
    [ddx1,ddx2]=circlebc(t,3);

else if bctype==2
    [x1,x2]=Penut(t,1);   
    [dx1,dx2]=Penut(t,2);
    [ddx1,ddx2]=Penut(t,3);
    else if bctype==3
            [x1,x2]=p_leaf(t,1);   
            [dx1,dx2]=p_leaf(t,2);
            [ddx1,ddx2]=p_leaf(t,3);
        else
            [x1,x2]=myrectangle(t,1);   
            [dx1,dx2]=myrectangle(t,2);
            [ddx1,ddx2]=myrectangle(t,3);
        end
    end
end

distance = sqrt( dx1.*dx1+dx2.*dx2 );

Nomegas = length(omegas);
Data = zeros(n_recv,n_src,2,2,Nomegas);
%% 对多频进行循环
for kk=1:Nomegas
    omega = omegas(kk);
    gamma = omega*gamma;
    cp = sqrt( lambda + 2*mu);
    cs = sqrt( mu );

    %% the corresponding wavenumber;
    kp = omega / cp;
    ks = omega / cs;
    tic
    A1 = singlelayar(omega, lambda, mu, n, bctype);
    A2 = singlelayar_traction(omega, lambda, mu, n, bctype)-1/2*diag([distance;distance]);
    toc
    disp('参数矩阵生成');
    %% 构造线性系统，方程左端的矩阵 %%
    

    A = 1i*gamma*diag([distance;distance])*A1+A2;
    tic
    %% 计算线性方程右端项
    f=zeros(4*n,2*n_src);
    % 同样的将循环去掉，用长向量代替
    xsource = reshape(repmat(source(:,:),2*n,1),2,2*n*n_src);
    xobs = repmat([x1,x2]',1,n_src);
    
    %Gs1 = Elastic_GreenTensor_Thalf_SIP4(omega,kp,ks,xsource,xobs);
    Gs1 = Elastic_GreenTensor_Thalf_SIP4(omega,kp,ks,xobs,xsource);
    Gs = - transpose(Gs1); 
    Ts1 = TractionNGreenTensor_2D_e1_e2(omega,kp,ks,xsource,xobs);
    Ts = - transpose(Ts1); 
    nx1= repmat(dx1,n_src,1);
    nx2= repmat(dx2,n_src,1);
    T1 = Ts(:,1).*nx2 - Ts(:,5).*nx1;
    T2 = Ts(:,2).*nx2 - Ts(:,6).*nx1;
    T3 = Ts(:,3).*nx2 - Ts(:,7).*nx1;
    T4 = Ts(:,4).*nx2 - Ts(:,8).*nx1;
    % 右端项为点源两个方向的组合，每一列代表一个方向
    f(1:2*n,:) = [reshape(T1,2*n ,n_src)+1i*gamma*diag(distance)*reshape(Gs(:,1),2*n ,n_src),reshape(T3,2*n ,n_src)+1i*gamma*diag(distance)*reshape(Gs(:,2),2*n ,n_src)];
    f(2*n+1:end,:) = [reshape(T2,2*n ,n_src)+1i*gamma*diag(distance)*reshape(Gs(:,3),2*n ,n_src),reshape(T4,2*n ,n_src)+1i*gamma*diag(distance)*reshape(Gs(:,4),2*n ,n_src)];
    toc
    disp('右端项计算完成');
    
    
    
    %% 解由单层位势表示的线性系统  
    %% U = S\phi
    
    
    % Phi是包含两个方向点源产生的核，前n_src代表一个方向，其余代表一个方向
    tic
    Phi = A\f;
    
    toc
    disp('得到线性方程解，即位势核');
    
   %% 计算样本点的散射数据
    % 同样的将循环去掉，用长向量代替
    tic
    if n_src==n_recv
        
        Gr = Gs1;
 
    else
        xreceiver = reshape(repmat(receiver(:,:),2*n,1),2,2*n*n_recv);
        xobs2= repmat([x1,x2]',1,n_recv);
        Gr = Elastic_GreenTensor_Thalf_SIP4(omega,kp,ks,xobs2,xreceiver);
    end
    
    G1=reshape(Gr(1,:),2*n,n_recv).';
    G2=reshape(Gr(2,:),2*n,n_recv).';
    G3=reshape(Gr(3,:),2*n,n_recv).';
    G4=reshape(Gr(4,:),2*n,n_recv).';
    distances = repmat(distance(:),1,2*n_src);
    
    toc
    disp('单层位势函数（障碍物到接收点）计算完成');
    
    tic
    
    temp1 = G1*(distances.*Phi(1:2*n,:)) + G3*(distances.*Phi(2*n+1:end,:));
    temp2 = G2*(distances.*Phi(1:2*n,:)) + G4*(distances.*Phi(2*n+1:end,:));
    
    
    
        Data(:,:,1,1,kk)=pi/n*temp1(:,1:n_src);
        Data(:,:,2,1,kk)=pi/n*temp2(:,1:n_src);
        Data(:,:,1,2,kk)=pi/n*temp1(:,(n_src+1):2*n_src);
        Data(:,:,2,2,kk)=pi/n*temp2(:,(n_src+1):2*n_src);
    
    toc
end

disp('数据合成结束');
end