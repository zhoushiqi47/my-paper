function Data = NyElasticWave_src_half_impedance(omegas, lambda, mu, n, bctype, n_src, n_recv, source, receiver,gamma)
%% the exact solution is given by the radiation solution
%% 2017 12 21 by zhou
%% gamma is impedance constanct



w = quad_weights(n);
R = zeros(2*n);
lamda=lambda;
for k=1:2*n
     idx=k:2*n;
     R(idx,k)=w(1:2*n-k+1);
     R(k,k)=R(k,k)/2;  %% for convinience
end
R=(R+R');

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


Ceuler = 0.577215665;
distance = sqrt( dx1.*dx1+dx2.*dx2 );

Nomegas = length(omegas);
Data = zeros(n_recv,n_src,2,2,Nomegas);
%% 对多频进行循环
for kk=1:Nomegas
    omega = omegas(kk);
    gamma = omega*gamma;
    cp = sqrt( lambda + 2*mu);
    cs = sqrt( mu );
    cp2 = lambda + 2*mu;
    %% the corresponding wavenumber;
    kp = omega / cp;
    ks = omega / cs;

    %% constant definition
    w2 = omega^2;
    w_2 = 1/omega^2 ;
    CC = mu/(lambda+2*mu)/(2*pi);
    
    alpha = 1/(2*pi)*( - 1/(4*w2)*(ks^2 + kp^2) );
    % beta1 = 1/(2*pi)*( 1/(32*w2)*( 3*ks^4 + kp^4 ) );
    % beta2 = 1/(2*pi)*( 1/(16*w2)*(   kp^4 - ks^4 ) );
    ka1 = -1/(4*pi*w2)*( ks^2*log(ks/2) + kp^2*log(kp/2) + 1/2*(ks^2 - kp^2) + (Ceuler - 1i*pi/2)*(ks^2+kp^2)  );
    ka2 = 1/(4*pi*w2)*( ks^2 - kp^2 );
    %% end of constant definiation

    M1 = cell(2,2);
    H1 = cell(2,2);
    M2 = cell(2,2);
    H2 = cell(2,2);
    for k=1:4
        M1{k}=zeros(2*n);
        H1{k}=zeros(2*n);
        M2{k}=zeros(2*n);
        H2{k}=zeros(2*n);
    end
    tic
    %% 计算矩阵green tensor
    % 将矩阵拉开成向量
    xl = repmat([x1,x2]',1,2*n);
    yl = reshape( repmat([x1,x2]',2*n,1),2,4*n*n);
    rGG  = Elastic_GreenTensor_Thalf_SIP_NY(omega,kp,ks,yl,xl);
    rT  = TractionNGreenTensor_2D(omega,kp,ks,yl,xl);


    
    
    
    
    for j=1:2*n
        dx1j=dx1(j); dx2j=dx2(j);
        for k=1:2*n
            if (j==k)
                dist = distance(k);
                
                G1 = rGG(:,(k-1)*2*n+j);
                G = reshape(G1,2,2); 

                T = rT(:,(k-1)*2*n+j);
                T1 = reshape(T(1:4),2,2); 
                T2 = reshape(T(5:8),2,2);
                Tn = (T1*dx2j-T2*dx1j)*dist;
                M1{1,1}(j,k) = alpha; M1{2,2}(j,k) = alpha;
                M1{1,2}(j,k) = 0; M1{2,1}(j,k) = 0;
            
                H1{1,1}(j,k) = 2*alpha*log(dist) + ka1 + ka2*dx1(j)*dx1(j)/dist^2+G(1,1);
                H1{2,2}(j,k) = 2*alpha*log(dist) + ka1 + ka2*dx2(j)*dx2(j)/dist^2+G(2,2);
            
                H1{1,2}(j,k) =                           ka2*dx1(j)*dx2(j)/dist^2+G(1,2);
                H1{2,1}(j,k) =                           ka2*dx2(j)*dx1(j)/dist^2+G(2,1);
                for l=1:2
                    for m=1:2
                        M1{l,m}(j,k) = dist*M1{l,m}(j,k);
                        H1{l,m}(j,k) = dist*H1{l,m}(j,k);
                    end
                end
                
                beta = 1/2*(ddx1(j)*dx2j-ddx2(j)*dx1j)/dist;
            
                H2{1,1}(j,k) = Tn(1,1)+(2*(lambda+mu)*dx1j^2/dist^2/mu+1)*beta*CC;
                H2{2,2}(j,k) = Tn(2,2)+(2*(lambda+mu)*dx2j^2/dist^2/mu+1)*beta*CC;
            
                H2{1,2}(j,k) = Tn(1,2)+(2*(lambda+mu)*dx1j*dx2j/dist^2/mu)*beta*CC-1/2*(ddx1(j)*dx1j+ddx2(j)*dx2j)/dist*CC;
                H2{2,1}(j,k) = Tn(2,1)+(2*(lambda+mu)*dx1j*dx2j/dist^2/mu)*beta*CC+1/2*(ddx1(j)*dx1j+ddx2(j)*dx2j)/dist*CC;
            else
                dist = distance(k);
                v = ([x1(j)-x1(k) x2(j)-x2(k)]); 
                z1 = v(1);
                z2 = v(2);
                lg4s = log(4*sin((t(j)-t(k))/2)^2) ; 
                rv = norm(v);
                G1 = rGG(:,(k-1)*2*n+j);
                G = reshape(G1,2,2);
                
                T = rT(:,(k-1)*2*n+j);
                T1 = reshape(T(1:4),2,2); 
                T2 = reshape(T(5:8),2,2);
                Tn = (T1*dx2j-T2*dx1j)*dist;
            
                
                Hp0 = besselj(0,kp*rv);
                Hp1 = besselj(1,kp*rv);

                Hs0 = besselj(0,ks*rv);
                Hs1 = besselj(1,ks*rv);
                
                phi1 = 1/(2*pi)*( -1/(2*mu)*Hs0 + 1/(2*w2*rv)* ( ks*Hs1 - kp*Hp1) ) ;
                phi2 = 1/(2*pi)*( 1/(2*w2)*( ks^2*Hs0 -2*ks/rv*Hs1 - kp^2*Hp0 + 2*kp/rv*Hp1  ) );
                
                r2 = rv*rv;
                r3= r2*rv;

                phis1 = -1/(4*pi)*w_2*(-ks^2*Hs0 + 2*ks*Hs1./rv);
                phip1 = -1/(4*pi)*w_2*(-kp^2*Hp0 + 2*kp*Hp1./rv);

                phis2 = -1/(4*pi)*w_2*(ks^3*Hs1 + 4*ks^2*Hs0./rv - 8*ks*Hs1./r2);
                phip2 = -1/(4*pi)*w_2*(kp^3*Hp1 + 4*kp^2*Hp0./rv - 8*kp*Hp1./r2);

                gs111 = 3*phis1.*z1./r2 + phis2.*z1.^3./r3;
                gs112 = phis1.*z2./r2 + phis2.*z1.^2.*z2./r3;
                gs122 = phis1.*z1./r2 + phis2.*z2.^2.*z1./r3;
                gs222 = 3*phis1.*z2./r2 + phis2.*z2.^3./r3;

                gp111 = 3*phip1.*z1./r2 + phip2.*z1.^3./r3;
                gp112 = phip1.*z2./r2 + phip2.*z1.^2.*z2./r3;
                gp122 = phip1.*z1./r2 + phip2.*z2.^2.*z1./r3;
                gp222 = 3*phip1.*z2./r2 + phip2.*z2.^3./r3;
                
                M2{1,1}(j,k)=(-2*mu*gs122 - cp2*gp111 - lamda*gp122)*dx2j-(mu*(gs112-gs222) - 2*mu*gp112)*dx1j;
                M2{2,1}(j,k)=(mu*(gs112-gs222) - 2*mu*gp112)*dx2j-(2*mu*gs122 - cp2*gp122 - lamda*gp111)*dx1j;
                M2{1,2}(j,k)=(2*mu*gs112 - cp2*gp112 - lamda*gp222)*dx2j-( mu*(gs122-gs111) - 2*mu*gp122)*dx1j;
                M2{2,2}(j,k)=(mu*(gs122-gs111) - 2*mu*gp122)*dx2j-( -2*mu*gs112 - cp2*gp222 - lamda*gp112)*dx1j;
                for l=1:2
                    for m=1:2
                        e=0;
                        if l==m
                            e=1;
                        end
                    
                        M1{l,m}(j,k) = dist*(phi1*e + phi2*v(l)*v(m)/rv^2);
                        H1{l,m}(j,k) = (G(l,m)*dist -   M1{l,m}(j,k)*lg4s);
                        M2{l,m}(j,k) = M2{l,m}(j,k)*dist;
                        H2{l,m}(j,k) = (Tn(l,m) -   M2{l,m}(j,k)*lg4s );
                    end
                end
           
            end
        end
    end

    toc
    disp('参数矩阵生成');
    %% 构造线性系统，方程左端的矩阵 %%
    
    A = zeros(4*n,4*n);

    A(1:2*n,1:2*n)  = R.*(1i*gamma*diag(distance)*M1{1,1}+M2{1,1}) + pi/n*(1i*gamma*diag(distance)*H1{1,1}+H2{1,1})- 1/2* diag(distance);
    A(1:2*n,2*n+1:end) = R.*(1i*gamma*diag(distance)*M1{1,2}+M2{1,2}) + pi/n*(1i*gamma*diag(distance)*H1{1,2}+H2{1,2});
    A(2*n+1:end,1:2*n) = R.*(1i*gamma*diag(distance)*M1{2,1}+M2{2,1}) + pi/n*(1i*gamma*diag(distance)*H1{2,1}+H2{2,1});
    A(2*n+1:end,2*n+1:end) = R.*(1i*gamma*diag(distance)*M1{2,2}+M2{2,2}) + pi/n*(1i*gamma*diag(distance)*H1{2,2}+H2{2,2})- 1/2* diag(distance);
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