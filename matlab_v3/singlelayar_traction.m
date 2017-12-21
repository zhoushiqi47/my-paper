function A=singlelayar_traction(omega, lambda, mu, n, bctype)

w = quad_weights(n);
R = zeros(2*n);

for k=1:2*n
     idx=k:2*n;
     R(idx,k)=w(1:2*n-k+1);
     R(k,k)=R(k,k)/2;  %% for convinience
end
R=(R+R');
lamda=lambda;
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

    cp = sqrt( lambda + 2*mu);
    cs = sqrt( mu );
    cp2 = lambda + 2*mu;
    %% the corresponding wavenumber;
    kp = omega / cp;
    ks = omega / cs;
    
    CC = mu/(lambda+2*mu)/(2*pi);

    %% constant definition
    w2 = 1/omega^2 ;

    M = cell(2,2);
    H = cell(2,2);
    for k=1:4
    M{k}=zeros(2*n);
    H{k}=zeros(2*n);
    end
    

    % 将矩阵拉开成向量
    xl = repmat([x1,x2]',1,2*n);
    yl = reshape( repmat([x1,x2]',2*n,1),2,4*n*n);
    rGG  = TractionNGreenTensor_2D(omega,kp,ks,yl,xl);



    for j=1:2*n
        dx1j=dx1(j); dx2j=dx2(j);
        for k=1:2*n
            if (j==k)
                dist = distance(k);
                
                T = rGG(:,(k-1)*2*n+j);
                T1 = reshape(T(1:4),2,2); 
                T2 = reshape(T(5:8),2,2);
                Tn = (T1*dx2j-T2*dx1j)*dist;
           
                
                beta = 1/2*(ddx1(j)*dx2j-ddx2(j)*dx1j)/dist;
            
                H{1,1}(j,k) = Tn(1,1)+(2*(lambda+mu)*dx1j^2/dist^2/mu+1)*beta*CC;
                H{2,2}(j,k) = Tn(2,2)+(2*(lambda+mu)*dx2j^2/dist^2/mu+1)*beta*CC;
            
                H{1,2}(j,k) = Tn(1,2)+(2*(lambda+mu)*dx1j*dx2j/dist^2/mu)*beta*CC-1/2*(ddx1(j)*dx1j+ddx2(j)*dx2j)/dist*CC;
                H{2,1}(j,k) = Tn(2,1)+(2*(lambda+mu)*dx1j*dx2j/dist^2/mu)*beta*CC+1/2*(ddx1(j)*dx1j+ddx2(j)*dx2j)/dist*CC;
            else
                
                v = ([x1(j)-x1(k) x2(j)-x2(k)]); 
                z1 = v(1);
                z2 = v(2);
                lg4s = log(4*sin((t(j)-t(k))/2)^2) ; 
                rv = norm(v);
                dist = distance(k);
                
                T = rGG(:,(k-1)*2*n+j);
                T1 = reshape(T(1:4),2,2); 
                T2 = reshape(T(5:8),2,2);
                Tn = (T1*dx2j-T2*dx1j)*dist;
                
                Hp0 = besselj(0,kp*rv);
                Hp1 = besselj(1,kp*rv);

                Hs0 = besselj(0,ks*rv);
                Hs1 = besselj(1,ks*rv);

                r2 = rv*rv;
                r3= r2*rv;

                phis1 = -1/(4*pi)*w2*(-ks^2*Hs0 + 2*ks*Hs1./rv);
                phip1 = -1/(4*pi)*w2*(-kp^2*Hp0 + 2*kp*Hp1./rv);

                phis2 = -1/(4*pi)*w2*(ks^3*Hs1 + 4*ks^2*Hs0./rv - 8*ks*Hs1./r2);
                phip2 = -1/(4*pi)*w2*(kp^3*Hp1 + 4*kp^2*Hp0./rv - 8*kp*Hp1./r2);

                gs111 = 3*phis1.*z1./r2 + phis2.*z1.^3./r3;
                gs112 = phis1.*z2./r2 + phis2.*z1.^2.*z2./r3;
                gs122 = phis1.*z1./r2 + phis2.*z2.^2.*z1./r3;
                gs222 = 3*phis1.*z2./r2 + phis2.*z2.^3./r3;

                gp111 = 3*phip1.*z1./r2 + phip2.*z1.^3./r3;
                gp112 = phip1.*z2./r2 + phip2.*z1.^2.*z2./r3;
                gp122 = phip1.*z1./r2 + phip2.*z2.^2.*z1./r3;
                gp222 = 3*phip1.*z2./r2 + phip2.*z2.^3./r3;
                
                M{1,1}(j,k)=(-2*mu*gs122 - cp2*gp111 - lamda*gp122)*dx2j-(mu*(gs112-gs222) - 2*mu*gp112)*dx1j;
                M{2,1}(j,k)=(mu*(gs112-gs222) - 2*mu*gp112)*dx2j-(2*mu*gs122 - cp2*gp122 - lamda*gp111)*dx1j;
                M{1,2}(j,k)=(2*mu*gs112 - cp2*gp112 - lamda*gp222)*dx2j-( mu*(gs122-gs111) - 2*mu*gp122)*dx1j;
                M{2,2}(j,k)=(mu*(gs122-gs111) - 2*mu*gp122)*dx2j-( -2*mu*gs112 - cp2*gp222 - lamda*gp112)*dx1j;



                for l=1:2
                    for m=1:2
                        M{l,m}(j,k) = M{l,m}(j,k)*dist;
                        H{l,m}(j,k) = (Tn(l,m) -   M{l,m}(j,k)*lg4s );
                    end
                end
            end
        end
    end

    %% 构造线性系统，方程左端的矩阵 %%
    
    A = zeros(4*n,4*n);

    A(1:2*n,1:2*n)  = R.*M{1,1} + pi/n*H{1,1};
    A(1:2*n,2*n+1:end) = R.*M{1,2} + pi/n*H{1,2} ;
    A(2*n+1:end,1:2*n) = R.*M{2,1} + pi/n*H{2,1} ;
    A(2*n+1:end,2*n+1:end) = R.*M{2,2} + pi/n*H{2,2};
end
