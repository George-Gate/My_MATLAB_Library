%  [George-Gate @2016-04-26]
%% method 1
% 考虑一个N维球面，可以用(N-1)个角度坐标来定位一个点
% (cos1*sin2*sin3*sin4*sin5,
%  sin1*sin2*sin3*sin4*sin5,
%       cos2*sin3*sin4*sin5,
%            cos3*sin4*sin5,
%                 cos4*sin5,
%                      cos5)
% 问题：概率分布不均匀
N=64;
theta=pi/2*rand(N-1,1);
dRho=ones(N,1);
dRho(1)=cos(theta(1));dRho(2)=sin(theta(1));
for i=2:N-1
    for k=1:i
        dRho(k)=dRho(k)*sin(theta(i));
    end
    dRho(i+1)=cos(theta(i));
end
rho=diag(dRho.^2);
if rand()<0.5
    U=randU(N);
else
    U=randFourier(N);
end
rho=U*rho*U';

%% method 2
% 将很多个随机的纯态加起来，得到一个混态
% 权重是一个等比数列
% 问题：貌似有些密度矩阵是永远都不可能取到的(概率超级低)
N=64;
rho=zeros(N);
ratio=rand();
for i=1:min(1000,max(2,ceil(log(1e-1*(1-ratio))/log(ratio))))
    psi=randPsi(N);
    rho=ratio*rho+(psi*psi');
end
rho=rho./trace(rho);

%% method 3
% 这个方法随出来的矩阵范围貌似最大
% 而且很快
N=64;
N0=round(N*rand()-0.5)+1;
dRho=[rand(N0,1);zeros(N-N0,1)];
dRho=dRho(randperm(N));
dRho=dRho/sum(dRho);
rho=diag(dRho);
if rand()<0.5
    U=randU(N);
else
    U=randFourier(N);
end
rho=U*rho*U';

%% method 4
% 这个方法随出来的矩阵范围很小，但是其它方法都不太容易随出这些矩阵
N=64;
rho=randRho(N);
