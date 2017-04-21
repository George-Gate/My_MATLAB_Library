%  [George-Gate @2016-04-26]
%% method 1
% ����һ��Nά���棬������(N-1)���Ƕ���������λһ����
% (cos1*sin2*sin3*sin4*sin5,
%  sin1*sin2*sin3*sin4*sin5,
%       cos2*sin3*sin4*sin5,
%            cos3*sin4*sin5,
%                 cos4*sin5,
%                      cos5)
% ���⣺���ʷֲ�������
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
% ���ܶ������Ĵ�̬���������õ�һ����̬
% Ȩ����һ���ȱ�����
% ���⣺ò����Щ�ܶȾ�������Զ��������ȡ����(���ʳ�����)
N=64;
rho=zeros(N);
ratio=rand();
for i=1:min(1000,max(2,ceil(log(1e-1*(1-ratio))/log(ratio))))
    psi=randPsi(N);
    rho=ratio*rho+(psi*psi');
end
rho=rho./trace(rho);

%% method 3
% �������������ľ���Χò�����
% ���Һܿ�
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
% �������������ľ���Χ��С������������������̫���������Щ����
N=64;
rho=randRho(N);
