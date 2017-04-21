% [George-Gate @2016-04-20]
% ����һ�������NάFourier Matrix
% ֻ��Fourier Matrix��һ���Ӽ����������ЩFourier Matrix��Զ��������֡�
% Fjk=1/sqrt(N)*exp(2i*pi*j*k/N)*exp(1i*(theta_j+phi_k));

function F=randFourier( N )
    [A,B]=meshgrid(1:N,1:N);
    C=exp(2i*pi*A.*B/N);
    A=diag( exp(2i*pi*rand(N,1)) );
    B=diag( exp(2i*pi*rand(N,1)) );
    F=A*C*B/sqrt(N);    
end

