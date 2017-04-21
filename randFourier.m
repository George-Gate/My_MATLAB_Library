% [George-Gate @2016-04-20]
% 产生一个随机的N维Fourier Matrix
% 只在Fourier Matrix的一个子集中随机，有些Fourier Matrix永远都不会出现。
% Fjk=1/sqrt(N)*exp(2i*pi*j*k/N)*exp(1i*(theta_j+phi_k));

function F=randFourier( N )
    [A,B]=meshgrid(1:N,1:N);
    C=exp(2i*pi*A.*B/N);
    A=diag( exp(2i*pi*rand(N,1)) );
    B=diag( exp(2i*pi*rand(N,1)) );
    F=A*C*B/sqrt(N);    
end

