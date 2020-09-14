function FockOperator(Nmodes, max_excite)
%% Description
% Generate following operators in fock representation: 
%     a1, a2
% Generate two mapping:
%     k2nn, nn2k
% Calc the dimension of Hilbert space: 
%     Dim
% Required input: (Should be set in workspace before calling this script)
%     N - maximum number of bosons
%
% [Encoding] 
%-------------------------------------------------------------------------------
%    [1;0;0;....] -> |0,0>  [0;1;0;...] -> |1,0> ... [0;...;0;1;0;...] -> |N,0>
%    [0;...;0;1;0;...] -> |0,1>  [0;...;0;1;0;...] -> |1,1> ... [...] -> |N-1,1>
%       (N+1 leading 0) 
%    [0;...;0;1;0;...] -> |0,2>  ...    [...] -> |N-2,2>
%       (2N+1 leading 0)
%    ......
%    [0;...;0;1;0;...] -> |0,k>  ...    [...] -> |N-k,k>    
%       ( k*N+(3k-k^2)/2 leading 0)
%-------------------------------------------------------------------------------
% Total Dim: (N+1)*(N+2)/2
% |n1,n2> --> psi( n2*N+(3*n2-n2^2)/2+n1+1 )=1   (n1+n2<=N)
%             psi( nn2k(n1+1,n2+1)  )=1
% psi( k )=1  --> |k2nn(1,k),k2nn(2,k)>
%
% Under some circumstance, a1 and a2' may not commute. For example, |k,N-k>,
% if we do a2' first, |k,N-k+1> has N+1 boson which is not included in our
% numerical representation. But if we do a1 first and then do a2', no error
% will occur. To avoid such circumstance, try to start the simulation from 
% N-1 bosons. Or always act annihilation before creation. 

Dim=(N+1)*(N+2)/2;
a1=spalloc(Dim,Dim,Dim);
a2=spalloc(Dim,Dim,Dim);
k2nn=zeros(2,Dim);
nn2k=zeros(N+1,N+1);
% make a1 and a2
for n1=0:N
    for n2=N-n1:-1:0
        i=n2*N+(3*n2-n2^2)/2+n1+1;
        j=(n2-1)*N+(5*n2-n2^2-4)/2+n1+1;
        k=n2*N+(3*n2-n2^2)/2+(n1-1)+1;
        if n1~=0 
            a1(k,i)=sqrt(n1); %#ok<SPRIX>
        end
        if n2~=0
            a2(j,i)=sqrt(n2); %#ok<SPRIX>
        end
        % record mapping & inverse mapping
        nn2k(n1+1,n2+1)=i;
        k2nn(1,i)=n1;
        k2nn(2,i)=n2;
    end
end
clear n1 n2 i j k;