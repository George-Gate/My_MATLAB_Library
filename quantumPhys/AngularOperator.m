function [ Lx,Ly,Lz,Lp,Lm ] = AngularOperator( l )
% Generate the matrix representation of angular momentum operators for a given quantum number l
%  [George-Gate @2017-08-13]
%
% [Usage]
%    [ Lx,Ly,Lz,Lp,Lm ] = AngularOperator( l )
% [Input]
%    l   -   The angular momentum quantum number, L^2=l*(l+1)
% [Output]
%    Lx, Ly, Lz  -  The matrix that represent the angular momentum operator  
%                   in different direction. Dim of them is (2*l+1).
%                   They satisfy the commutation relation:
%                        [Li,Lj]=i*e(i,j,k)*Lk
%                   Where we set hbar=1
%    Lp, Lm  -  The matrix of upper/lower operator L_+ and L_-
%                   Lp = Lx + i * Ly
%                   Lm = Lx - i * Ly
%
% [State Coding]
%    There are (2*l+1) states involving: 
%           |l,ml>,  ml=-l, -l+1, ..., l
%    Where Lz*|l,ml>=ml*|l,ml>
%    We encode 
%         |l,l> as [1;0;0;...]
%         |l,l-1> as [0;1;0;...]
%         ......
%         |l,-l> as [0;0;...;0;1]
%    So Lz=diag([l, l-1 ,l-2, ..., -l]);
%
    Dim=2*l+1;
    eigs=l:-1:-l;
    Lz=diag(eigs);
    Lp=zeros(Dim);Lm=zeros(Dim);
    % Make Lp and Lm
    for i=1:Dim
        ml=eigs(i);
        if (i>1)
            Lp(i-1,i)=sqrt(l*(l+1) - ml*(ml+1));
        end
        if (i<Dim)
            Lm(i+1,i)=sqrt(l*(l+1) - ml*(ml-1));
        end
    end
    % Calc Ly, Lz from Lp, Lm
    Lx=(Lp+Lm)/2;
    Ly=(Lp-Lm)/2i;
end

