function [ gp_x, gw ] = getGaussPts( ngp )
%Find the gauss points and weigths for Gauss¨CLegendre quadrature
% [George-Gate @2017-04-21]
%
% --------------- Generate gauss points and weights --------------
% gp_x: gauss points; gw: weights
    m=(1:ngp-1)';
    c=m./sqrt(4*m.*m-1);
    Jm=spdiags([[c;0],[0;c]],[-1,1],ngp,ngp);
    [Vec,Dig]=eigs(Jm,ceil(ngp/2),'sa');
    gp_x=diag(Dig);
    % normalize V
    Vec=Vec./repmat(sqrt(sum(Vec.*Vec,1)),ngp,1);
    gw=2*abs(Vec(1,:).*Vec(1,:))';
    if mod(ngp,2)==1
        gp_x=[gp_x;-gp_x(end-1:-1:1)];
        gw=  [gw  ; gw(  end-1:-1:1)];
    else
        gp_x=[gp_x;-gp_x(end:-1:1)];
        gw=  [gw  ; gw(  end:-1:1)];
    end
end

