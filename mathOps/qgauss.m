function I = qgauss( func, x1, x2, ngp )
%Perform Gauss¨CLegendre quadrature for 1-D function
% [George-Gate @2017-04-21]
%********************************************************************************
%     Returns the SINGLE integral of the function (of ONE VARIABLE) "func"
% between x1 and x2 by N-point Gauss-Legendre integration. The function
% is evaluated exactly N times at interior points in the range of
% integration.  
%********************************************************************************
    [gp_x, gw] = getGaussPts( ngp );
    
    I=(x2-x1)/2*sum(func((gp_x+1)/2*(x2-x1)+x1).*gw);

end

