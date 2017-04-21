function [ cI, cx, cy, cz ] = pauliDecomposite( M )
%Give a 2 x 2 matrix M, decomposite it on the pauli matrix base.
%  [George-Gate @2017-04-10]
%   cI = 0.5 * trace(M)
%   cx = 0.5 * trace(sigma_x * M)
%   cy = 0.5 * trace(sigma_y * M)
%   cz = 0.5 * trace(sigma_z * M)

    if ( ~isequal(size(M),[2,2]) )
        error('Matrix dimension invalid!');
    end
    sx=[0,1;1,0];sy=[0,-1i;1i,0];sz=[1,0;0,-1];
    cI=0.5 * trace(M);
    cx=0.5 * trace(sx * M);
    cy=0.5 * trace(sy * M);
    cz=0.5 * trace(sz * M);
    display(['[cI, cx, cy, cz]=',num2str(cI),', ',num2str(cx),', ',num2str(cy),', ',num2str(cz)]);
end

