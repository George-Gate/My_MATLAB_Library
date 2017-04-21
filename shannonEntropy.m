% [George-Gate @2016-03-02]
%
%  E = shannonEntropy( p, base )
%    Calculate the Shannon Entropy of a given probability distribution.
%    'base' is the base of the log. Default value is 2.
%    The program will normalize the distribution first.

function E = shannonEntropy( p, base )
    if (nargin<2)
        base=2;
    end
    p=p(:);
    p=p(abs(p)>eps);
    p=p/sum(p);
    E=-sum(p.*log(p))/log(base);
end

