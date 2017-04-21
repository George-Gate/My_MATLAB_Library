function mat = vec2full( iList, jList, vList, m, n )
% Convert a vector list to a full matrix
% [ George-Gate @2017-04-21]
%   The usage is very similar to sparse
%
%   mat=vec2full(i,j,v,m,n);
%   i,j,v is col vectors
%   mat is a (m x n) matrix with all non-zero elements given as (i,j,v)
%   if there is two lines that have the same (i,j), the program will sum those v to the same position.
    mat=zeros(m,n);
    if ~(iscolumn(iList) && iscolumn(jList) && iscolumn(vList) && length(iList)==length(jList) && length(jList)==length(vList))
        error('i,j,v must be column vectors of the same size.');
    end
    if max(iList)>m || min(iList)<1 || max(jList)>n || min(jList)<1
        error('iList/jList out of range.');
    end
    for i=1:length(iList)
        mat(iList(i),jList(i))=mat(iList(i),jList(i))+vList(i);
    end
end

