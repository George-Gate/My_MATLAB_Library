function [ varargout ] = matStat( A, name )
%% Display the number of non-zero elements in a Matrix
%  [George-Gate @2017-06-14]
%  [Usage]
%       matDensity=matStat( A )
%   
    if nargin==1
        name=inputname(1);
    end
    nnzA=nnz(A);
    matDensity=nnzA/numel(A);
    disp(['Number of non-zero elements in ',name,': ',num2str(nnzA),',  ',num2str(100*matDensity,'%.2f'),'%'...
         ,',  ',num2str(nnzA/size(A,1),'%.1f'),' per row on average.']);
    if nargout ==1
        varargout{1}=matDensity;
    end
end

