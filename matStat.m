function [ varargout ] = matStat( A )
%% Display the number of non-zero elements in a Matrix
%  [George-Gate @2017-05-06]
%  [Usage]
%       matDensity=matStat( A )
%   

    matDensity=nnz(A)/numel(A);
    disp(['Number of non-zero elements in ',inputname(1),': ',num2str(nnz(A)),',  ',num2str(100*matDensity),'%']);
    if nargout ==1
        varargout{1}=matDensity;
    end
end

