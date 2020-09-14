function newGrid = refineGrid_1D( oldGrid, ins_pts)
% Refine a 1D grid. i.e. make the grid more dense.
% [ George-Gate @ 2020-09-15]
%
%   newGrid = refineGrid_1D( oldGrid, ins_pts );
%   
%       oldGrid: old 1D grid.
%       ins_pts:   an integer, will insert ins_pts points between each two grid points of oldGrid.
%                  
%                length(newGrid) - 1 = ( length(oldGrid) - 1 ) * ( ins_pts + 1 )
%
    p = inputParser;
    p.addRequired('oldGrid', @(x)validateattributes(x,{'numeric'},{'vector','finite', 'nonnan', 'real'}));
    p.addRequired('ins_pts', @(x)validateattributes(x,{'numeric'},{'scalar','integer', 'positive'}));
    p.parse(oldGrid, ins_pts);
    
    sec_pts = ins_pts + 1;
    sec_grid = arrayfun(@(a,b)linspace(a,b,sec_pts+1).', oldGrid(1:end-1), oldGrid(2:end), 'UniformOutput', false);
    sec_grid = cell2mat(sec_grid);
    if isrow(oldGrid)
        newGrid = [reshape(sec_grid(1:sec_pts,:),1,[]), oldGrid(end)];
    else
        newGrid = [reshape(sec_grid(1:sec_pts,:),[],1); oldGrid(end)];
    end
    
end

