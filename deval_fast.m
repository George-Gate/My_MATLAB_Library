function [Sxint,Spxint] = deval_fast(sol,xint)
%DEVAL_FAST 用于加速deval的执行
%    当数组 xint 非常大时，deval_fast(sol,xint) 会显著快于 deval(sol,xint)
%    
% Author: GeorgeGate
% Last Revision: 2020/09/28

    if ~isa(sol,'struct')  
      % Try  DEVAL(XINT,SOL)
      temp = sol;
      sol  = xint;
      xint = temp;
    end
    
    % Evaluate the first derivative?
    Spxint_requested = (nargout > 1);   
    
    if length(xint) < 5e4
        if Spxint_requested
            [Sxint,Spxint] = deval(sol,xint);
        else
            Sxint = deval(sol,xint);
        end
        return;
    end
        
    try
      t = sol.x;    
    catch
      error(message('MATLAB:deval:SolNotFromDiffEqSolver', inputname( 1 )));
    end
    
    % Make tint a row vector and if necessary, 
    % sort it to match the order of t.
    tint = reshape(xint, 1, []);  
    tdir = sign(t(end) - t(1));
    had2sort = any(tdir*diff(tint) < 0);
    if had2sort
        if tdir>0
            [tint,tint_order] = sort(tint, 'ascend');
        else
            [tint,tint_order] = sort(tint, 'descend');
        end
    end
    
    % 给tint分组，并调用deval
    max_it = findMaxID(t, 1, length(t), tint(end), 0.5, tdir) + 1;
    max_it = min(max_it, length(t));
    
    id_start = 1;
    id_end = 1;
    xint_grouped = cell(1, max_it-1 );
    nCell = 1;
    for it = 2:max_it
        % find the max id_end, such that tint(id_end) <= t(it)
        stepSizeRatio = min(0.8,1.05*(t(it) - t(it-1)) / (  t(max_it) - t(it-1)  ));  % 根据当前区间占总区间的长度来估计答案的位置，这里假设了xint是均匀的。
        id_end = findMaxID(tint, id_end, length(tint), t(it), stepSizeRatio, tdir);
        xint_grouped{nCell} = reshape(tint(id_start:id_end),1,[]);
        id_start = id_end+1;
        if ~isempty(xint_grouped{nCell})
            nCell = nCell + 1;
        end
        if id_start > length(tint)
            break;
        end
    end
    xint_grouped = xint_grouped(1:nCell-1);

    
    % call deval
    if Spxint_requested
        [Sxint,Spxint] = cellfun(@(x)deval(sol,x), xint_grouped, 'UniformOutput', false);
        Spxint = cell2mat(Spxint);
    else
        Sxint = cellfun(@(x)deval(sol,x), xint_grouped, 'UniformOutput', false);
    end
    Sxint = cell2mat(Sxint);
    
    if had2sort   % Restore the order of tint in the output.
        Sxint(:,tint_order) = Sxint;
        if Spxint_requested
            Spxint(:,tint_order) = Spxint;
        end
    end
    
end

function ID = findMaxID(array, iS, iE, target, stepSizeRatio, direction)

    cmpFun = @(a,b) direction*a <= direction*b;

    step = ceil(stepSizeRatio * (iE-iS));
    if step < 1   % 结束条件： iE <= iS
        ID = iE;
        return;
    end
    iMid = iS + step;
    if iMid > iE
        error('iMid > iE, which means stepSizeRatio > 1.');
    end
    if abs(iMid - iE)<0.1    % 结束条件：iMid == iE
        ID = iMid;
        return;
    end
    iMid_status = cmpFun(array(iMid), target);
    if iMid_status && ~cmpFun(array(iMid+1), target)    % 结束条件：array(iMid) 刚好满足条件
        ID = iMid;
        return;
    end
    if iMid_status
        ID = findMaxID(array, iMid, iE, target, stepSizeRatio, direction);   % 一开始，先以输入的 stepSizeRatio 往前找到一个 cmpFun(array(iE), target) == false 的范围 [iS, iE]
    else
        ID = findMaxID(array, iS, iMid-1, target, 0.5, direction);           % 然后开始二分搜索
    end
end