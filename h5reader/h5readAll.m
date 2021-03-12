% 将h5文件中的内容读取到matlab内存，并以struct的形式返回。
% filename: h5文件的文件名
% rootPath: 要读取的内容在h5文件中的路径。输入 rootPath = '/' 完整整个h5文件的全部内容。
%           默认值：'/'
%
% Author: GeorgeGate @ 2021-03-12
%
%
function fulldata = h5readAll(filename, rootPath)
    if nargin < 2
        rootPath = '/';
    end
    if ~strcmp(rootPath(end) , '/')
        rootPath = [rootPath '/'];
    end
    info = h5info(filename, rootPath);
    fulldata = struct();
    % read datasets
    for ii = 1:length(info.Datasets)
        name = validateName(info.Datasets(ii).Name);
        res = h5read(filename,[rootPath '/' name]);
        if iscell(res) && isscalar(res) && (isstring(res{1}) || ischar(res{1}))
            fulldata.(name) = res{1};
        else
            fulldata.(name) = res;
        end
    end
    % read groups
    for ii = 1:length(info.Groups)
        fullPath = info.Groups(ii).Name;
        name = split(fullPath, '/');
        name = validateName(name{end});
        fulldata.(name) = h5readAll(filename, fullPath);
    end

end

function validName = validateName(name)
    validName = replace(name,' ','_');
    validName = matlab.lang.makeValidName(validName);
end