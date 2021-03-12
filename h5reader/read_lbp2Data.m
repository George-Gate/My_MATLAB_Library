% 读取 Newport LBP2 系列 光束质量分析仪 软件保存的图像数据。
% 该图像数据的文件扩展名为 ".lbp2Data"，本质上是一个h5文件。
%
% Author: GeorgeGate @ 2021-03-12
%
function mapData = read_lbp2Data(filename, showfigure)
    if nargin < 2
        if nargout == 0
            showfigure = true;
        else
            showfigure = false;
        end
    end
    info = h5info(filename);
    numFigs = length(info.Groups(1).Groups);
    mapData.BG_SETUP = h5readAll(filename, '/BG_SETUP');
    for ii = 1:numFigs
        mapData.BG_DATA(ii) = h5readAll(filename, ['/BG_DATA/' num2str(ii)]);
    end
    w = double(mapData.BG_DATA(1).RAWFRAME.WIDTH);
    h = double(mapData.BG_DATA(1).RAWFRAME.HEIGHT);
    dx = mapData.BG_DATA(1).RAWFRAME.PIXELSCALEXUM;
    dy = mapData.BG_DATA(1).RAWFRAME.PIXELSCALEYUM;
    maps = cell(1,numFigs);
    for ii = 1:numFigs
        maps{ii} = reshape(mapData.BG_DATA(ii).DATA, w, h);
    end
    mapData.intensity = maps;
    [mapData.x, mapData.y] = ndgrid(0:dx:dx*(w-1), 0:dy:dy*(h-1));
    
    if showfigure
        figure();
        for ii = 1:numFigs
            plotIntensity(subplot(4,4,ii),mapData.x, mapData.y, maps{ii}, sprintf('Intensity Map %d', ii));
        end
    end
end

function plotIntensity(ax, x,y,z,tag)
    s = pcolor(ax, x,y,z);
    s.EdgeColor = 'none';
    colormap(ax, 'hot');
    xlabel(ax, 'x / um');ylabel(ax, 'y / um');
    title(ax, tag);
end
