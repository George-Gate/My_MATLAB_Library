function saveToCSV( filename, varList )
% George-Gate @ 2016-03-15
% 以CSV格式将varList中的变量保存至filename所指的文件中
% varList为cell类型，其中只需要列出变量名
% 目前只支持数字及文本变量
    fid=fopen(filename,'wt');
    % write title
    fprintf(fid,'varName,values\n');
    % print vars
    for i=1:length(varList)
        fprintf(fid,'%s,%s\n',varList{i},evalin('caller',['num2str(',varList{i},',16)']));
    end
    fclose(fid);
end

