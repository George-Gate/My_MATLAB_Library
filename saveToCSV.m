function saveToCSV( filename, varList )
% George-Gate @ 2016-03-15
% ��CSV��ʽ��varList�еı���������filename��ָ���ļ���
% varListΪcell���ͣ�����ֻ��Ҫ�г�������
% Ŀǰֻ֧�����ּ��ı�����
    fid=fopen(filename,'wt');
    % write title
    fprintf(fid,'varName,values\n');
    % print vars
    for i=1:length(varList)
        fprintf(fid,'%s,%s\n',varList{i},evalin('caller',['num2str(',varList{i},',16)']));
    end
    fclose(fid);
end

