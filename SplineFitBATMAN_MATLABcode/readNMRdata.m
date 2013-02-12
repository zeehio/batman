function Data = readNMRdata(filedir)
% written 120213 Dr. Jie Hao, Imperial College London

fid = fopen(filedir,'r');
tLines = fgets(fid);
delimiter = sprintf('\t','');
numCols = numel(strfind(tLines,delimiter)) + 1;
FormatString=repmat('%f',1,numCols);
InputText=textscan(fid,FormatString,'delimiter','\n');
Data=cell2mat(InputText);
fclose(fid);
