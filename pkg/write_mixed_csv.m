function write_mixed_csv(filedir, strings)
% written 120213 Dr. Jie Hao, Imperial College London
fid = fopen(filedir,'w');

fmtString = [repmat('%s,',1,size(strings,2)-1),'%s\n'];
for i = 1:size(strings,1)
    for j = 1:size(strings,2)
        st{j,i} = ['"' strings{i,j} '"'];
    end
end
fprintf(fid,fmtString,st{:});
fclose(fid);