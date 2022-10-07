function writematrix(matrix, filename)
[m,n] = size(matrix);

OutputFileName = strcat(filename);
[fid, msg] = fopen(OutputFileName, 'wt');
if fid < 0
    error('Could not open file "%s" because "%s"', fid, msg);
end
for i = 1:1:m
    for j = 1:1:n-1
        fprintf(fid, '%f,',matrix(i,j));
    end
    fprintf(fid, '%f',matrix(i,n));
    fprintf(fid, '\n');
end

fclose(fid);
return