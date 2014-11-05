function customwrite(filename,output,delimiter,mode)
% original (cell2csv.m) by Sylvain Fiedler, KA, 2004
% modified by Rob Kohr, Rutgers, 2005 - changed to english and fixed delimiter
% modified and adapted by Matthias Scharfe, 2014 simpler code

%mode: 'w': create new file and overwrite
%'a': append

% Open the file
fid = fopen(filename,mode);

if fid == -1
    % could not open the file
    error('Error opening file')
end

for z = 1:size(output,1)
    for s = 1:size(output,2)
        
        var = output{z,s};
        
        if size(var,1) == 0
            var = '';
        end
        
        if isnumeric(var) == 1
            var = num2str(var);
        end
        
        fprintf(fid,var);
        
        if s ~= size(output,2)
            fprintf(fid,delimiter);
        end
    end
    fprintf(fid,'\n');
end
fclose(fid);

end

