function minpar = readfittinglog(fname)

if nargin < 1 
    fname = 'fitting_log.txt';
end

delimiter = sprintf('\t','');
fid = fopen(fname,'rt');
tLines = fgets(fid);
numCols = numel(strfind(tLines,delimiter));
frewind(fid);
content = textscan(fid,repmat('%f',1,numCols),'HeaderLines',1);
fclose(fid);
content = [content{:}];
% find the minimum, residual is in the last column
[~,Imin] = min(content(:,numCols));
minpar = content(Imin,:)';
