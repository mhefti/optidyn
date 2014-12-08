function minpar = readfittinglog(fname)

if nargin < 1 
    fname = 'fitting_log.txt';
end

content = dlmread(fname,'',2,0);
[~,Imin] = min(content(:,end));
minpar = content(Imin,:)';
