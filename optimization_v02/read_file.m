function struc = read_file(fileName,headerlinesIn,type)

fid = fopen(fileName,'r');

switch type
    
    case 'Tcol'
        sc = textscan(fid,repmat('%f',1,2),...
            'delimiter','\t','HeaderLines',headerlinesIn);
        fclose(fid);
        trel = sc{2}(:);
        Temp = sc{2}(:);
        struc.t = trel;
        struc.T = Temp;
        
    case 'means'
        sc = textscan(fid,repmat('%f',1,10),...
            'delimiter','\t','HeaderLines',headerlinesIn);
        fclose(fid);
        flow = sc{1}(:);
        Press = sc{2}(:);
        yf = sc{9}(:);
        struc.Q = flow;
        struc.P = Press;
        struc.yfeed = yf;
        
    case 'column'
        sc = textscan(fid,repmat('%f',1,8),...
            'delimiter','\t','HeaderLines',headerlinesIn);
        fclose(fid);
        struc.time = sc{1}(:);
        struc.press1 = sc{2}(:);
        struc.press1 = sc{3}(:);
        struc.flow = sc{4}(:);
        struc.temp = sc{5}(:);
        struc.flowliq = sc{6}(:);
        struc.pH2O = sc{7}(:);
        struc.yH2O = sc{8}(:);
        
    case 'dP'
        sc = textscan(fid,repmat('%f',1,9),...
            'delimiter','\t','HeaderLines',headerlinesIn);
        fclose(fid);
        struc.pmeasure = sc{1}(:);
        struc.ppredict = sc{2}(:);
        struc.pcolin = sc{3}(:);
        struc.pcolout = sc{4}(:);
        struc.pcol = sc{5}(:);
        struc.yfeed = sc{6}(:);
        struc.rHfeed = sc{7}(:);
        struc.Tcol = sc{8}(:);
        struc.flow = sc{9}(:);
       
end
        
        
        