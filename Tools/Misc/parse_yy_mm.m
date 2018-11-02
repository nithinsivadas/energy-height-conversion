function arrayYYMM=parse_yy_mm(minTimeStr, maxTimeStr,varargin)
    p = inputParser;
    addParameter(p,'Format','dd-mm-yyyy',@(x) isstring(x)||ischar(x));
    addRequired(p,'minTimeStr',@(x) isstring(x)||ischar(x));
    addRequired(p,'maxTimeStr',@(x) isstring(x)||ischar(x));
    parse(p,minTimeStr,maxTimeStr,varargin{:});
    format = p.Results.Format;
    
    startyy = year(minTimeStr, format);
    startmm = month(minTimeStr, format);
    endyy = year(maxTimeStr, format);
    endmm = month(maxTimeStr, format);
    k = 1;
    imm=startmm;
    while datetime(startyy,imm,1)<=datetime(endyy,endmm,1)
        arrayYYMM(1,k) = year(datetime(startyy,imm,1));
        arrayYYMM(2,k) = month(datetime(startyy,imm,1));
        imm = imm+1;
        k = k+1;
    end
end