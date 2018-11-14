function site=parse_thg_location_xls(fileName)
    % Parse the location of themis GBO sites from an xls table found online
    % here: http://themis.ssl.berkeley.edu/images/ASI/THEMIS_ASI_Station_List_Nov_2011.xls
    if ispc
        [~,site.code] = xlsread(fileName,'A2:A25');
        [~,site.name] = xlsread(fileName,'F2:F25');
        [site.glat] = xlsread(fileName,'D2:D25');
        [site.glon] = xlsread(fileName,'E2:E25');
        [site.mlat] = xlsread(fileName,'G2:G25');
        [site.mlon] = xlsread(fileName,'H2:H25');  
    else
        [num,txt] = xlsread(fileName);
        site.code = cell(txt(2:end,1));
        site.name = cell(txt(2:end,6));
        site.glat = num(2:end,1);
        site.glon = num(2:end,2);
        site.glat = num(2:end,1);
        site.mlon = num(2:end,4);
        site.mlat = num(2:end,5);
    end
end