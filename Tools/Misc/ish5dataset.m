function [status,info,ME]=ish5dataset(h5FileStr,dataSetPath)
        % Finds if the particular dataset exists or not

        info = 0;
        ME = 'No error';
        try
            info=h5info(h5FileStr,dataSetPath);
        catch ME
        end

        if isnumeric(info)
            if info==0
            status = false;
            end
        elseif isstruct(info)
            status = true; 
        else
            status = 'unknown';
        end
        
end
