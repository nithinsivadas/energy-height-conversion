function write_h5_dataset(h5FileStr, datasetPath, varValue,...
    timeDim, setAppend, setComment)
%% write_h5_dataset Writes a variable containing data to hdf5 file.
%  Input:
%         h5FileStr   - Name/Path of h5 file
%         datasetPath - dataset path and name inside the hdf5 file
%         varValue    - variable containing the data
%         timeDim     - The index of the time dimension (if it is 0, then it
%                       is considered that there is no time dimension)
%         setAppend   - If true, appends on the existing data set
%         setComment  - If true, will fprintf stages of the process.
if nargin < 6
    setComment = false;
end

if nargin < 5
    setAppend = true;
end

if nargin < 4 || isempty(timeDim)
    timeDim = 0;
end

setAppendOriginal = setAppend;

[status, info, ~] = ish5dataset(h5FileStr, datasetPath);

dataSize = size(varValue);
if timeDim > 0
    nTime = dataSize(timeDim);
    nDim = length(dataSize);
    varValue = permute(varValue,fliplr(1:1:nDim));
    dataSize = size(varValue);
    timeDim = find(dataSize == nTime);
else
    nDim = length(dataSize);
    varValue = permute(varValue,fliplr(1:1:nDim));
    dataSize = size(varValue);
end

if timeDim>0

    if ~status
        setAppend = true;
        comment('Creating a new variable \n',setComment);
        dataMaxSize = dataSize;
        dataMaxSize(timeDim) = Inf;
        chunkSize = dataSize;
        chunkSize(timeDim) = ceil(chunkSize(timeDim)/4);
        h5create(h5FileStr, datasetPath, dataMaxSize, 'ChunkSize', chunkSize,...
            'Deflate', 9);
        [status, info, ~] = ish5dataset(h5FileStr, datasetPath);
        dataMaxSizeOld = info.Dataspace.MaxSize;
        dataSizeOld = info.Dataspace.Size;
        nTimeOld = 0;
    else
        dataMaxSizeOld = info.Dataspace.MaxSize;
        dataSizeOld = info.Dataspace.Size;
        nTimeOld = dataSizeOld(timeDim);
        nonTimeDim = find(dataSizeOld ~= nTimeOld);
        if ~isequal(dataSizeOld(nonTimeDim),dataSize(nonTimeDim))
            error('The non-time dimensions of hdf5 data do not match with input');
        end
    end
    % Now Dataset exists
    nDimOld =  length(dataSizeOld);
    if nDim~=nDimOld
        error('The dimensions of the existing data and varValue do not match');
    end
    if ~any(isinf(dataMaxSizeOld))
        warning('Cannot append');
        setAppend = false;
    end

    if setAppend
        startIndx = ones(1,nDim);
        countIndx = dataSize;
        countIndx(timeDim) = nTime;
        startIndx(timeDim) = nTimeOld+1;
        h5write(h5FileStr,datasetPath,...
                        varValue,startIndx,countIndx);
    elseif ~setAppendOriginal
        comment('Rewriting the variable, since requested \n',setComment);
        h5write(h5FileStr,datasetPath,varValue);
    else
        error ('Cannot append because either time or varValue datasets do not allow it');
    end
else
    if ~status
        comment('Creating a new variable \n',setComment);
        dataMaxSize = dataSize;
        chunkSize = ceil(dataSize/4);
        h5create(h5FileStr, datasetPath, dataMaxSize, 'ChunkSize', chunkSize, 'Deflate', 9);
    end
    comment('Rewriting the variable, since variable time independent \n',setComment);
    h5write(h5FileStr,datasetPath,varValue);
end

comment(['Succeeded in writing ',datasetPath,' \n'],setComment);

end

function comment(str,setComment)

    if setComment
        fprintf(str);
    end

end
