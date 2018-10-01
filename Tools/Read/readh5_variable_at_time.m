%% Function to read h5 variable at a particular time index
function varValue=readh5_variable_at_time(h5FileStr,varStr,address,thisTimeIndx)
    % readh5_variable_at_time.m Reads variable name ( in varStr) of the hdf5
    % file ****_energyFlux.h5. 
    %----------------------------------------------------------------------
    % Input
    %--------     
    %   h5FileStr   - The path of the hdf5 file
    %   varStr      - Variable name
    %   address     - H5 file address of the variable
    %   thisTimeIndx- Index of the time of interest
    %
    %----------------------------------------------------------------------
    % Output
    %--------
    %   varValue    - Data within the variable for that time instant
    %----------------------------------------------------------------------
    % Last Updated: 30th Sep 2018
    % Create by   : Nithin Sivadas
    
    array1DwithTimeVars = {'message'};
    matrix2DwithTimeVars = {'MSE','maxIter','time','BmagEq'};
    matrix3DwithTimeVars = {'dEnergyFlux','energyFlux','flux','qError',...
        'qInput','qInverted','Ne','dNeFrac','ASI','BgeoEq','gradBmagEq',...
        'magEqCoordGEO'};
    matrix4DwithTimeVars = {'diffBEq'};
    matrix2D = {'azCalData','elCalData','az','el','lat','lon',...
        'sensorloc','A','alt','energyBin','projectionAlt',...
        'lowAzGradientFilter','minElFilter','ionosphereCoordGDZ'};
    matrix3D = {'magCartENU','magGeodeticLatLonAlt'};
    totalAddress=[address,varStr];
    if sum(strcmp(varStr,matrix2DwithTimeVars))==1
        tempInfo=h5info(h5FileStr,totalAddress);
        countIndx = tempInfo.Dataspace.Size;
        if thisTimeIndx <= countIndx(2)
            countIndx(2) = 1;
            startIndx = [1 thisTimeIndx];
        else
            error('Time index exceeds the size of records in h5 file.');
        end
    elseif sum(strcmp(varStr,array1DwithTimeVars))==1
        tempInfo=h5info(h5FileStr,totalAddress);
        countIndx = tempInfo.Dataspace.Size;
        if thisTimeIndx<=countIndx
            countIndx = 1;
            startIndx = thisTimeIndx;
        else
            error('Time index exceeds the size of records in h5 file.');
        end
    elseif sum(strcmp(varStr,matrix3DwithTimeVars))==1
        tempInfo=h5info(h5FileStr,totalAddress);
        countIndx = tempInfo.Dataspace.Size;
        if thisTimeIndx <= countIndx(3)
            countIndx(3) = 1;
            startIndx = [1 1 thisTimeIndx];
        else
            error('Time index exceeds the size of records in h5 file.');
        end
    elseif sum(strcmp(varStr,matrix4DwithTimeVars))==1
        tempInfo=h5info(h5FileStr,totalAddress);
        countIndx = tempInfo.Dataspace.Size;
        if thisTimeIndx <= countIndx(4)
            countIndx(4) = 1;
            startIndx = [1 1 1 thisTimeIndx];
        else
            error('Time index exceeds the size of records in h5 file.');
        end
    elseif sum(strcmp(varStr,matrix2D))==1
        tempInfo=h5info(h5FileStr,totalAddress);
        countIndx = tempInfo.Dataspace.Size;
        startIndx = [1 1];
    elseif sum(strcmp(varStr,matrix3D))==1
        tempInfo=h5info(h5FileStr,totalAddress);
        countIndx = tempInfo.Dataspace.Size;
        startIndx = [1 1 1];
    else
        error('Variable not found in the HDF5 file');
    end
    varValue = h5read(h5FileStr,totalAddress,startIndx,countIndx);

end