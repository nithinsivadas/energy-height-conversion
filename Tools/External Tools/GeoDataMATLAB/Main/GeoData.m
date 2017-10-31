classdef GeoData <matlab.mixin.Copyable%handle
    % Class that holds data for different sensors.
    % This is a class to hold geophysical data from different types of 
    % sensors to study near earth space physics.
    properties
        data% struct
        coordnames % type of coordinates.
        dataloc % location of data points
        sensorloc% sensor location in lla
        times% times in posix formats
    end
    
    methods
        %% Init
        function self = GeoData(readmethod,varargin)
            % This function will be the contructor for the GeoData class.
            % The two inputs are a file handle and the set of inputs for
            % the file handle. The outputs must follow the output structure
            % below.
            [self.data,self.coordnames,self.dataloc,self.sensorloc,self.times] = readmethod(varargin{:});
            self.timefix()
        end
        %% == ~=
        function out = eq(self,GD2)
            % This is the == operorator for the GeoData class 
            % Checks to see if all of the locations, times, parameters
            % and data sets are the same.
            proplist = properties(self);
            proplist2 = properties(GD2);
            
            if ~isequal(proplist,proplist2)
                error('GD and/or GD2 is not a GeoData object');
            end
            
            for k = 1:length(proplist)
                prop1 = self.(proplist{k});
                prop2 = GD2.(proplist{k});
                if ~isequaln(prop1,prop2)
                    out=false;
                    return
                end
            end
            out=true;
        end
        
        function out=ne(self,GD2)
            % This is the ~= operorator for the GeoData class
            % Checks to see if all of the locations, times, parameters
            % and data sets are not the same.
            out = ~(self==GD2);
        end
        %% datanames
        function dnames = datanames(self)
            % This will output a cell array of strings which hold the data
            % names.
            dnames = fieldnames(self.data);
        end
        %% Change Data
        function changedata(self,dataname,newname,func,varargin)
            % This will change the data points based on a given function.
            % changedata(dataname,newname,func)
            % changedata(dataname,newname,func,params)
            % changedata(dataname,newname,func,params,rm_old)
            %
            % This function will take one of the elements in the data
            % matrix and apply the transform in function handle given. It
            % will then assign a new name in the struct.
            % Inputs
            % dataname: a string, the name of the field you want to operate
            % on.
            % newname: a string that will be name of the resultant field.
            % func: a function handle that will augment the data. The
            % first argument must be the data from the instance of the
            % class.
            % params: (default ={}) A list of parameters for the function.
            % rm_old: (default= true) a bool to determine if you should
            % remove the old field.
            if nargin==4
                params={};
                rm_old=true;
            elseif nargin==5
                params=varargin{1};
                rm_old=true;
            elseif nargin==6
                params=varargin{1};
                rm_old=varargin{2}; 
            end
            
            self.data.(newname)= func(self.data.(dataname),params{:});
            if rm_old
                self.data=rmfield(self.data,dataname);
            end
            
        end
        %% Reduce time
        function  timereduce(self,timelist,varargin)
            % This pull out time instances.
            % Give a list of time instances that the user wants to keep and
            % discard the rest.
            % Inputs
            % timelist - A vector of desired times
            % type - Either 1 or 2, if 1 then just give the list of indices
            % for the times. If 2 then give the actual time values.
            if nargin==2
                type = 1;
            elseif nargin==3
                type=varargin{1};
            end
            assert(any(type==[1,2]),'Type must be 1 or 2');
            
            if type ==1
                listkeep = timelist;
            elseif type==2
                [~,~,listkeep] = intersect(self.times,timelist);
            end
            if ndims(self.times)==1
                self.times=self.times(listkeep);
            elseif ismatrix(self.times)
                self.times=self.times(listkeep,:);
            end
            
            datafields = self.datanames();
            
            for ifield = 1:length(datafields)
                self.data.(datafields{ifield})= self.data.(datafields{ifield})(:,listkeep);
            end
            
        end
        %% Interpolate
        function interpolate(self,new_coords,newcoordname,varargin)
            % This will interpolate the data in a GeoData object.
            % Uses standard coordinate change routines avalible in the tool
            % box for GeoData.
            % Inputs
            % new_coords - A Nx3 array of coordinates to interpolate the
            % data to.
            % newcoordname - A string for the new type of coordinates to
            % interpolate to.
            
            % method - A string of 'linear', 'nearest', 'cubic','natural'
            % to determine how to interpolate the data. The default is
            % 'natural'.
            % extrapmethod - A string that determines the extrapolation 
            % method. Default is none, which fills extrapolated values with nans. 
            if nargin <4
                method = 'linear';
                extrapmethod = 'none';
                twodinterp = false;
            elseif nargin <5
                method = varargin{1};
                extrapmethod = 'none';
                twodinterp = false;
            elseif nargin <6
                method = varargin{1};
                extrapmethod = varargin{2};
                twodinterp = false;
            else
                method = varargin{1};
                extrapmethod = varargin{2};
                twodinterp = varargin{3};
            end
            
            curavalmethods = {'linear', 'nearest', 'cubic','natural'};
            interpmethods = {'linear', 'nearest', 'cubic','natural'};
            assert(any(strcmp(curavalmethods,method)),...
                ['Must be one of the following methods: ', strjoin(curavalmethods,', ')]);
            
            Nt = length(self.times);
            ONlocs = size(self.dataloc,1);
            NNlocs = size(new_coords,1);
            curcoords = self.changecoords(newcoordname);
            
            origloc = self.dataloc(1,:);
            loclog = origloc(ones(ONlocs,1),:)==self.dataloc;
            dimrm = all(loclog,1);
            
            origloc2 = new_coords(1,:);
            loclog2 = origloc2(ones(NNlocs,1),:)==new_coords;
            dimrm2 = all(loclog2,1);
            
            dimrm = dimrm|dimrm2;
            newcoordshold = new_coords;
            if twodinterp
                new_coords = new_coords(:,~dimrm);
                curcoords = curcoords(:,~dimrm);
            end
            
            
            % Loop through parameters and create temp variable
            paramnames = fieldnames(self.data);
            for iparam =1:length(paramnames)
                fprintf('Interpolating parameter %s, %d of %d\n',paramnames{iparam},iparam,length(paramnames));
                iparamn = paramnames{iparam};
                New_param = zeros(NNlocs,Nt);
                for itime = 1:Nt
                    fprintf('\tInterpolating time %d of %d\n',itime,Nt);
                    curparam =self.data.(iparamn)(:,itime);
                    paramnum = ~isnan(curparam);
                    curparam = curparam(paramnum);
                    curcoordstemp = curcoords(paramnum,:);
                    if any(strcmp(method,interpmethods))
                        F = scatteredInterpolant(curcoordstemp,curparam,method,extrapmethod);
                        intparam = F(new_coords);
                        if isempty(intparam)
                            New_param = nan(NNlocs,1);
                        else
                            New_param(:,itime)= intparam;
                        end
                    end
                end
                self.data.(iparamn) = New_param;
                
            end
            self.coordnames=newcoordname;
            self.dataloc = newcoordshold;       
        end
        
        function interpolate_customcoords(self,new_coords,newcoordname,curcoords,varargin)
            % This will interpolate the data in a GeoData object given a
            % custom set of coordinates of the original data.
            % Inputs
            % new_coords - A Nx3 array of coordinates to interpolate the
            % data to.
            % newcoordname - A string for the new type of coordinates to
            % interpolate to.
            % curcoords - An alternative set of coordinates to for the 
            % original data to be used instead of its current coordinates.
            % method - A string of 'linear', 'nearest', 'cubic','natural'
            % to determine how to interpolate the data. The default is
            % 'natural'.
            % extrapmethod - A string that determines the extrapolation 
            % method. Default is none, which fills extrapolated values with nans. 
            if nargin <5
                method = 'linear';
                extrapmethod = 'none';
            elseif nargin <6
                method = varargin{1};
                extrapmethod = 'none';
            else
                method = varargin{1};
                extrapmethod = varargin{2};
            end
            
            curavalmethods = {'linear', 'nearest', 'cubic','natural'};
            interpmethods = {'linear', 'nearest', 'cubic','natural'};
            assert(any(strcmp(curavalmethods,method)),...
                ['Must be one of the following methods: ', strjoin(curavalmethods,', ')]);
            
            Nt = length(self.times);
            ONlocs = size(self.dataloc,1);
            NNlocs = size(new_coords,1);
            % Loop through parameters and create temp variable
            paramnames = fieldnames(self.data);
            for iparam =1:length(paramnames)
                fprintf('Interpolating parameter %s, %d of %d\n',paramnames{iparam},iparam,length(paramnames));
                iparamn = paramnames{iparam};
                New_param = zeros(NNlocs,Nt);
                for itime = 1:Nt
                    fprintf('\tInterpolating time %d of %d\n',itime,Nt);
                    curparam =self.data.(iparamn)(:,itime);
                    paramnum = ~isnan(curparam);
                    curparam = curparam(paramnum);
                    curcoordstemp = curcoords(paramnum,:);
                    if any(strcmp(method,interpmethods))
                        F = scatteredInterpolant(curcoordstemp,curparam,method,extrapmethod);
                        intparam = F(new_coords);
                        if isempty(intparam)
                            New_param = nan(NNlocs,1);
                        else
                            New_param(:,itime)= intparam;
                        end
                    end
                end
                self.data.(iparamn) = New_param;
                
            end
            self.coordnames=newcoordname;
            self.dataloc = new_coords;
            
        end
        %% Time registration
        function outcell = timeregister(self,self2)
            % Create a cell array which shows the overlap between two
            % instances of GeoData.
            % Inputs
            % self2 - A GeoData object.
            % Outputs
            % outcell - A cellarray of vectors the same length as the time
            % vector in self. Each vector will have the time indecies from 
            % the second GeoData object which overlap with the time indicie 
            % of the first object.
            times1 = self.times;
            times2 = self2.times;
            
            if ndims(times1)==1
                timeahead = times1(2:end);
                timeahead = timeahead(:);
                avediff = mean(diff(timeahead));
                timesahead = [timeasahead;timeahead(end)+avediff];
                times1 = [times1(:),timesahead];
            end
            
            if ndims(times2)==1
                timeahead = times2(2:end);
                timeahead = timeahead(:);
                avediff = mean(diff(timeahead));
                timesahead = [timeasahead;timeahead(end)+avediff];
                times2 = [times2(:),timesahead];
            end
            outcell = cell(1, size(times1,1));
            for k =  1:size(times1,1)
                l = times1(k,:);
                
                ind1 = find(l(1)>times2(:,1),1,'last');
                ind2 = find(l(2)<times2(:,2),1,'first');
                outcell{k}=ind1:ind2;
            end
        end
        function timefix(self)
            if size(self.times,2)==2
                return
            end
            
            time1 = self.times(:);
            avdiff = mean(diff(time1));
            time2 = circshift(time1,-1);
            time2(end) = time2(end-1)+avdiff;
            self.times = [time1,time2];
        end
        %% Coordinate change
        function oc = changecoords(self,newcoordname)
            % Changes the coordinates
            d2r = pi/180;
            cc = self.dataloc;
            if strcmpi(self.coordnames,'spherical')&&strcmpi(newcoordname,'enu')
                [x,y,z] = sphere2cart(cc(:,1),cc(:,2),cc(:,3));
                oc = [x,y,z];
            elseif strcmpi(self.coordnames,'spherical')&&strcmpi(newcoordname,'cartesian')
                [x,y,z] = sphere2cart(cc(:,1),d2r*cc(:,2),d2r*cc(:,3));
                oc = [x,y,z];
            elseif strcmpi(self.coordnames,'enu')&&strcmpi(newcoordname,'spherical')
                [r,az,el] = cart2sphere(cc(:,1),cc(:,2),cc(:,3));
                oc = [r,az,el];
            elseif strcmpi(self.coordnames,'wgs84')&&strcmpi(newcoordname,'enu')
                ECEF_COORDS = wgs2ecef(self.dataloc.');
                locmat = repmat(self.sensorloc',[1,size(ECEF_COORDS,2)]);
                oc = ecef2enul(ECEF_COORDS,locmat);
                oc = oc.';
            elseif strcmpi(self.coordnames,'enu')&&strcmpi(newcoordname,'cartesian')
                oc = cc*1e-3;
            elseif strcmpi(self.coordnames,'cartesian')&&strcmpi(newcoordname,'enu')
                oc = cc*1e3;
            elseif strcmpi(self.coordnames,newcoordname)
                oc = cc;
            end
        end
        %% Write out
        function write_h5(self,filename)
            % This will write out the h5 file in our defined format.
            % Input filename must be a string.
            proplist = properties(self);
%             file_id  = H5F.create(filename, 'H5F_ACC_TRUNC', ...
%                              'H5P_DEFAULT', 'H5P_DEFAULT');
%             H5F.close(file_id);            
            for k = 1:length(proplist)
                prop1 = self.(proplist{k});
                if isa(prop1,'struct')
                    fnames = fieldnames(prop1);
                    for l = 1:length(fnames)
                        value = prop1.(fnames{l});
                        rvalue = permute(value,ndims(value):-1:1);
                        location = ['/',proplist{k},'/',fnames{l}];
                        h5create(filename,location,size(rvalue));
                        h5write(filename,location,rvalue)
                    end
                else
                    location = ['/',proplist{k}];
                    % TODO make this into a seperate function
                    % For some god damn reason matlab can not write strings
                    % to HDF files so for now we have this bull shit.
                    if ischar(prop1)
%                        
                        file_id = H5F.open(filename,'H5F_ACC_RDWR','H5P_DEFAULT');
                        space_id = H5S.create('H5S_SCALAR');
                        stype = H5T.copy('H5T_C_S1'); 
                        sz = numel(prop1);  
                        H5T.set_size(stype,sz);
                        
                        dataset_id = H5D.create(file_id,proplist{k}, ...
                            stype,space_id,'H5P_DEFAULT');
                        H5D.write(dataset_id,stype,'H5S_ALL','H5S_ALL','H5P_DEFAULT',prop1);
                        H5D.close(dataset_id)
                        H5S.close(space_id)
                        H5F.close(file_id);
                    else % for none char values.
                        rprop = permute(prop1,ndims(prop1):-1:1);
                        h5create(filename,location,size(rprop));
                        h5write(filename,location,rprop);
                    end
                end
            end
        end
    end
    
end

