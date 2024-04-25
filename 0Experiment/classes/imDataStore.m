classdef imDataStore<handle
    %IMDATASTORE is a class for managing "imagedatastore" objects
    
    properties
        ds          (1,1) matlab.io.datastore.ImageDatastore = imageDatastore(pwd)
        time        (1,:) duration
        maxIndex    (1,1) double = 1
        filters     table
        imgroups    (:,1) cell

        mask                        % initial mask
        roi         (1,:) double    % region of interest
        rot         (1,1) double=0  % 90 degree rotation value
    end
    properties (Hidden)
        imgroup_h   double
        imgroup_s   struct          % imgroup metadata
        mask_data   struct          % roi metadata
        time0       (1,:) duration
        hasds0      logical = false
        ds0files    (1,:) cell
        ds0         (1,1) matlab.io.datastore.ImageDatastore = imageDatastore(pwd)
    end
    %% Set functions and callbacks
    methods
        function set.mask(edi,val)
            edi.mask = edi.setmask(val);
        end
        function set.roi(edi,val)
            edi.roi = edi.setroi(val);
        end
        function set.rot(edi,val)
            edi.rot = edi.setrot(val);
        end
        function set.ds(imds,val)
            if ~or(isa(val,'matlab.io.datastore.ImageDatastore'),isa(val,'imMerge'))
                error('Invalid datastore type')
            end
            imds.ds = val;
            imds.setds
        end
        function set.filters(comp,val)
           comp.filters = val;
           comp.filterds
        end
        function resetFilter(imds,val)
            if isempty(imds.filters)
                return
            end
            logic = imds.filters.logic;
            imds.filters(contains(logic,{'imtime','imgroup'}),:)=[];
            imds.ds.Files   = imds.ds0.Files;
            imds.time       = imds.time0;
            imds.imgroups   = {};
            imds.imgroup_h  = [];
            imds.imgroup_s  = [];
        end
    end
    methods (Hidden)
        function checkFiles(imds)
            imds.maxIndex   = numel(imds.ds.Files);
            imds.sortfiles
        end
    end
    methods (Access = private)
        function val = setroi(edi,val)
            try %#ok<*TRYNC> 
                if isa(val,'cell')
                    val = cell2mat(val);
                end
                if contains(char(val),{'get','reset','g','r'})
                    val = char(val);
                end
            end
            switch class(val)
                case 'double'
                    return
                case {'string','char'}
                otherwise
                    error('Not a valid argument class for "roi"')
            end

            switch val
                case 'reset'
                    val = [];
                case 'get'
                    im  = imds.ds.loadimage(edi.index,true,true);
                    val = getROI(im,edi.roi);
                    if prod(val)==numel(im)
                        val = [];
                    end
            end
        end
        function val = setrot(~,val)
            val = mod(val,4);
        end
        function val = setmask(imds,val)
            if isa(val,'logical')||isempty(val)
                return
            end
            switch val
                case 'get'
                    f = findall(0,'type','figure','name','Masking');
                    if ~isempty(f)
                        uim = findall(f.Children,'type','uiDrawMask');
                        f.focus
                    else
                        uim = uiDrawMask;
                    end
                    set(uim,...
                        'imds',imds.ds, ...
                        'mask_data',imds.mask_data);
                    uim.MaskCreatedFcn = @updateMask;
                case 'reset'
                    val = true;
            end
            function updateMask(e,~)
                imds.mask        = e.mask;
                imds.mask_data   = e.mask_data;
            end
        end
        function setds(imds)
            imds.checkFiles
            if ~imds.hasds0
                imds.ds0files   = imds.ds.Files';
                imds.ds0        = imageDatastore(unique(fileparts(imds.ds.Files)));
                imds.ds0.Files  = imds.ds0.Files(ismember(imds.ds0.Files,imds.ds.Files));
                imds.time0      = imds.time;
                imds.hasds0     = true;
            end
        end
        function sortfiles(imds)
            imds.gettime
            [~,map] = sort(imds.time);
            imds.time = imds.time(map);
            imds.ds.Files = imds.ds.Files(map);
        end
        function gettime(imds)
            d = cellfun(@dir,imds.ds.Files);
            t = datetime({d.date});
            imds.time = t-t(1);
        end
        function filterds(imds)
            if isempty(imds.filters)
                return
            end
            imds.file_filter
        end
    end
    %% Constructor
    methods 
        function imds = imDataStore(ds)
            if nargin<1
                return
            end
            imds.ds = ds;
        end
    end
    %% Image loading
    methods
        function im     = loadimage(imds,index,rotation,cropping,masking,useDS0)
            if nargin<3
                rotation = imds.rot;
            end
            if nargin<4
                cropping = imds.roi;
            end
            if nargin<5
                masking = imds.mask;
            end
            if nargin<6
                useDS0 = false;
            end
            [rotation,cropping,masking] = checkInputs(imds,rotation,cropping,masking);

            index = imds.valid_index(index);
            if numel(index)>1
                im = arrayfun(@(index) imds.loadimage(index,rotation,cropping,masking), ...
                    index,'UniformOutput',false);
                im = cat(3,im{:});
                return
            end

            im = handleMerged(imds,index);
            im = transformImage(im);

            function im=transformImage(im)
                if ~isempty(masking)
                    im(masking) = nan;
                end
                if all(cropping)
                    if ~isempty(cropping)
                        im = imcrop(im,cropping);
                    end
                end
                if rotation
                    im = rot90(im,rotation);
                end

                try %#ok<*TRYNC>
                    im = rgb2gray(im);
                end
            end
            function im     = handleMerged(imds,index)
                % Switch between "imMerge" class (two indepedendent datastores)
                % and "ImageDatastore" class

                switch class(index)
                    case 'double'
                        % limit index
                        index(index<1)=1;
                        index(index>imds.maxIndex)=[];
                        if isempty(index)
                            index = imds.maxIndex;
                        end
                    case {'char','string'}
                        index = 1:floor(imds.maxIndex/str2double(index)):imds.maxIndex;
                end

                % switch image reading mode (merged image/multistream image or
                % singular image)
                condition = isa(imds.ds,'imMerge');

                im      = switchReadMode(imds,index);

                function im = switchReadMode(imds,index)
                    datastore = imds.ds;
                    if useDS0
                        datastore = imds.ds0;
                    end
                    switch condition
                        case true
                            im = datastore.readImageAtIndex(index);
                        case false
                            im = datastore.readimage(index);
                    end

                    im = im2double(im);
                end
            end
            function [rotation,cropping,masking] = checkInputs(imds,rotation,cropping,masking)
                if isa(rotation,'logical')
                    if numel(rotation)==1 && rotation
                        rotation = imds.rot;
                    end
                end
                if isa(cropping,'logical')
                    if numel(cropping)==1 && cropping
                        cropping = imds.roi;
                    end
                end
                if isa(masking,'logical')
                    if numel(masking)==1 && masking
                        masking = imds.mask;
                    end
                end
            end
        end
      end
    methods (Access = protected)
        function index  = valid_index(imds,index)
            index(index<1)=[];
            if isempty(index)
                index = 1;
            end
            index(index>imds.maxIndex)=[];
            if isempty(index)
                index = imds.maxIndex;
            end

        end
    end
    
    %% Image management
    methods
        function [idx,h] = image_identify(imds,index,method,threshold,silent)
            %%IMAGE_IDENTIFY
            if nargin<3
                method = 'hist';
            end
            if nargin<4
                threshold = 0.2;
            end
            if nargin<5
                silent = false;
            end
            if numel(index)==1
                index = 1:index:imds.maxIndex;
            end
            index = imds.valid_index(index);

            state = isReScan(imds,method,index);
            imds.imgroup_s  = struct('method',method,'threshold',threshold,'index',index);

            if ~state
                h = arrayfun(@(index) idImage(index,imds,method),index,'UniformOutput',false);
                h = cat(2,h{:});
                switch method
                    case 'hist'
                        [~,i]   = sortrows(h');
                        h       = sum(h(:,i)./mean(h,2),'omitnan');
                    case 'sum'
                        [~,i] = sort(h);
                        h = h(i);
                    case 'fft'
                        [h,i] = sort(sum(h./mean(h,2),'omitnan'));
                end
                imds.imgroup_s.index = index(i);
            else
                h = imds.imgroup_h;
            end
            index = imds.imgroup_s.index;
            

            % find groups
            h0 = h;
            h = rescaleMax(abs(gradient(h)));
            idx = double(h>threshold);
            idx(1) = 1;
            idx = find(gradient(idx)>0);

            if ~isempty(idx)
                idx = [idx(1) diff(idx)];
                idx = cumsum(idx);
                if numel(idx)>1
                    idx = [1 idx([idx(1) diff(idx)]>1)];
                    idx = [idx;circshift(idx,-1)];
                    idx(1,2:end) = idx(1,2:end)+1;
                    idx(end) = numel(h);
                else
                    idx = [1 idx+1;idx numel(h)];
                end

                idx = arrayfun(@(j) index(idx(1,j):idx(2,j)),1:width(idx),'UniformOutput',false);
                idx = cellfun(@sort,idx,'UniformOutput',false);
            else
                idx = {1:numel(h)};
            end

            imds.imgroups   = idx;
            imds.imgroup_h  = h0;

            if ~nargout
                clear idx h
            end

            function state = isReScan(imds,method,index)
                state = false;
                gs = imds.imgroup_s;
                if isempty(gs)
                    return
                end
                if ~strcmp(gs.method,method)
                    return
                end
                if all(~ismember(gs.index,index))
                    return
                end
                state = true;


            end
            function h = idImage(i,imds,method)
                if ~silent
                    displayProgress('Identifying',i,index(1),index(end))
                end
                im = imds.loadimage(i);%,false,false,false,true);
                switch method
                    case 'hist'
                        h = imhist(im,100);
                    case 'sum'
                        h = sum(im,[1 2]);
                    case 'fft'
                        h = reshape(fftshow(fft2(im)),[],1);
                end

            end
        end
    end
    
    %% File management (copying, moving, renaming/labelling)
    methods
        function file_filter(imds,filter)
            if nargin<2
                filter = imds.filters;
            end
            if isempty(filter)
                return
            end

            files   = imds.ds0.Files;
            grp     = contains(filter.logic,{'imgroup','imtime'});
            imgrp   = filter(grp,:).string(:);
            files   = filter_time(files,imgrp,filter);
            files   = filter_group(files,imgrp,filter);
            files   = filter_string(files,grp,filter);
            imds.ds.Files = files;
            imds.checkFiles
            
            function files = filter_time(files,grp,filter)
                valid = contains(filter.logic,'imtime');
                if all(~valid)
                    return
                end
                files = files(grp{valid});
            end
            function files = filter_group(files,grp,filter)
                valid = contains(filter.logic,'imgroup');
                if all(~valid)
                    return
                end
                files = files(grp{contains(filter.logic,'imgroup')});
            end
            function files = filter_string(files,grp,filter)
                if isempty(filter(~grp,:))
                    return
                end
                filter  = filter(~grp,:);
                f = filter(~grp,:);
                f = sortrows(f,'logic');
                val = cellfun(@(str,logic) contains(files,str),f.string,'UniformOutput',false);
                val = cat(2,val{:});
                val = val.*~contains(f.logic,'not')';
                lgc = val(:,1);
                for i=2:width(val)-1
                    lgc = parseLogic(lgc,val(:,i),f.logic{i});
                end
                files = files(logical(lgc));

                function lgc = parseLogic(array1,array2,logic) %#ok<*INUSD>
                    if strcmp(logic,'not')
                        logic = 'and';
                    end
                    lgc = eval(sprintf('%s(array1,array2);',logic));
                end
            end
        end
        function file_copy(imds,from,to)
        end
        function file_filterApp(imds)
            uiims = uiImageDataStoreSettings('imds',imds);
        end
    end
end

