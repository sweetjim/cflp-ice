classdef edIce<handle& matlab.mixin.SetGet
    %EDICE 'edge-detection: ice' 
    %% Global properties
    properties
        %IDS matlab.io.datastore.ImageDatastore = imageDatastore(pwd)
        imds            imDataStore = imDataStore(imageDatastore(pwd))
        maxIndex        (1,1) double
        labExperiment   labExperiment
    end
    %% Edge-detection properties
    properties
        ds          (1,1) double = 4.2700e-04; % pixels to meters
        mask                        % initial mask
        roi         (1,:) double    % region of interest
        rot         (1,1) double=0  % 90 degree rotation value
        sequence    table           % sequence for edge-detection algorithm with filter info
        
        % algorithm filters blocks  | elements
        s_image     struct          % vignette, clim, blur, image difference
        s_edge      struct          % edge-detection, gradient threshold
        s_dilation  struct          % dilation, median filtering, pixel rejection
        s_conn      struct          % connections, smoothing
        
        % algorithm output
        a_output    struct
    end
    properties (Hidden)
        time        duration        % file time 
        s_all       struct          % all algorithm filter blocks
        savepath    {mustBeFolder(savepath)} = pwd % default is root of image directory
        savelabel   string          % label for saved files (<savepath>\ediData_<label>.mat)
        index       (1,1) double=1  % image index
        roirot      (1,1) double=0  % roi set by rot
        d_output    (:,:)           % floating filter output
        h_output    (:,:) double    % floating detected edge output
        h_ref       (:,:) double    % floating reference detected edge output
        mask_data   struct          % roi metadata
    end
    %% Set functions
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
        function set.imds(edi,val)
            edi.imds = val;
            edi.setimds
        end
        function set.s_all(edi,val)
            if isempty(val)
                val = get_s_all(edi);
            end
            edi.set_s_all(val)
        end
    end
    %% Get functions
    methods (Hidden)
        function out = get_s_all(edi)
            S_all.s_image       = edi.s_image;
            S_all.s_edge        = edi.s_edge;
            S_all.s_dilation    = edi.s_dilation;
            S_all.s_conn        = edi.s_conn;
            edi.s_all   = S_all;
            out         = S_all;
        end
    end
    methods (Access = private)
        function out = get_output(edi,index)
            if nargin<2
                index = edi.index;
            end
            out = edi.d_output;
            if isempty(out)
                out = edi.loadimage(index);
            end
        end
    end
    %% Set-callbacks
    methods (Access = private)
        function val = setmask(edi,val)
            if isa(val,'logical')
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
                        'imds',edi.imds.ds, ...
                        'index',edi.index,...
                        'mask_data',edi.mask_data);
                    uim.MaskCreatedFcn = @updateMask;
                case 'reset'
                    val = true;
            end
            edi.updateimds
            function updateMask(e,~)
                edi.mask        = e.mask;
                edi.mask_data   = e.mask_data;
            end
        end
        function set_s_all(edi,val)
            edi.s_image     = val.s_image;
            edi.s_edge      = val.s_edge;
            edi.s_dilation  = val.s_dilation;
            edi.s_conn      = val.s_conn;
        end
        function val = setrot(~,val)
            val = mod(val,4);
            edi.updateimds
        end
        function rotateroi(edi)
            if isempty(edi.roi)
                edi.roirot = edi.rot;
                return
            end
            if edi.rot==edi.roirot
                return
            end
            im  = edi.loadimage(1,true,true);
            dr  = edi.rot-edi.roirot;

            % create mask
            m   = zeros(size(im));
            m(edi.roi(2):sum(edi.roi([2 4])),edi.roi(1):sum(edi.roi([1 3])))=1;

            % rotate mask to 0-rot and find new roi via bounding box
            n       = rot90(m,dr);
            p       = regionprops(n);
            edi.roi = p.BoundingBox;
            
            edi.roirot = edi.rot;
        end
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
                    im = edi.loadimage(edi.index,true,true);
                    val = getROI(im,edi.roi);
                    if prod(val)==numel(im)
                        val = [];
                    end
            end
            edi.roirot = edi.rot;
            edi.updateimds
        end
        function setimds(edi)
            edi.imds.maxIndex  = edi.imds.maxIndex;
            edi.setsavepath
            edi.sortfiles
            edi.updateimds
        end
        function updateimds(edi)
            edi.imds.roi = edi.roi;
            edi.imds.rot = edi.rot;
            edi.imds.mask = edi.mask;
            edi.imds.mask_data = edi.mask_data;
        end
        function sortfiles(edi)
            edi.gettime
            [~,map] = sort(edi.time);
            edi.time = edi.time(map);
            edi.imds.ds.Files = edi.imds.ds.Files(map);
        end
        function gettime(edi)
            d = cellfun(@dir,edi.imds.ds.Files);
            t = datetime({d.date});
            edi.time = t-t(1);
        end
        function setsavepath(edi)
            if ~isempty(edi.imds.ds.Folders)
                [r,lbl] = fileparts(edi.imds.ds.Folders);
            else
                r = unique(fileparts(edi.imds.ds.Files));
                [~,lbl] = fileparts(r);
            end
            if numel(r)>1
                r   = unique(fileparts(r));
                lbl = '+';
            end
            if isa(r,"cell")
                r = r{:};
            end
            edi.savepath    = string(r);
            edi.savelabel   = string(lbl);
        end
    end
    %% Default filter blocks and sequence
    methods
        function set_default_sequence(edi)
            blocks = {'image','edge','dilation','conn'};
            state  = false(size(blocks));
            edi.sequence = array2table(state,'VariableNames',blocks,'RowNames',{'state'});
        end
        function set_default_fblocks(edi)
            edi.s_image         = edi.get_default_fblocks('image');
            edi.s_edge          = edi.get_default_fblocks('edge');
            edi.s_dilation      = edi.get_default_fblocks('dilation');
            edi.s_conn          = edi.get_default_fblocks('conn');
        end
        function s_item = get_default_fblocks(~,item)
            
            S_image.vignette    = 0;        % double    [0,100] flatfield correction
            S_image.imdiff      = 0;        % int       [0,inf] index separation for image difference
            S_image.clim        = [];       % double    [0,1]   
            S_image.blurx       = 0;        % int       [0,inf] horizontal blur strength
            S_image.blury       = 0;        % int       [0,inf] vertical blur strength
            
            
            S_edge.gradient     = 0;        % double    [0,inf] gv/px threshold
            S_edge.algorithm    = 'canny';  % char      {'sobel','prewitt','roberts','canny','log'} edge-detection algorithm modifier
            S_edge.threshold    = 0.05;     % double    [0,inf] edge-detection threshold strength
            S_edge.sigma        = 1;        % double    [0,1]   edge-detection sigma strength
            

            S_dilation.shape    = 'cross';  % char      {'line','rect','cross','circ'} dilation shape modifier
            S_dilation.arg1     = 0;        % int       [0,inf] length argument 
            S_dilation.arg2     = 0;        % int       [0,inf] rotation argument
            S_dilation.medfilt  = 0;        % double    [0,inf] median filter strength
            S_dilation.reject   = 0;        % int       [0,inf] clumped-pixels rejection threshold
            

            S_conn.connect      = 1;        % int       [0,inf] number of longest connections threshold
            S_conn.smooth       = 0;        % int       [0,inf] smoothing operation strength
            S_conn.ref          = 0;        % int       [0,maxIndex] reference edge index (0 negates operation)
            S_conn.direction    = 'right';  % char      {'left','right','up','down'} direction to take longest connection
            

            try
                items = split(item,'.');
                s_item = eval(sprintf('S_%s.%s',items{1},items{2}));
                return
            end
            try
                switch item
                    case 'image'
                        s_item = S_image;
                    case 'edge'
                        s_item = S_edge;
                    case 'dilation'
                        s_item = S_dilation;
                    case {'conn','connectivity','connect'}
                        s_item = S_conn;
                end
            end
        end
        function s_item = get_default_ranges(edi,item)
            %GET_DEFAULT_RANGES is used for component creation
            if nargin<2
                item = [];
            end
            S_image.vignette    = [0,100];        % double    [0,100] flatfield correction
            S_image.imdiff      = [0,inf];        % int       [0,inf] index separation for image difference
            S_image.clim        = [0,1];         % double    [0,255] gray-value (gv) limits 256bit
            S_image.blurx       = [0,10];        % int       [0,inf] horizontal blur strength
            S_image.blury       = [0,10];        % int       [0,inf] vertical blur strength
            
            S_edge.gradient     = [0,1]*1e-2;        % double    [0,inf] gv/px threshold
            S_edge.algorithm    = {'sobel','prewitt','roberts','canny','log'};% edge-detection algorithm modifier
            S_edge.threshold    = [0,1];     % double    [0,inf] edge-detection threshold strength
            S_edge.sigma        = [0,1];        % double    [0,1]   edge-detection sigma strength

            S_dilation.shape    = {'line','rect','cross','circ'};% dilation shape modifier
            S_dilation.arg1     = [0,10];        % int       [0,inf] length argument 
            S_dilation.arg2     = [0,10];        % int       [0,inf] rotation argument
            S_dilation.medfilt  = [0,100];        % double    [0,inf] median filter strength
            S_dilation.reject   = [0,1e5];        % int       [0,inf] clumped-pixels rejection threshold

            S_conn.connect      = [0,10];        % int       [0,inf] number of longest connections threshold
            S_conn.smooth       = [0,5e2];        % int       [0,inf] smoothing operation strength
            S_conn.ref          = [0,edi.imds.maxIndex];        % int       [0,maxIndex] reference edge index (0 negates operation)
            S_conn.direction    = {'left','right','up','down'};
           
           
            try
                items = split(item,'.');
                s_item = eval(sprintf('S_%s.%s',items{1},items{2}));
                return
            end
            try
                switch item
                    case 'image'
                        s_item = S_image;
                    case 'edge'
                        s_item = S_edge;
                    case 'dilation'
                        s_item = S_dilation;
                    case {'conn','connectivity','connect'}
                        s_item = S_conn;
                end

            end
        end
    end
    %% CONTSTRUCTOR
    methods
        function edi = edIce(imdatastore)
            if nargin<1
                imdatastore = pwd;
            end
            edi.imds        = imdatastore;
            edi.set_default_fblocks
            edi.set_default_sequence
        end
        function im  = loadimage(edi,index,norot,nocrop)
            if nargin<3
                norot = false;
            end
            if nargin<4
                nocrop = false;
            end
            index = edi.valid_index(index);
            if numel(index)>1
                im = arrayfun(@(index) edi.loadimage(index),index,'UniformOutput',false);
                im = cat(3,im{:});
                return
            end
            im = readimage(edi.imds.ds,index);
            im = im2double(im);
            if ~isempty(edi.mask)
                im(edi.mask) = nan;
            end
            if ~nocrop
                if ~isempty(edi.roi)
                    im = imcrop(im,edi.roi);
                end
            end
            if ~norot
                im = rot90(im,edi.rot);
            end
           
            try
                im = rgb2gray(im);
            end
        end
    end
    %% SAVING
    methods
        function saveData(edi)
            save( ...
                fullfile(edi.savepath,join(['ediData' edi.savelabel],'_')), ...
                "edi")
        end
    end
    %% METHODS Edge-detection sequencer
    methods
        function [out,h,z,t]    = parse_ed(edi,range,varargin)

            p = inputParser;
            addOptional(p,'silent',false)
            addOptional(p,'parallel',false)
            parse(p,varargin{:})

            track = ~p.Results.silent;
            useParallel = p.Results.parallel;
           
            if nargin<2
                range = 1:edi.imds.maxIndex;
            end
            if numel(range)==1
                if ~isinf(range)
                    range = 1:range:edi.imds.maxIndex;
                else
                    range = [1 edi.imds.maxIndex];
                end
            end
           
            range = edi.valid_index(range);

            if ~isempty(gcp('nocreate'))&&useParallel
                q = parallel.pool.DataQueue;
                parfevalOnAll(@warning, 0,'off','all');
                afterEach(q,@displayProgressparfor)
                displayProgressparfor([],numel(range))

                h = edi.parse_sequence(1);
                h = nan(numel(h),numel(range));
                parfor i=1:numel(range)
                    h(:,i) = edi.parse_sequence(range(i)); %#ok<PFBNS> 
                    send(q,[])
                end
            else
                h = arrayfun(@(i) pFcn(edi,i),range,'UniformOutput',false);
                if any(cellfun(@isrow,h))
                    h = cat(1,h{:})';
                else
                    h = cat(2,h{:});
                end
            end

            [z,t]   = edi.getDimensions;
            t       = t(range);
            s_out.h = h*edi.ds;
            s_out.z = z;
            s_out.t = t;
            out = s_out;

            edi.a_output = s_out;
            function h=pFcn(edi,i)
                h=edi.parse_sequence(i);
                if track
                    displayProgress('Processing',i,range(1),range(end))
                end
            end
        end
        function h      = parse_sequence(edi,index,overlay,ax)
            if nargin<2
                index = 1;
            end
            if nargin<3
                overlay = false;
            end
            if nargin<4
                if nargout
                    ax = gca;
                end
            end
            s = edi.sequence;
            if (all(not([s.edge,s.dilation,s.dilation])))
                overlay = false;
            end
            edi.index = edi.valid_index(index);
            im  = edi.loadimage(index);
            edi.d_output = im;
            h = edi.get_sequence;

            if ~nargout
                updateAxes(ax)
                clear h
                return
            end

            if isempty(h)
                h = nan(1,size(im,1));
            end

            function updateAxes(ax)
                imd = findall(ax,'type','image');
                pd  = findall(ax,'type','line','tag','h');
                if isempty(imd)
                    if overlay
                        if ~edi.sequence.conn
                            imd = imagesc(ax,im);
                            try
                                imd.AlphaData = ~edi.d_output;
                            end
                        else
                            imagesc(ax,im)
                        end
                    else
                        imagesc(ax,edi.d_output)
                    end
                else

                    imd.AlphaData = true;
                    if overlay
                        imd.CData = im;
                        if ~edi.sequence.conn
                            try
                                imd.AlphaData = ~edi.d_output;
                            end
                        end
                    else
                        imd.CData  = edi.d_output;
                    end
                end
                if isempty(pd)
                    hold(ax,'on')
                    try
                        plot(ax,h,1:size(im,1),'r','tag','h')
                    end
                    hold(ax,'off')
                else
                    
                    pd.XData = h;
                    if isempty(h)
                        pd.YData = [];
                    else
                        pd.YData = 1:size(im,1);
                    end
                end

                 set(ax,'Color','w')
                if overlay
                    set(ax,'Color','b')
                end
            end
        end
    end
    methods (Access = private)
        function h = get_sequence(edi)
            h = [];
            s = edi.sequence;
            if s.image;edi.filter_image;end
            if s.edge;edi.filter_edge;end
            if s.dilation;edi.filter_dilation;end
            if s.conn;h=edi.filter_conn;end
        end
    end
    %% METHODS Edge-detection filters functions
    methods
        function out = filter_image(edi,in)
            filters = edi.s_image;
            if nargin<2
                in = get_output(edi);
            end
            out = f_imdiff(in);
            out = f_vignette(out);
            out = f_clim(out);
            out = f_blur(out);
            edi.d_output = out;

            function out = f_imdiff(in)
                out = in;
                if ~filters.imdiff
                    return
                end
                try
                    idx     = edi.index;
                    imd     = filters.imdiff;
                    r       = [idx+imd idx idx-imd];
                    r       = r(1:2);

                    im          = edi.loadimage(r);
                    [~,~,dtdi]  = gradient(im);
                    out         = mean(dtdi,3,'omitnan');
                end
            end
            function out = f_vignette(in)
                out = in;
                try
                    map = isnan(in);
                    in(map)=mean(in,"all",'omitnan');
                    out = imflatfield(in,filters.vignette);
                    out(map)=nan;
                end
            end
            function out = f_clim(in)
                out = in;
                try
                    lims    = filters.clim;
                    if isinf(lims(2))
                        lims(2) = max(in,[],'all');
                    end
                    imLow   = in>lims(1);
                    imHigh  = in<lims(2);
                    imRange = imLow.*imHigh;
                    out      = imRange.*in + ~imLow*lims(1) + ~imHigh*lims(2);
                end
            end
            function out = f_blur(in)
                out = in;
                try
                    map = isnan(in);
                    in(map)=mean(in,"all");
                    if filters.blurx==0||filters.blury==0
                        out = in;
                        return
                    else
                        out = imgaussfilt(in,[filters.blurx filters.blury]);
                    end
                    out(map)=nan;
                end
            end
        end
        function out = filter_edge(edi,in)
            if nargin<2
                in = get_output(edi);
            end
            filters = edi.s_edge;
            out = f_gthres(in);
            out = f_edgedetect(out);
            edi.d_output = out;
            function out = f_gthres(in)
                out = in;
                try
                    if filters.gradient>0
                        in  = gradient(in);
                        out = in>max(in,[],'all')*filters.gradient*1e-2;
                    end
                end
            end
            function out = f_edgedetect(in)
                out = in;
                try
                    switch lower(filters.algorithm)
                        case {'sobel','prewitt','roberts'}
                            out = edge(in, ...
                                filters.algorithm, ...
                                filters.threshold);
                        case {'canny','log'}
                            if filters.threshold>=1
                                filters.threshold=1-1e-3;
                            end
                            out = edge(in, ...
                                filters.algorithm, ...
                                filters.threshold, ...
                                filters.sigma);
                    end
                end
            end
        end
        function out = filter_dilation(edi,in)
            if nargin<2
                in = get_output(edi);
            end
            filters = edi.s_dilation;
            out     = f_dilate(in);
            out     = f_medfilt(out);
            out     = f_reject(out);

            edi.d_output = out;

            function out = f_dilate(in)
                out = in;
                try
                    if all([filters.arg1 filters.arg2])
                        dilation = struct('factor',[filters.arg1 filters.arg2]);
                        if numel(dilation.factor)==1
                            line1 = strel('line',dilation.factor,45);
                            line2 = strel('line',dilation.factor,-45);
                            dilate_matrix = imbinarize(line1.Neighborhood+line2.Neighborhood);
                        else
                            dilate_matrix = dilateMatrix(filters.shape,...
                                dilation.factor(1),dilation.factor(2));
                        end

                        if numel(dilate_matrix)<9
                            dilate_matrix = true(dilation.factor);
                        end
                        out = imdilate(in,dilate_matrix);
                        edi.s_dilation.dilationMatrix = dilate_matrix;
                    end
                end
            end
            function out = f_medfilt(in)
                out = in;
                try
                    skip = false;
                    if filters.medfilt==0
                        out = gather(in);
                        skip = true;
                    end

                    if ~skip
                        if ...
                                within(filters.medfilt,[3 15])...       % GPU supports only odd-length sides from 3-15
                                &&(edi.useGPU)...                       % GPU is active
                                &&(mod(filters.medfilt,2)~=0)           % odd
                            out = gather(medfilt2(in,[1 1].*filters.medfilt));
                        else % convert to CPU
                            out = medfilt2(gather(in),[1 1].*filters.medfilt);
                        end
                    end
                end
            end
            function out = f_reject(in)
                out = in;
                try
                    out     = bwareaopen(in,round(filters.reject),4);
                end
            end
        end
        function out = filter_conn(edi,in)
            filters = edi.s_conn;
            if nargin<2
                in = edi.get_output;
            end
            out = f_con(in);
            out = f_smooth(out);
            out = f_ref(out);
            out(out<0) = nan;
            edi.h_output = out;

            function out = f_con(in)
                out = [];
                try
                    biggest_n   = filters.connect;
                    CC          = bwconncomp(in,8);
                    numPixels   = cellfun(@numel,CC.PixelIdxList);
                    [~,idx]     = sort(numPixels,'descend');

                    if isempty(idx)
                        x   = NaN;
                        y   = x;
                        return
                    end

                    if biggest_n>length(idx)
                        biggest_n=length(idx);
                    end

                    list    = idx(1:biggest_n);

                    for i=1:CC.NumObjects
                        if ~sum(i==list)>0
                            in(CC.PixelIdxList{i}) = 0;
                        end
                    end
                   
                    tp = false;
                    switch filters.direction
                        case 'left'
                            dir = 'last';
                        case 'up'
                            tp = true;
                            dir = 'last';
                        case 'right'
                            dir = 'first';
                        case 'down'
                            tp = true;
                            dir = 'first';
                    end
                    if tp
                        in = in';
                    end
                    scan = 1:size(in,1);
                    out = zeros(size(scan));
                    for i=scan
                        [~,idx] = find(in(i,:),1,dir);
                        if isempty(idx)
                            idx = nan;
                        end
                        out(i)=idx;
                    end
                    out(out<0)=nan;
                    out(out>size(in,2))=size(in,2);
                end
            end
            function out = f_smooth(in)
                out = in;
                try
                    % Apply smoothing
                    if filters.smooth
                        out = smooth(in,filters.smooth);
                    end

                    out = out-size(edi.dilationMatrix,2)/2-shift;

                    out(out<=0) = nan;
                    if isempty(out)
                        out = zeros(size(in,1),1);
                    end
                end
            end
            function out = f_ref(in)
                out = in;
                try
                    if isempty(edi.h_ref)
                        if filters.ref==edi.index
                            edi.h_ref = in;
                        end
                    end
                    map = out>edi.h_ref;
                    out(map) = edi.h_ref(map);
                end
            end
        end
    end
    %% METHODS File selection
    methods
        function files_endmembers(edi)
            r = unique(fileparts(edi.imds.ds.Files));
            f = cellfun(@(r) getEndMembers(edi.imds.ds,r),r,'UniformOutput',false);
            f = cat(1,f{:});

            edi.imds.ds.Files      = f;
            edi.imds.maxIndex        = numel(f);
            edi.index           = 1;
            function f = getEndMembers(imds,r)
                f = imds.Files(contains(imds.Files,r));
                d = cellfun(@dir,f);
                d = datetime({d.date});

                [~,imin]    = min(d);
                [~,imax]    = max(d);
                f           = f([imin imax]);
            end
        end
    end
    %% METHODS Misc. functions
    methods
        function [z,t]=getDimensions(edi)
            d=cellfun(@dir,edi.imds.ds.Files);
            t=datetime({d.date});
            t=t-t(1);
            z=(0:size(edi.loadimage(1),1)-1)*edi.ds;
        end
        function [m,h,t,f] = getMeltingData(edi,varargin)
            p = inputParser;
            addParameter(p,'normalize',false)
            addParameter(p,'type','integrate',@(x) any(validatestring({'integrate','bulk'})))
            parse(p,varargin{:})
            type      = p.Results.type;
            normalized= p.Results.normalize;
            if nargin<2
                type = 'integrate';
            end
            out = edi.a_output;
            t   = days(out.t);
            h   = out.h;
           
            switch type
                case 'integrate'
                    % depth-integrated fraction
%                     h = h(:,[1:5 end-4:end]);
%                     t = t([1:5 end-4:end]);

                    % omit if 10% of column data is nan
                    omit       = sum(isnan(h))/size(h,1)>.1;
                    h(:,omit)  = [];
                    t(omit)    = [];

                    h       = sum(h,'omitnan')*edi.ds;
                    if normalized
                        h = h/h(1);
                    end
                    
%                     h       = 1-h/max(h);
%                     omit    = h>1/exp(1);
%                     h(omit) = [];
%                     t(omit) = [];
%                     h       = h-h(1);
%                     t(h<0)  = [];
%                     h(h<0)  = [];

                    f  = fit(t',h','poly1','Robust','Bisquare');
                    m  = coeffvalues(f);
                    m  = -m(1);
                case 'bulk'
                    % bulk difference
                    dh = mean(mean(h(:,1:5),'omitnan')-mean(h(:,end-4:end),'omitnan'));
                    dt = mean(t(1:5))-mean(t(end-4:end));
                    m  = -dh/dt;
            end

        end
    end
    methods (Access = private)
        function index=valid_index(edi,index)
            index(index<1)=[];
            if isempty(index)
                index = 1;
            end
            index(index>edi.imds.maxIndex)=[];
            if isempty(index)
                index = edi.imds.maxIndex;
            end
            
        end
    end
    %% METHODS Plotting
    methods
        function plot_timeseries(edi,type,varargin)
            p = inputParser;
            addParameter(p,'ax',gca)
            addOptional(p,'normalized',false)
            addParameter(p,'color','k',@(x) any(validatestring(x,{'k','b','r','g','c','m','w','y'})))
            parse(p,varargin{:})
            ax      = p.Results.ax;
            color   = p.Results.color;
            normalized = p.Results.normalized;

            out = edi.a_output;
            if isempty(out)
                error('No timeseries output found')
            end
            if nargin<2
                type = '';
            end
            % plot
            t = days(out.t);
            h = out.h;
            switch type
                case ''
                    plot(ax,h,out.z)
                    addColorbar('ax',ax, ...
                        'lineplot', ...
                        'reverse', ...
                        'cmap','ice', ...
                        'linemap',t, ...
                        'title','$t$ (h)','latex')
                    addlabels(ax,'x','$h$ (m)','y','$z$ (m)')
                case 'pcolor'
                    pcolor(ax,t,out.z,h)
                    shading(ax,'flat')
                    addlabels(ax,'x','$t$ (d)','y','$z$ (m)')
                case 'sum'
                    [~,h,t,f] = getMeltingData(edi, ...
                        'normalize',normalized);
                    %mb = edi.getMeltingData('bulk');
                    plot(ax,t,h,color, ...
                        t,polyval(coeffvalues(f),t),'--b')
                    %,...
                     %   t,polyval(coeffvalues(f).*[0 1]-[mb 0],t),'--r')
                    addlabels(ax,'x','$t$ (d)','y','$\phi H$')
                case 'm'
                    [~,h,t] = getMeltingData(edi, ...
                        'normalize',normalized);
                    h2 = smooth(h,50);
%                     h2(gradient(h2)<-1e-4)=nan;
                    mb = edi.getMeltingData('bulk');
                    plot(ax,t,-gradient(h2,t),color)
                    line(ax,xlim(ax),mb*[1 1],'LineStyle','--','Color','r')
                    addlabels(ax,'x','$t$ (d)','y','$\partial_t (\phi H)$')
                    ylim(ax,[0 inf])
            end
        end
    end
    %% Class linking
    methods
        function link_labexperiment(edi)
            %%LINK_LABEXPERIMENT attempts to find .mat files in the same
            %%'savepath' labelled as 'experimentdata.mat'
            mf = what(edi.savepath);
            if isempty(mf.mat)
                warning('No .mat files in %s',edi.savepath)
                return
            end
            valid = contains(mf.mat,'experimentdata');
            if isempty(valid)
                warning('No "experimentdata.mat" files in %s',edi.savepath)
                return
            end
            edi.labExperiment = load(string(fullfile(mf.path,mf.mat(valid)))).le;
        end
    end
end