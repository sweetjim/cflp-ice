classdef imMerge<handle
    %IMMERGE merges two ImageDataStore objects
    
    properties (Access = public)
        ds1
        ds2
        dates1 datetime
        dates2 datetime
        Folders cell
        Files cell
        maxIndex (1,1) double
        mergeSettings struct
    end
    properties (Access = private,Hidden)
        files1
        files2
    end
    
    methods
        function imm = imMerge(datastore1,datastore2)
           condition = all([
               isa(datastore1,'matlab.io.datastore.ImageDatastore')
               isa(datastore2,'matlab.io.datastore.ImageDatastore')]);
           if ~condition
               error('Both arguments must be ImageDatastore types')
           end
           imm.ds1 = datastore1;
           imm.ds2 = datastore2;
           imm.Folders = [datastore1.Folders datastore2.Folders];
           imm.getFiles;
           imm.getDates;
        end
        function [im,varargout] = readImageAtIndex(imm,index,applyTransform)
            if nargin<3
                applyTransform = true;
            end
            fc1 = numel(imm.ds1.Files);
            fc2 = numel(imm.ds2.Files);

            if numel(index)>1
                % read multiple files simultaenously
                [];
            end

            try 
                im1 = readimage(imm.ds1,index);
            catch
                im1 = readimage(imm.ds1,fc1);
            end
            try 
                im2 = readimage(imm.ds2,index);
            catch
                im2 = readimage(imm.ds2,fc2);
            end
            if applyTransform
                im1 = transformImage(imm,im1,'left');
                im2 = transformImage(imm,im2,'right');
            end
            if nargout>1
                im = im1;
                varargout{1}=im2;
                return
            end
            im = imm.stitchImages(im1,im2);
            
        end
        function [im,varargout] = readImageAtTime(imm,time,applyTransform)
            if nargin<3
                applyTransform = true;
            end
            condition = [isa(time,'duration') isa(time,'double') isa(time,'datetime')];
            if ~any(condition)
                error('Time argument value must be a valid duration, double, or datetime')
            end
            if any(condition([1 2]))
                if time<0
                    error('Time argument must be positive')
                end
            end
            if condition(2)
                time = seconds(time);
            end

            if numel(time)>1
                % read multiple files simultaenously
                [];
            end

            t0 = min([min(imm.dates1) min(imm.dates2)]);
            idx1 = find(imm.dates1>=t0+time,1,'first');
            idx2 = find(imm.dates2>=t0+time,1,'first');
            idx1 = find(ismember(imm.ds1.Files,imm.files1(idx1)));
            idx2 = find(ismember(imm.ds2.Files,imm.files2(idx2)));
            im1 = readimage(imm.ds1,idx1);
            im2 = readimage(imm.ds2,idx2);

            if applyTransform
                im1 = transformImage(imm,im1,'left');
                im2 = transformImage(imm,im2,'right');
            end
            if nargout>1
                im = im1;
                varargout{1}=im2;
                return
            end
            im = imm.stitchImages(im1,im2);
        end
        function getDates(imm)
            [imm.dates1,imm.files1,map1] = fileDatetime(imm.Folders{1},'children');
            [imm.dates2,imm.files2,map2] = fileDatetime(imm.Folders{2},'children');

            omit1 = contains(imm.files1,imm.ds1.SupportedOutputFormats);
            omit2 = contains(imm.files2,imm.ds1.SupportedOutputFormats);
            
            imm.dates1  = imm.dates1(omit1);
            imm.dates2  = imm.dates2(omit2);
            imm.files1  = imm.files1(omit1);
            imm.files2  = imm.files2(omit2);
            map1        = map1(omit1);
            map2        = map2(omit2);

            valid1 = ismember(imm.Files,imm.files1);
            valid2 = ismember(imm.Files,imm.files2);
%             imm.Files(valid1) = imm.files1(map1);
%             imm.Files(valid2) = imm.files2(map2);
        end
        function getFiles(imm)
            imm.Files       = cat(1,imm.ds1.Files,imm.ds2.Files);
            imm.maxIndex    = min([numel(imm.ds1.Files) numel(imm.ds2.Files)]);

            % look for timeindex tag 
            % filename format: cam<id>_C<count>T<hh_mm_ss>
            [~,f1] = fileparts(imm.ds1.Files);
            [~,f2] = fileparts(imm.ds2.Files);
            if all(~contains(f1,'cam'))
                return
            end
            f1C = cellfun(@str2double,extractBetween(f1,'C','T'));
            f2C = cellfun(@str2double,extractBetween(f2,'C','T'));
            [~,map1] = sort(f1C);
            [~,map2] = sort(f2C);
            imm.ds1.Files = imm.ds1.Files(map1);
            imm.ds2.Files = imm.ds2.Files(map2);
        end
        
    end
    methods (Access=public,Hidden)
        function im = stitchImages(imm,im1,im2,direction)
            if numel(im1)~=numel(im2)
                im = [];
                return
            end
            if nargin<4
                direction = 1;
            end
            ms = imm.mergeSettings;
            if ~isempty(fieldnames(ms))
                direction = ms.catDirection;
            end
            try
                im = cat(direction,im1,im2);
            catch
                im = cat(direction,im1,im2);
            end
        end
        function image = transformImage(imm,image,side)
            ms = imm.mergeSettings;
            if isempty(fieldnames(ms))
                return
            end
            switch side
                case 'left'
                    ud = ms.leftSide;
                case 'right'
                    ud = ms.rightSide;
            end
            if isempty(ud)
                return
            end
            if isempty(image)
                return
            end

            if ud.ShiftX
                image = circshift(image,ud.ShiftX,1);
            end
            if ud.ShiftY
                image = circshift(image,ud.ShiftY,2);
            end
            
            image = rot90(image,ud.Rot);
            if getlogical(ud.FlipLR)
                image = fliplr(image);
            end
            if getlogical(ud.FlipUD)
                image = flipud(image);
            end
        end
    end
end

