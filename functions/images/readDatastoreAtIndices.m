function images = readDatastoreAtIndices(datastore,indices,outputtype)
%READDATASTOREATINDICES 
if nargin<3
    outputtype = 'matrix';
end
hasParallel = false;%isvalid(gcp('nocreate'));
switch class(datastore)
    case 'imMerge'
        % combine datastores
        ds = combine(datastore.ds1,datastore.ds2);

        % get datastore with correct indices
        ds = subset(ds,indices);

        % read images
        images = cellfun(@(ds) readall(ds,UseParallel=hasParallel),ds.UnderlyingDatastores, ...
            'ErrorHandler',@(~,~) uint8(nan), ...
            'UniformOutput',false);

        % attempt transform operations
        images = cellfun(@(idx,side) ...
                    cellfun(@(image) datastore.transformImage(image,side),images{idx},'UniformOutput',false), ...
                 {1,2},{'left','right'},'UniformOutput',false);

        % array stacking
         images = [images{:}];
        % stitch images together
        images = arrayfun(@(idx) datastore.stitchImages(images{idx,1},images{idx,2}),1:height(images),'UniformOutput',false)';

    case 'matlab.io.datastore.ImageDatastore'
        ds = subset(datastore,indices);
        images = readall(ds,UseParallel=hasParallel);
    otherwise
        error('Datastore argument incorrect class type')
end

% convert to double
images = cellfun(@im2double,images,'UniformOutput',false);

switch outputtype
    case 'cell'
        return
    case 'matrix'
        im = zeros([size(images{1}),height(images)]);
        for i=1:size(im,3)
            im(:,:,i) = images{i};
        end
        images = im;
end

    function readAllProgress(ds)
    end
end

