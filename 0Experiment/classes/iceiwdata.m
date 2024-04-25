classdef iceiwdata<handle
    %ICEIWDATA Summary of this class goes here
    %   Detailed explanation goes here
    
    properties
        id icedata
        iw iwdata
        N double = 0.3
    end
    
    methods
        function iciw = iceiwdata(iceData,iwData)
            iciw.id                 = iceData;
            if nargin<2
                iwData = iwdata;
            end
            iciw.iw                 = iwData;
            iciw.id.Dimensions.ds   = 1/iwData.ds;
            if iciw.iw.N>0
                iciw.N = iciw.iw.N;
            end
        end
    end

    %% Melting rates
    methods
        function [refs,M,A,OMEGA] = melting_decompose(iciw)
            %%
            data        = iciw.id;
            [H,~,T]     = data.parseAnalysis;
            H(H<0)      = nan;
            H(H>.2)     = nan;

            afcn        = @(fcn,array) arrayfun(fcn,array,'UniformOutput',false);
            cfcn        = @(fcn,array) cellfun(fcn,array,'UniformOutput',false);
            cats        = @(refs,newrefs,newreflbl) catstruct(refs,cell2struct(newrefs,newreflbl));

            % Get global time
            gtime       = data.times;

            % Find wave frequencies+amplitudes from filenames
            amp         = strrep(strrep(extractBetween({data.Folder.Files.name},'amp','C'),'-',''),'_','');
            amps        = str2double(amp);
            om          = strrep(extractBetween({data.Folder.Files.name},'omegaN','amp'),'-','');
            om          = cfcn(@(om) str2double(join(split(om,'_'),'.')),om);

            % Decompose edge-detection data by amplitudes
            tmp         = find(abs(gradient(amps))>0);
            idx         = tmp([diff(tmp)==1 false]);
            amps0       = amps(idx);

            % Find absent 0-amplitude pictures in global time
            [];

            %%
            om          = om(idx);
            om(amps0==0)={nan};

            % Get index range
            idx         = [circshift(idx+1,1); idx]';
            idx(1,1)    = 1;
            ranges      = afcn(@(i) idx(i,1):idx(i,2),1:height(idx));
            amps0       = afcn(@(x) x,amps0);

            refs        = cell2struct([amps0;ranges;om],{'amp','Ti','omegaN'});

            % Get depth-integrated melt fraction
            phi         = 1-sum(H,'omitnan')/sum(H(:,1),'omitnan');
            phiS        = afcn(@(x) phi(x.Ti),refs)';

            % Relative to first index
            phiS0       = cfcn(@(x) x-x(1),phiS);
            Tr          = afcn(@(x) T(x.Ti),refs)';
            T0          = afcn(@(x) T(x.Ti)-T(x.Ti(1)),refs)';
            refs        = cats(refs,[phiS;phiS0;Tr;T0],{'phi','phi0','t','t0'});
            
            % Relative to maximum time of unforced segments
            tmax        = min(arrayfun(@(x) max(x.t0),refs([refs.amp]==0)));
            if isempty(tmax)
                tmax = min(cell2mat(afcn(@(x) max(x.t0),refs)));
            end
            tmaxi       = afcn(@(x) find(round(x.t0/10)*10>=round(tmax/10)*10,1),refs)';
            refs        = cats(refs,tmaxi,'tmax');

            % Relative melting rate between 0<t<tmax
            m           = afcn(@(x) x.phi0(x.tmax)/x.t0(x.tmax),refs)';
            refs        = cats(refs,m,'dt0dphi0');

            % Outputs
            M           = [refs.dt0dphi0];
            A           = [refs.amp];
            OMEGA       = [refs.omegaN]*iciw.iw.N;

        end
    end
end

