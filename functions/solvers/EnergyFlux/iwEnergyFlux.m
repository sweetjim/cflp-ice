classdef iwEnergyFlux<handle
    %ENERGYFLUX Summary of this class goes here
    %   Detailed explanation goes here

    properties  % measurable
        N       double = 0              % buoyancy frequency    [rad/s]
        dt      duration = seconds(1)   % time step             [s/frame]
        x       (:,1) double            % horizontal vector     [m]
        z       (:,1) double            % vertical vector       [m]
        X       double                  % horizontal mesh       [m]
        Z       double                  % vertical mesh         [m]
    end
    properties (Hidden) % computational
        derivativeMethod char = 'hs'    % time derivative method
        useBuffer   logical = false     % window buffering state
        useDiffuse  logical = false     % diffusive window buffering state
        modeRange   double = [1 inf]    % spectral mode range (default is auto)
    end
    properties (Hidden) % program
        root    {mustBeFolder(root)} = pwd
        ic      logical = false         % internal command boolean
    end

    properties (Access = private)
        rho0    double = 999.8          % reference density     [kg/m3]
        g       double = 9.81           % gravitational acc.    [m/s2]
        dx      double = 1              % horizontal resolution [m/px]
        dz      double = 1              % vertical resolution   [m/px]
    end

    %% CONSTRUCTOR
    methods
        function ef = EnergyFlux(root)
            ef.root = root;
        end
    end
    %% CALCULATIONS
    methods
        function [Jx,Jz] = calculateEnergyFlux(ef,rho)

            % velocity field calculation
            [u,w] = calculateVelocityField(ef,rho);

            % pressure field calculation
            p = calculatePressureField(ef,rho);

            % energy fluxes
            Jx = u.*p;
            Jz = w.*p;

        end
        function p = calculatePressureField(ef,rho)
            if ef.useBuffer == 1
                p = nestedCalculation(ef,rho);
                p = p(ind(1):ind(2),ind(3):ind(4),:);
            else
                p = nestedCalculation(ef,rho);
            end

            function p = nestedCalculation(ef,rho)
                DZ  = ef.dz;
                X0  = ef.X;
                Z0  = ef.Z;
                N0  = ef.N;

                modeStart   = ef.modeRange(1);
                modeEnd     = ef.modeRange(2);

                [~,dRho_dZ] = gradient(rho,DZ);

                a = N0^2/2/ef.g;

                % Set domain parameters

                % Set minimium points to zero
                xGrid0 = X0 - min(X0(:));
                yGrid0 = Z0 - min(Z0(:));

                % Domain length
                l = max(xGrid0(:));

                % Domain height
                h = max(yGrid0(:));

                % Create vertical coordinates z and z'

                z0 = yGrid0(:,1);
                z_prime = z0;

                [z0,z_prime] = meshgrid(z0,z_prime);

                % Calculate the source term

                f_source = (N0^2.*rho + ef.g*dRho_dZ).*exp(a*yGrid0);
                f = fft(f_source,[],2)/(size(f_source,2))*2;

                if isempty(modeStart) == 1
                    modeStart = 1;
                end

                fProfile = mean(abs(f));

                % calculate last meaningful spectral mode contribution (<1%)
                if isempty(modeEnd)||isinf(modeEnd)
                    [mx, mi] = max(fProfile);
                    possibleEndIndex = find(fProfile < mx*.01);
                    f_end = find(possibleEndIndex > mi, 1, 'first');
                    modeEnd = possibleEndIndex(f_end);
                    if isempty(modeEnd)
                        modeEnd = floor(size(f,2)/2);
                    end
                end

                Fki_store_fast = -imag(f(:,2:modeEnd+1));
                Fkr_store_fast = real(f(:,2:modeEnd+1));


                % Create base pressure
                pGrid_gf = 0*xGrid0;

                % Set mode number
                %modeCount = modeEnd-modeStart+1;

                for n = modeStart:modeEnd
                    % Create mode specific parameters: kx, kappa, gamma

                    kx      = 2*pi*n/l;
                    kappa   = sqrt(kx^2 + N0^4/4/ef.g^2);
                    gamma   = -4*kx^2*kappa*sinh(kappa*h);

                    % Calculate G_k pos = z>z', neg = z<z'
                    Gk_pos_term1 = (kappa + a)^2 * exp(kappa * (z0+z_prime-h));
                    Gk_neg_term1 = (kappa + a)^2 * exp(kappa * (z_prime+z0-h));

                    Gk_pos_term2 = 2*kx^2 * cosh( kappa * (z0-z_prime-h));
                    Gk_neg_term2 = 2*kx^2 * cosh( kappa * (z_prime-z0-h));

                    Gk_pos_term3 = (kappa - a)^2 * exp(kappa * (-(z0+z_prime)+h));
                    Gk_neg_term3 = (kappa - a)^2 * exp(kappa * (-(z_prime+z0)+h));

                    Gk_pos = (Gk_pos_term1 + Gk_pos_term2 + Gk_pos_term3)/gamma;
                    Gk_neg = (Gk_neg_term1 + Gk_neg_term2 + Gk_neg_term3)/gamma;

                    Gk = Gk_pos;

                    f = find(z0<z_prime);

                    Gk(f) = Gk_neg(f);

                    % Calculate the fourier components for the given wave number

                    Fkr = Fkr_store_fast(:,n);
                    Fki = Fki_store_fast(:,n);


                    % Calculate the convolution value (G_k * F)

                    Gkintegrandr = Gk.*(ones(length(Fkr),1)*Fkr');
                    Gkintegrandi = Gk.*(ones(length(Fki),1)*Fki');


                    term1_fast = trapz(yGrid0(:,1),Gkintegrandr')'*cos(kx*xGrid0(1,:));
                    term2_fast = trapz(yGrid0(:,1),Gkintegrandi')'*sin(kx*xGrid0(1,:));

                    deltaP_gf = -exp(-N0^2/2/ef.g*yGrid0).*(term1_fast+term2_fast);

                    pGrid_gf = pGrid_gf + deltaP_gf;

                end

                p = pGrid_gf;
            end
        end
        function [u,w] = calculateVelocityField(ef,rho)
            if size(rho,3)==1
                error('Density perturbation field must be non-unitary')
            end
            
            DZ      = ef.dz;
            DX      = ef.dx;
            DT      = seconds(ef.dt);
            N0      = ef.N;
            G       = ef.g;
            RHO0    = ef.rho0;


            % Calculate density perturbation time derivative
            nSlices = size(rho,3);
            dRho_dt = zeros(size(rho,1),size(rho,2),nSlices);

            for j = 1:nSlices
                if j == 1
                    dRho_dt(:,:,j) = (-rho(:,:,1) + rho(:,:,2))/DT;
                elseif j == nSlices
                    dRho_dt(:,:,j) = (-rho(:,:,end-1) + rho(:,:,end))/DT;
                elseif j == 2
                    dRho_dt(:,:,j) = (-rho(:,:,1) + rho(:,:,3))/DT/2;
                elseif j == nSlices-1
                    dRho_dt(:,:,j) = (-rho(:,:,end-2) + rho(:,:,end))/DT/2;
                else
                    dRho_dt(:,:,j) = ( rho(:,:,j-2)/12 - 2*rho(:,:,j-1)/3 + 2*rho(:,:,j+1)/3 - rho(:,:,j+2)/12 )/DT;
                end
            end

            % Calculation of the vertical velocity

            w = G./N0.^2/RHO0.*dRho_dt;

            % Vertical velocity gradient

            dw_dz_exp = permute(fourth_grad(permute(w,[2 1 3]),DZ),[2 1 3]);

            % Horizontal integration

            for j = 1:nSlices
                u1(:,:) = -cumsum(dw_dz_exp(:,:,j),2)*DX;
                u2(:,:) =  cumsum(dw_dz_exp(:,end:-1:1,j),2)*DX; 
                u2      =  u2(:,end:-1:1);


                % Identify where the boundary is good to start integration from

                rStartMean  = zeros(1,size(rho,1));
                rEndMean    = zeros(1,size(rho,1));
                for i = 1:size(rho,1)

                    ind = i-2:i+2;
                    f = find(ind>0 & ind<=size(rho,1));
                    ind = ind(f);

                    rStartMean(i) = mean(abs(rho(ind,1,j)));
                    rEndMean(i) = mean(abs(rho(ind,end,j)));

                end

                indStart = find(rStartMean<rEndMean);
                indEnd = find(rEndMean<=rStartMean);

                % Set the horizontal velocity

                u(indStart,:,j) = u1(indStart,:);
                u(indEnd,:,j) = u2(indEnd,:);

                % Smooth the velocity fields

                for k = 1:3
                    for i = 1:size(u,2)
                        u(:,i,j) = smooth(u(:,i,j),5);
                    end
                end

            end
            function df_dx = fourth_grad(f,dx)
                df_dx = 0*f;
                for jj = 1:size(f,3)
                    df_dx(:,1,jj)        = -3/2*f(:,1,jj) + 2*f(:,2,jj) - f(:,3,jj)/2;
                    df_dx(:,2,jj)        = -f(:,1,jj)/2 + f(:,3,jj)/2;
                    df_dx(:,3:end-2,jj)  = 1/12*f(:,1:end-4,jj) - 2/3*f(:,2:end-3,jj) + 2/3*f(:,4:end-1,jj) - 1/12*f(:,5:end,jj);
                    df_dx(:,end-1,jj)    = -f(:,end-2,jj)/2 + f(:,end,jj)/2;
                    df_dx(:,end,jj)      = 3/2*f(:,end,jj) - 2*f(:,end-1,jj) + f(:,end-2,jj)/2;
                end
                df_dx = df_dx/dx;
            end
        end
        function [xGrid,yGrid,rho,ind] = bufferData(ef,rho, gpX, gpY)
            DX      = ef.dx;
            xGrid   = ef.X;
            yGrid   = ef.Z;

            % set buffer-window size (fraction of image size)
            if nargin<3
                gpX = .1;
            end
            if nargin<4
                gpY = .1;
            end


            bufferSizeX = (max(xGrid(:))-min(xGrid(:)))*gpX;
            nPx = ceil(bufferSizeX/DX);
            %bufferSizeX = nPx*DX;
            bufferX = DX*[1:nPx]; %#ok<*NBRAK1>

            bufferSizeY = (max(yGrid(:))-min(yGrid(:)))*gpY;
            nPy = ceil(bufferSizeY/DX);
            %bufferSizeY = nPy*DX;
            bufferY = DX*[1:nPy];

            xBuffer = [min(xGrid(:))-bufferX(end:-1:1) xGrid(1,:) max(xGrid(:))+bufferX];
            yBuffer = [min(yGrid(:))-bufferY(end:-1:1) yGrid(:,1)' max(yGrid(:))+bufferY];

            [xGridWindow,zGridWindow] = meshgrid(xBuffer,yBuffer);

            [s1, s2] = size(xGrid);
            f_grid_out = boarder_points(1,1,s1,s2)';

            [s1, s2] = size(xGridWindow);
            f_out_toInterpolate = boarder_points(nPy, nPx, s1, s2);

            xI = [xGrid(f_grid_out); xGridWindow(1,2:end-1)'; xGridWindow(end,2:end-1)'; xGridWindow(:,1); xGridWindow(:,end)];
            yI = [yGrid(f_grid_out); zGridWindow(1,2:end-1)'; zGridWindow(end,2:end-1)'; zGridWindow(:,1); zGridWindow(:,end)];

            rpg_Window = zeros(size(xGridWindow,1),size(xGridWindow,2),5);

            % h_bar = waitbar(0,'Buffering Data: 0%');

            if ef.useDiffuse == 1
                warning('off','MATLAB:scatteredInterpolant:DupPtsAvValuesWarnId');
                for k = 1:size(rho,3)
                    rpg = rho(:,:,k);
                    rhoPertGridWindow = zeros(s1,s2);


                    rI = [rpg(f_grid_out); 0*xGridWindow(1,2:end-1)'; 0*xGridWindow(end,2:end-1)'; 0*xGridWindow(:,1); 0*xGridWindow(:,end)];

                    rhoPertGridWindow(f_out_toInterpolate) = griddata(xI,yI,rI,xGridWindow(f_out_toInterpolate),zGridWindow(f_out_toInterpolate));

                    rhoPertGridWindow(1+nPy:end-nPy,1+nPx:end-nPx) = rho(:,:,k);

                    for i = 1:10


                        for j = 1:size(rhoPertGridWindow,1)
                            rhoPertGridWindow(j,:) = smooth(rhoPertGridWindow(j,:),max([1 ceil(nPx/3)]));
                        end
                        for j = 1:size(rhoPertGridWindow,2)
                            rhoPertGridWindow(:,j) = smooth(rhoPertGridWindow(:,j),max([1 ceil(nPy/3)]));
                        end

                        rhoPertGridWindow(1+nPy:end-nPy,1+nPx:end-nPx) = rho(:,:,k);

                    end

                    rpg_Window(:,:,k) = rhoPertGridWindow;
                    waitbar(k/size(rho,3),h_bar,sprintf('Buffering Data: %d%%',k/size(rho,3)*100))

                end

                warning('on','MATLAB:scatteredInterpolant:DupPtsAvValuesWarnId');

            else
                rpg_Window = zeros(size(xGridWindow,1),size(xGridWindow,2),size(rho,3));

                for k = 1:size(rho,3)
                    rpg_Window(1+nPy:end-nPy,1+nPx:end-nPx,k) = rho(:,:,k);
                    % waitbar(k/size(rho_pert,3),h_bar,sprintf('Buffering Data: %d%%',k/size(rho_pert,3)*100))
                end

            end

            % close(h_bar)

            xGrid = xGridWindow;
            yGrid = zGridWindow;
            rho = rpg_Window;

            ind = [nPy+1 size(rho,1)-nPy nPx+1 size(rho,2)-nPx];

            function f_out = boarder_points(nPy, nPx, s1, s2)

                [ind1, ind2] = meshgrid(1+nPy:s1-nPy,1+nPx:s2-nPx);
                f = sub2ind([s1 s2],ind1,ind2);
                f_out = setdiff(1:s1*s2,f(:));

            end
        end
        function [fx,fz,ft] = getDerivatives(~,im1,im2,method)
            if nargin<3
                method = 'hs';
            end

            if size(im2,1)==0
                im2=zeros(size(im1));
            end
            useconv     = @(im1,im2,kern) conv2(im1,kern,'same') + conv2(im2, kern,'same');
            useimfilt   = @(im1,im2,kern) imfilter((im1+im2)/2, kern, 'symmetric',  'same');
            kernt       = 1/4*ones(2);
            switch method
                case 'hs'
                    % Horn-Schunck original method
                    kern = 1/4*...
                        [-1 1
                        -1 1];
                    fx = useconv(im1,im2,kern);
                    fz = useconv(im1,im2,kern');
                    ft = conv2(im1,kernt,'same') + conv2(im2,-kernt,'same');
                case 'barron'
                    % derivatives as in Barron (not advisable)
                    kern = 1/12*...
                        [-1 8 0 -1 1];
                    fx= conv2(im1,kern,'same');
                    fz= conv2(im1,kern','same');
                    ft = conv2(im1,kernt,'same') + conv2(im2,-kernt,'same');
                case 'diff'
                    % An alternative way to compute the spatiotemporal derivatives is to use simple finite difference masks.
                    fx = conv2(im1,[1 -1],'same');
                    fz = conv2(im1,[1; -1],'same');
                    ft= im2-im1;
                case 'ls'
                    % derivates as in Liu-Shen
                    kern = 1/2*...
                        [0  0  0
                        0 -1 -1
                        0  1  1];
                    kernt = 1/4*...
                        [0  0  0
                        0  1  1
                        0  1  1];
                    %
                    fx = useimfilt(im1,im2,kern');
                    fz = useimfilt(im1,im2,kern);
                    ft = imfilter(im2-im1, kernt, 'symmetric',  'same');
            end
            fx = fx/ef.dx;
            fz = fz/ef.dz;
            ft = ft/ef.dt;
        end
    end
    %% SET and ROUTINES
    methods
        function set.x(ef,val)
            ef.x = val;
            getDimensions(ef)
        end
        function set.z(ef,val)
            ef.z = val;
            getDimensions(ef)
        end
        function set.X(ef,val)
            ef.X = val;
            setDimensions(ef)
        end
        function set.Z(ef,val)
            ef.Z = val;
            setDimensions(ef)
        end
    end
    methods (Hidden, Access = private)
        function getDimensions(ef)
            if ef.ic; return;end
            if any([isempty(ef.x) isempty(ef.z)])
                return
            end
            ef.ic = true;
            [ef.X,ef.Z] = meshgrid(ef.x,ef.z);
            ef.ic = false;
        end
        function setDimensions(ef)
            if ef.ic; return;end
            if any([isempty(ef.X) isempty(ef.Z)])
                return
            end
            ef.ic = true;
            ef.x = min(ef.X);
            ef.z = min(ef.Z);
            ef.ic = false;
        end
    end
end

