classdef MushyLayerOutput < ChomboOutput
    properties
        subcycled
    end
    
    methods
        function m = MushyLayerOutput(dim, frame, output_dir, plot_prefix, subcycled, interpolationMethod)
            if nargin < 6
                % other options bilinear, bicubic
                interpolationMethod = 'nearest';
            end
            
            m@ChomboOutput(dim, frame, output_dir, plot_prefix, interpolationMethod);
            
            m.subcycled = subcycled;
            
            
            
        end
        
        
        % Rewriting as AMR version
        function [Nu_av, Nu_base] = nusseltNumber(obj, dir)
            if nargin < 2
                dir = 1;
            end
            
            levels = obj.levelArray;
            
            Nu_av = 0; Nu_base = 0;
            num_levs = length(levels);
            
            %Force code to only use coarse level by fixing num_levs = 1
            %num_levs = 1;
            for lev_i = num_levs:-1:1
                
                level = levels(lev_i);
                
                dx = level.dx;  dz = level.dx;
                
                if obj.subcycled
                    theta =level.dataForComp(obj.components.temperature);
                    Vz = level.dataForComp(obj.components.yDarcyVel);
                    Vx = level.dataForComp(obj.components.xDarcyVel);
                else
                    
                    theta =level.dataForComp(obj.components.theta);
                    Vz = level.dataForComp(obj.components.yFluidVel);
                    Vx = level.dataForComp(obj.components.xFluidVel);
                end
                
                [X, Z] = level.grid_no_ghost();
                X_arr = X(1,:);
                Z_arr = Z(:, 1);
                
                %                [dtheta_dx, dtheta_dz] = gradient2order(theta, dx,dx);
                [dtheta_dx, dtheta_dz] = gradient(theta, dx,dx);
                
                L_x = X_arr(end) - X_arr(1);
                L_z = Z_arr(end) - Z_arr(1);
                
                % dir = 1 is vertical heat transport
                if dir == 1
                    integrand = (-dtheta_dz.' + (Vz.*theta).')*dx;
                    
                    
                    %integrand = 0*integrand + 1;
                    
                    % Get x and y the right way around
                    
                    % Remove regions covered by finer levels
                    for finer_lev_i = lev_i + 1:num_levs
                        finer_lev =  levels(finer_lev_i);
                        level_coarsen = finer_lev.dx/level.dx;
                        for finer_box_i = 1:length(finer_lev.boxArray)
                            finer_box = finer_lev.boxArray(finer_box_i);
                            [lo_i, lo_j, hi_i, hi_j] = finer_box.get_extent_no_ghost(level_coarsen);
                            
                            integrand(lo_i:hi_i, lo_j:hi_j) = NaN*ones(hi_i+1-lo_i, hi_j+1-lo_j);
                        end
                    end
                    
                    
                    y_num = size(integrand, 2);
                    x_num = size(integrand, 1);
                    
                    Nu = zeros(y_num, 1);
                    
                    % Iterate over y-direction
                    for y_i = 1:y_num
                        
                        % Integrate in the x-direction
                        
                        % Need to only integrate over valid (not-nan)
                        % regions. To do this, construct vector of valid
                        % values
                        valid_integrand = [];
                        for x_i =1:x_num
                            val = integrand(x_i, y_i);
                            if isnan(val)
                                % If we've reached a nan value, integrate
                                % over the values we already have
                                Nu_valid_region = trapz(valid_integrand)*dx;
                                Nu(y_i) = Nu(y_i) + Nu_valid_region;
                                
                                % Reset our vector of valid values
                                valid_integrand = [];
                            else
                                valid_integrand(end+1) = val;
                            end
                        end
                        
                        % Integrate over final valid region (if applicable)
                        Nu_valid_region = trapz(valid_integrand)*dx;
                        Nu(y_i) = Nu(y_i) + Nu_valid_region;
                        
                        
                    end
                    
                    %Calculate Nu(y) by integrating in the x direction for each y point
                    %Nu = trapz(integrand.')*dx;
                    
                    %Now average over each y level - should do this in a
                    %clever weighted way
                    Nu_av = Nu_av + mean(Nu);
                    
                    %integrand = - dtheta_dz.'/L_x;
                    integrand = (- dtheta_dz.' + (Vz.*theta).' )/(L_x);
                    % Calculate integrand differently:
                    % Just calculate dtheta_dz at z= 0 (as v=0 at z=0)
                    % Need to use one sided difference approx'
                    
                    integrand_base = zeros(x_num, 1);
                    
                    bcval = 1;
                    
                    for x_i = 1:x_num
                        farVal = theta(2, x_i);
                        nearVal = theta(1, x_i);
                        ghostVal  = (8/3)*bcval + (1/3)*farVal - 2*nearVal;
                        %integrand_base(x_i) = ( -2*theta(end, x_i) + 3*theta(end-1, x_i) - theta(end-2, x_i) ) / (dx);
                        integrand_base(x_i) = (nearVal - ghostVal)/dx;
                    end
                    
                    %integrand = ...;
                    Nu = trapz(integrand)*dx;
                    Nu_base = Nu(1);
                    Nu_base = -trapz(integrand_base)*dx/L_x;
                    Nu_base = -sum(integrand_base)*dx/(L_x+dx);
                    Nu_av = sum(Nu)/length(Nu);
                    
                elseif dir == 0
                    integrand = (-dtheta_dx + Vx.*theta)/L_z;
                    
                    %Calculate Nu(y) by integrating in the x direction for each y point
                    Nu = trapz(integrand.')*dz;
                    
                    %Now average over each y level
                    Nu_av = Nu_av + mean(Nu);
                    
                    % integrand = - dtheta_dx.'/L;
                    integrand = (- dtheta_dx.' + (Vx.*theta).' )/(L_z);
                    Nu = trapz(integrand)*dz;
                    Nu_base = Nu_base + Nu(1);
                end
                
                
                
            end
            
            
            
            
        end
        
        function [psi] = getStreamfunction(obj, solverIterations, level)
            if nargin < 3
                level = 1;
            end
            
            if nargin < 2
                level = 1;
                solverIterations = 10000;
            end
            
            compList = obj.components;
            if isfield(compList, 'streamfunction')
                psi = obj.dataForComp(obj.components.streamfunction).';
                return
            end
            
            if (obj.subcycled)
                V = obj.levelArray(level).dataForComp(obj.components.yAdvectionvelocity);
                U = obj.levelArray(level).dataForComp(obj.components.xAdvectionvelocity);
            else
                V = obj.levelArray(level).dataForComp(obj.components.yFluidVel);
                U = obj.levelArray(level).dataForComp(obj.components.xFluidVel);
            end
            
            if nargin < 3
                [X,Y] = obj.grid();
                %V =V.';
                % U = U.';
            else
                [X,Y]= obj.levelArray(level).grid_no_ghost();
                
            end
            
            x = X(1, :);
            y = Y(:, 1);
            
            dx = x(2) - x(1);
            
            %smooth = smoothn({U, V}, 50, 'robust');
            %U = smooth{1}; V = smooth{2};
            
            % [vorticity, cav] = curl(X, Y, U, V);
            vorticity = curl(X, Y, U, V);
            vorticity = -dx*vorticity; %ensure correct scaling
            psi = solvePoisson(vorticity, x, y, solverIterations);
        end
        
        
        
        
        function plotStreamArrows(obj, contours, color, atGrad)
            xc = contours(1,:); yc = contours(2,:);
            hold on;
            lengthGrow = 0;
            
            %atGrad = 0.0;
            tolerance = 0.06;
            
            i=2;
            
            headWidth = 10;
            headLength = 10;
            LineLength = 0.001;
            
            while i < (length(contours)-1)
                gradient =  (yc(i+1)-yc(i))/(xc(i+1)-xc(i));
                if (xc(i) < 3 && (gradient < atGrad+tolerance && gradient > atGrad-tolerance && yc(i) > 0.25) ...
                        || xc(i) > 3 && (gradient < -atGrad+tolerance && gradient > -atGrad-tolerance && yc(i) > 0.25))
                    if xc(i) < xc(i+1)
                        xinit = xc(i) - lengthGrow;
                        xfinal = xc(i+1) + lengthGrow;
                    else
                        xinit = xc(i+1) - lengthGrow;
                        xfinal = xc(i) +lengthGrow;
                    end
                    
                    % Do this to get arrow data
                    hq = quiver(xinit,yc(i),(xfinal-xinit),(yc(i+1)-yc(i)),0);
                    
                    hq.Visible = 'off';
                    
                    % Get arrow data
                    Uquiv = hq.UData;
                    Vquiv = hq.VData;
                    Xquiv = hq.XData;
                    Yquiv = hq.YData;
                    
                    xSign = -1;
                    ySign = -1;
                    %color='black';
                    if Xquiv < 3 % chimney position
                        xSign = 1;
                    end
                    
                    ah = annotation('arrow',...
                        'headStyle','cback1','HeadLength',headLength,'HeadWidth',headWidth, 'Color', color);
                    set(ah,'parent',gca);
                    
                    set(ah,'position',[Xquiv Yquiv xSign*LineLength*Uquiv ySign*LineLength*Vquiv]);
                    
                    % Move away from this region
                    i=i+20;
                end
                
                i=i+1;
            end
            hold off;
            
        end
        
        % Compute horizontally averaged porosity (z)
        function [averaged] = horizontallyAverage(obj, comp, makeAbs)
            field = obj.dataForComp(comp);
            
            if makeAbs
                field = abs(field);
            end
            
            
            
            probDomain = obj.problemDomain;
            dx = probDomain.dxCoarse;
            
            numx = (probDomain.domainExtent.hi_i - probDomain.domainExtent.lo_i) + 1;
            numy = (probDomain.domainExtent.hi_j - probDomain.domainExtent.lo_j) + 1;
            
            x = linspace(double(probDomain.domainExtent.lo_i), double(probDomain.domainExtent.hi_i), numx)*dx;
            y = linspace(double(probDomain.domainExtent.lo_j), double(probDomain.domainExtent.hi_j), numy)*dx;
            
            averaged = [];
            
            for j = 1:numy
                thisSection = squeeze(field(:, j));
                thisMean = mean(mean(thisSection));
                averaged(j) = thisMean;
                
            end
        end
        
        % Compute horizontally averaged porosity (z)
        function [maxVal] = horizontallyMax(obj, field, makeAbs)
            %field = obj.dataForComp(comp);
            
            if makeAbs
                field = abs(field);
            end
            
            
            
            probDomain = obj.problemDomain;
            dx = probDomain.dxCoarse;
            
            numx = (probDomain.domainExtent.hi_i - probDomain.domainExtent.lo_i) + 1;
            numy = (probDomain.domainExtent.hi_j - probDomain.domainExtent.lo_j) + 1;
            
            x = linspace(double(probDomain.domainExtent.lo_i), double(probDomain.domainExtent.hi_i), numx)*dx;
            y = linspace(double(probDomain.domainExtent.lo_j), double(probDomain.domainExtent.hi_j), numy)*dx;
            
            maxVal = [];
            
            for j = 1:numy
                thisSection = squeeze(field(:, j));
                thisMax = max(max(thisSection));
                maxVal(j) = thisMax;
                
            end
        end
        
        
        function [width, depth, totalNumChannels, mushHeight] = channelGeometry(obj)
            
            probDomain = obj.problemDomain;
            dx = probDomain.dxCoarse;
            
            numx = (probDomain.domainExtent.hi_i - probDomain.domainExtent.lo_i) + 1;
            numy = (probDomain.domainExtent.hi_j - probDomain.domainExtent.lo_j) + 1;
            
            x = linspace(double(probDomain.domainExtent.lo_i), double(probDomain.domainExtent.hi_i), numx)*dx;
            y = linspace(double(probDomain.domainExtent.lo_j), double(probDomain.domainExtent.hi_j), numy)*dx;
            
            [X,Y] = meshgrid(x, y);
            
            porosity = obj.dataForComp(obj.components.Porosity);
            
            % Do some smoothing on porosity
            refinement = 32;
            porosityRefined = resizem(porosity, refinement, 'bicubic');
            Xref = resizem(X, refinement, 'bilinear');
            Yref = resizem(Y, refinement, 'bilinear');
            
            [numx, numy] = size(porosityRefined);
            dxRefined = dx/refinement;
            
            liquidLimit = 1.0 - 1e-2;
            [liquidx, liquidy] =  find(porosityRefined > liquidLimit);
            [mushx, mushy] =  find(porosityRefined < liquidLimit);
            
            high_j = max(liquidy);
            
            for x_i = 1:numx
                chiForx = porosityRefined(x_i, :);
                [~, liquidVals] = find(chiForx > liquidLimit);
                mlBoundaryVals(x_i) = max(liquidVals);
            end
            
            low_j = median(mlBoundaryVals);
            
            channelDepth = (high_j - low_j)*dxRefined;
            
            [eutecticx, eutecticy]= find(porosityRefined < 0.001);
            eutectic_j = min(min(eutecticy));
            mushHeight = (eutectic_j - low_j)*dxRefined;
            
            %Channel width: measure at a few places and average
            
            
            %maxG = 8;
            % calculate width over middle 50%
            start_j = round(low_j + (high_j - low_j)*0.25);
            end_j = round(low_j + (high_j - low_j)*0.75);
            
            numChannels = 0*ones(1, end_j);
            thisChanWidths = NaN*ones(1, end_j);
            for this_j = start_j:1:end_j
                % this_j = round(low_j + (g/maxG)*(high_j - low_j));
                numChannels(this_j) = 0;
                
                % Width at this j:
                liquidCount = 0;
                maxLiquidCount = 0;
                for trace_i = 1:numx
                    if porosityRefined(trace_i, this_j) > liquidLimit
                        liquidCount = liquidCount + 1;
                    else
                        maxLiquidCount = max(maxLiquidCount, liquidCount);
                        
                        if liquidCount > 0
                            numChannels(this_j) = numChannels(this_j) + 1;
                        end
                        
                        
                        liquidCount = 0;
                        
                        
                    end
                end
                
                % Count the last channel
                if liquidCount > 0
                    numChannels(this_j) = numChannels(this_j) + 1;
                end
                
                maxLiquidCount = max(maxLiquidCount, liquidCount);
                
                thisChanWidths(end+1) = maxLiquidCount*dxRefined;
            end
            
            totalNumChannels = max(numChannels);
            
            if length(thisChanWidths) > 0
                width = nanmean(thisChanWidths);
                %width = thisChanWidths(round(end/2)); % just take the middle point
            else
                width = NaN;
            end
            depth = channelDepth;
            %optimalChannelDepth(i) =
            
            %             fprintf('width = %f, depth = %f \n', optimalChannelWidth(i), optimalChannelDepth(i));
            %             subplot(1, length(RaCs), i);
            %
            %             hold on;
            %             imagesc(x, y, porosity.');
            %             [c, h] = contour(Xref, Yref, porosityRefined.', [-0.6 liquidLimit]);
            %
            %             h.LineWidth = 2;
            %             h.LineColor = 'r';
            %
            %
            %
            %             hold off;
            %
            %             ax = gca;
            %
            %             %ax.XTick = [];
            %
            %             if i > 1
            %                 ax.YTick = [];
            %             end
            %             set(ax,'clim',[0 1]);
            %             set(ax, 'Ydir', 'normal');
            %
            %             title(num2str(RaC));
            % %
            %              pause;
            
        end
        
        function [maxU, maxV, maxSpeed] = maxMushyVel(obj)
            probDomain = obj.problemDomain;
            dx = probDomain.dxCoarse;
            
            numx = (probDomain.domainExtent.hi_i - probDomain.domainExtent.lo_i) + 1;
            numy = (probDomain.domainExtent.hi_j - probDomain.domainExtent.lo_j) + 1;
            
            x = linspace(double(probDomain.domainExtent.lo_i), double(probDomain.domainExtent.hi_i), numx)*dx;
            y = linspace(double(probDomain.domainExtent.lo_j), double(probDomain.domainExtent.hi_j), numy)*dx;
            
            [X,Y] = meshgrid(x, y);
            
            porosity = obj.dataForComp(obj.components.Porosity);
            u = obj.dataForComp(obj.components.xAdvectionvelocity);
            v = obj.dataForComp(obj.components.yAdvectionvelocity);
            
            speed = sqrt(u.^2 + v.^2);
            
            maxU = 0; maxV = 0;
            
            % Remove liquid values
            liquidPorosityLimit = 0.98; % not entirely sure what should be
            u(porosity > liquidPorosityLimit) = 0;
            v(porosity > liquidPorosityLimit) = 0;
            speed(porosity > liquidPorosityLimit) = 0;
            
            maxU = max(max(abs(u)));
            maxV = max(max(abs(v)));
            maxSpeed = max(max(abs(speed)));
        end
        
        function [maxU, maxV, maxSpeed] = maxLiqVel(obj)
            probDomain = obj.problemDomain;
            dx = probDomain.dxCoarse;
            
            numx = (probDomain.domainExtent.hi_i - probDomain.domainExtent.lo_i) + 1;
            numy = (probDomain.domainExtent.hi_j - probDomain.domainExtent.lo_j) + 1;
            
            x = linspace(double(probDomain.domainExtent.lo_i), double(probDomain.domainExtent.hi_i), numx)*dx;
            y = linspace(double(probDomain.domainExtent.lo_j), double(probDomain.domainExtent.hi_j), numy)*dx;
            
            [X,Y] = meshgrid(x, y);
            
            porosity = obj.dataForComp(obj.components.Porosity);
            u = obj.dataForComp(obj.components.xAdvectionvelocity);
            v = obj.dataForComp(obj.components.yAdvectionvelocity);
            
            speed = sqrt(u.^2 + v.^2);
            
            maxU = 0; maxV = 0;
            
            % Remove mushy values
            liquidPorosityLimit = 0.99; % not entirely sure what should be
            u(porosity < liquidPorosityLimit) = 0;
            v(porosity < liquidPorosityLimit) = 0;
            speed(porosity < liquidPorosityLimit) = 0;
            
            maxU = max(max(abs(u)));
            maxV = max(max(abs(v)));
            maxSpeed = max(max(abs(speed)));
        end
        
        
        function [channelDepthPorosity, averagePorosity, depth] = porosityMetrics(obj, porosityQuery)
            probDomain = obj.problemDomain;
            dx = probDomain.dxCoarse;
            
            numx = (probDomain.domainExtent.hi_i - probDomain.domainExtent.lo_i) + 1;
            numy = (probDomain.domainExtent.hi_j - probDomain.domainExtent.lo_j) + 1;
            
            x = linspace(double(probDomain.domainExtent.lo_i), double(probDomain.domainExtent.hi_i), numx)*dx;
            y = linspace(double(probDomain.domainExtent.lo_j), double(probDomain.domainExtent.hi_j), numy)*dx;
            
            [X,Y] = meshgrid(x, y);
            
            porosity = obj.dataForComp(obj.components.Porosity);
            
            %Average porosity in mushy layer is easy
            mushyPorosity = porosity;
            mushyPorosity(mushyPorosity > 0.999) = NaN;
            mushyPorosity(mushyPorosity < 0.001) = NaN;
            averagePorosity = nanmean(nanmean(mushyPorosity));
            
            % Need to find average porosity at depth of channel
            channel_j = 0;
            for channel_j = 1:numy
                maxChi = max(porosity(:, channel_j));
                if maxChi < 0.99
                    % We've found the top of the channel if there are no
                    % liquid regions
                    break;
                end
            end
            
            channelDepthPorosity = mean(mean(porosity(:, channel_j)));
            
            % Now let's find the depth at which mean porosity = porosityQuery
            
            
            % First - find boundary position
            liquidLimit = 0.99;
            [~, liquidy] =  find(porosity > liquidLimit);
            % [mushx, mushy] =  find(porosity < liquidLimit);
            
            for x_i = 1:numx
                chiForx = porosity(x_i, :);
                [~, liquidVals] = find(chiForx > liquidLimit);
                mlBoundaryVals(x_i) = max(liquidVals);
            end
            
            boundary_j = median(mlBoundaryVals);
            
            if porosityQuery < 0
                porosityQuery = averagePorosity;
            end
            
            for j = boundary_j:numy
                chiFory = porosity(:, j);
                avChi = mean(mean(chiFory));
                if avChi < porosityQuery
                    break;
                end
            end
            
            depth = double(j - boundary_j)*dx;
            
            
        end
        
        function [saltSrc] = maxMushyAdvectionSrcTerm(obj)
            probDomain = obj.problemDomain;
            dx = probDomain.dxCoarse;
            
            numx = (probDomain.domainExtent.hi_i - probDomain.domainExtent.lo_i) + 1;
            numy = (probDomain.domainExtent.hi_j - probDomain.domainExtent.lo_j) + 1;
            
            x = linspace(double(probDomain.domainExtent.lo_i), double(probDomain.domainExtent.hi_i), numx)*dx;
            y = linspace(double(probDomain.domainExtent.lo_j), double(probDomain.domainExtent.hi_j), numy)*dx;
            
            [X,Y] = meshgrid(x, y);
            
            salt = abs(obj.dataForComp(obj.components.Saltsrctermgodunov));
            porosity = obj.dataForComp(obj.components.Porosity);
            
            %Average porosity in mushy layer is easy
            %salt(porosity > 0.95) = NaN; % 0.95 because trying to avoid the chimney
            %salt(porosity < 0.7) = NaN; % trying to avoid top of domain
            %saltSrc = max(max(salt));
            saltSrc = median(median(salt)); % trying to ignore massive values near chimney
            
        end
        
        function [channelDepthPerm, averagePerm] = permMetrics(obj)
            probDomain = obj.problemDomain;
            dx = probDomain.dxCoarse;
            
            numx = (probDomain.domainExtent.hi_i - probDomain.domainExtent.lo_i) + 1;
            numy = (probDomain.domainExtent.hi_j - probDomain.domainExtent.lo_j) + 1;
            
            x = linspace(double(probDomain.domainExtent.lo_i), double(probDomain.domainExtent.hi_i), numx)*dx;
            y = linspace(double(probDomain.domainExtent.lo_j), double(probDomain.domainExtent.hi_j), numy)*dx;
            
            [X,Y] = meshgrid(x, y);
            
            permeability = obj.dataForComp(obj.components.Permeability);
            
            %Average porosity in mushy layer is easy
            mushyPermeability = permeability;
            maxPerm = max(max(permeability));
            mushyPermeability(mushyPermeability > 0.99*maxPerm) = NaN;
            mushyPermeability(mushyPermeability < 0.001) = NaN;
            averagePerm = nanmean(nanmean(mushyPermeability));
            
            % Need to find average porosity at depth of channel
            channel_j = 0;
            for channel_j = 1:numy
                maxChi = max(permeability(:, channel_j));
                if maxChi < 0.99*maxPerm
                    % We've found the top of the channel if there are no
                    % liquid regions
                    break;
                end
            end
            
            channelDepthPerm = mean(mean(permeability(:, channel_j)));
            
        end
        
        
        
        function flux = computeVerticalSoluteFlux(obj, frameAdv, Le, katzUnits)
            if nargin < 4
                katzUnits = true;
            end
            
            %T = obj.dataForComp(obj.components.Temperature).';
            chi = obj.dataForComp(obj.components.Porosity).';
            
            S = obj.dataForComp(obj.components.Bulkconcentration).';
            Sl = obj.dataForComp(obj.components.Liquidconcentration).';
            
            V = obj.dataForComp(obj.components.yAdvectionvelocity).';
            
            
            % Do some conversion to Wells Units
            %CR = CR - 1; % Always do this, as CR comes from the folder name
            
            if katzUnits
                
                % H = H - 1;
                %T = T - 1;
                S = S + 1;
                Sl = Sl + 1;
            end
            
            [X, Y] = obj.grid();
            dx = X(1, 2) - X(1, 1);
            
            flux = computeSaltFlux(S, V, Sl, chi, dx, frameAdv, Le);
            
            
        end
        
    end
    
    
    
end
