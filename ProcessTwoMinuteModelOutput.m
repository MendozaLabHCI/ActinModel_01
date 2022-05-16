function DATA = ProcessTwoMinuteModelOutput(SimData)

%-----------------------------------------------------------------------------
%-----------------------------------------------------------------------------
% Protrusion Speed vs. Time --------------------------------------------------

    dt = diff(SimData.TimeVec(1,1:2));
    nT = size(SimData.TimeVec,2);
    
    nodes = SimData.Segments(:,1); % Use the left node position of each membrane segment to calculate velocity
    % Velocity for each segment and timepoint:
    MembraneVelocity = squeeze( (SimData.Nodes(nodes,2,2:end) - SimData.Nodes(nodes,2,1:end-1))/dt ); 
    
    % Resample Time Vector and take averages of Membrane Velocity for every second
    nTh = round(1/dt);
    RSV = 1:nTh:nT; % ReSampleVector
    Time = SimData.TimeVec(1,RSV);
    
    % Calculate the mean values over 1 second intervals
    SampledMembraneVelocity = zeros(size(MembraneVelocity,1),size(Time,2)-1);
    for n = 1:(size(RSV,2)-1)
        indices = RSV(1,n):(RSV(1,n)+nTh-1);
        SampledMembraneVelocity(:,n) = mean( MembraneVelocity(:,indices),2 );
    end
    
%-----------------------------------------------------------------------------
%-----------------------------------------------------------------------------
% Fil-Adhesion Tension vs. Time

    nSeg =  size(SimData.Segments,1);
    N = size(SimData.DATA,1);
    
    % Accumulate values to calculate 95th percentile ---------------------------------- 
%     tensions = [];
%     for n = 2:N % For each timepoint n (model interation)
%         tensions = [tensions; SimData.DATA{n,1}.FAConnections.Tensions];
%     end
%    %----------------------------------------------------------------------------------
%    lim1 = prctile(tensions,[2.5,97.5]); % Calculate median centered 95th percentile limits
    
    FAtensions = nan(nSeg,N-1);
    for n = 2:N % For each timepoint n (model interation)
        TempValues = cell(nSeg,1); % Reset this to a blank cell array at each time point
        for m = 1:length(SimData.DATA{n,1}.FAConnections.Tensions) % Go through each tension at timepoint n, and sort it into it's appropriate region
            Tension = double( SimData.DATA{n,1}.FAConnections.Tensions(m,1) );
            R  = SimData.DATA{n,1}.FAConnections.Regions(m,1);
            if ~isnan(R) % && Tension > lim1(1) && Tension < lim1(2)   % && tension < 500
                TempValues{R,1} = [TempValues{R,1}; Tension]; % accumulate values for each regions at timepoint n.
            end
        end
        FAtensions(:,n-1) = cellfun(@mean,TempValues); % Calculate mean values for each region at timepoint n.
    end

    % Calculate the mean values over 1 second intervals
    SampledFAtensions = nan(size(FAtensions,1), size(Time,2)-1);
    for n = 1:(size(RSV,2)-1)
        indices = RSV(1,n):(RSV(1,n)+nTh-1);
        SampledFAtensions(:,n) = mean(FAtensions(:,indices),2,'omitnan');
    end
    %SampledFAtensions(isnan(SampledFAtensions)) = 0;
    %SampledFAtensions = ReplaceNANwithLocalMean(SampledFAtensions);
    
%-----------------------------------------------------------------------------
%-----------------------------------------------------------------------------
% Retro Speed  vs. Time 

%     retrospeeds = [];
%     for n = 2:N % For each timepoint n (model interation)
%         retrospeeds = [retrospeeds; SimData.DATA{n,1}.FilamentTips.YSpeed];
%     end  
%     %----------------------------------------------------------------------------------
%     lim2 = prctile(retrospeeds,[2.5,97.5]); % Calculate median centered 95th percentile limits

    RetroSpeeds = nan(nSeg,N-1);
    
    for n = 2:N % For each timepoint n (model interation)
        TempValues = cell(nSeg,1); % Reset this to a blank cell array at each time point
        for m = 1:length(SimData.DATA{n}.FilamentTips.YSpeed) % Go through each retro speed at timepoint n, and sort it into it's appropriate region
            Rspeed =  double( SimData.DATA{n}.FilamentTips.YSpeed(m,1) );
            R =      SimData.DATA{n}.FilamentTips.Region(m,1);
            if ~isnan(R) && ~isequal(R,0)% && Rspeed > lim2(1) && Rspeed < lim2(2) %abs(speed) < 1000
                TempValues{R,1} = [TempValues{R,1}; Rspeed]; % accumulate values for each regions at timepoint n.
            end
        end
        RetroSpeeds(:,n-1) = cellfun(@mean,TempValues); % Calculate median values for each region at timepoint n.
    end
   
    % Calculate the mean values over 1 second intervals
    SampledRetroSpeeds = nan(size(RetroSpeeds,1),size(Time,2)-1);
    for n = 1:(size(RSV,2)-1)
        indices = RSV(1,n):(RSV(1,n)+nTh-1);
        SampledRetroSpeeds(:,n) = mean(RetroSpeeds(:,indices),2,'omitnan');
    end
    %SampledRetroSpeeds = ReplaceNANwithLocalMean(SampledRetroSpeeds);
    

    
% SAVE -------------------------------------------------------------------
    DATA.MembraneVelocity = single(MembraneVelocity);
    DATA.SampledMembraneVelocity = single(SampledMembraneVelocity);
    DATA.FAtensions = single(FAtensions);
    DATA.SampledFAtensions = SampledFAtensions;
    DATA.RetroSpeeds = single(RetroSpeeds);
    DATA.SampledRetroSpeeds = single(SampledRetroSpeeds);

end