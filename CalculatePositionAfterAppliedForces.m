function [Filaments, Membrane, Adhesions, FAConnections, Tensions, kBreaks, Data] = ...
                        CalculatePositionAfterAppliedForces(Filaments,Membrane,Adhesions,FAConnections,ModelParameters)

    % Find all the main filaments of connected filament structures (group of attached filaments)
    idxMF = find(Filaments.Parent == 0);
    nMF = length(idxMF); % Total number of filament structures
    dt = ModelParameters.TimeStep;
    nS = size(Membrane.Segments,1); % Number of membrane segments
    FilamentForce = ModelParameters.FilamentForce; %1; % Generic filament force along filament direction (pN)
    k  = ModelParameters.MembraneSpringConstant;
    D  = ModelParameters.SpringWidth;
    
    % Go through each filament structure and test each of its filaments to
    % see if it is hitting the membrane and calculate the total force acting on the structure.
    Tensions = [];
    kBreaks = [];
    Data.Adhesions = Adhesions;
    Data.FAConnections.Tensions = [];
    Data.FAConnections.Regions  = [];
    Data.FAConnections.FilmentName = [];
    Data.FAConnections.AdhesionIdx = [];
    Data.FilamentTips.YSpeed    = [];
    Data.FilamentTips.Region    = [];
    Data.FilamentTips.FilamentName = [];
    Data.FilamentTips.XYPosition = [];
    
    % FILAMENT SECTION
    if ~isempty(Filaments.Name)
        for MF = 1:nMF
            Offset = 0;
            L = 0; % Total length of Filament
            idx1 = find( Filaments.MainIndex == Filaments.Name(idxMF(MF)) ); % Find all filaments attached to same structure
            % Calculate total length of filament structure
                for f = idx1';  L = L + length(Filaments.MonomerIndices{f}); end
            % CALCULATE FORCES DUE TO FILAMENT-ADHESION CONNECTIONS and break connections where force is greater than threshold
               [FAx,FAy,Adhesions,FAConnections,Ts,kBs,Regions,AttachedFilament,AdhesionIndex] = CalculateForceDueToFAConnections(idx1,Filaments,Adhesions,FAConnections,ModelParameters);
                % Concatentate all tensions and data about adhesions attached to filaments in this structure.
                Tensions = [Tensions;Ts];
                kBreaks  = [kBreaks; kBs];
                Data.FAConnections.Tensions = [Data.FAConnections.Tensions; Ts];
                Data.FAConnections.Regions  = [Data.FAConnections.Regions; Regions];
                Data.FAConnections.FilmentName = [Data.FAConnections.FilmentName; AttachedFilament];
                Data.FAConnections.AdhesionIdx = [Data.FAConnections.AdhesionIdx; AdhesionIndex];
            % CALCULATE FORCES ACTING ON FILAMENT FROM MEMBRANE
            % first check if and which filaments of the current structure are hitting the membrane ---------------
                [FilInd,MemInd] = FindWhichFilamentsAreHittingMembraneAndWhere(idx1,Filaments,Membrane,ModelParameters);
                % Calculate Membrane Spring and Boundary Force on Filamentary structure
                Fm = TotalForceExertedByMembraneOnFilamentStructure(FilInd,MemInd,Filaments,Membrane,k,D,nS,ModelParameters);
            % Create Containment Force pointing at [0,0] from both sides in the X-dimension
                Fc = 0;% Fc = CreateContainmentForce(Filaments,idx1,ModelParameters);
            % CALCULATE NEW POSITION FOR FILAMENT STRUCTURE
                gamma1 = (4*pi*3.01*L)/(0.84+log(L/0.007))/1000; % gamma for random fluid motion (pN*S/nm) 
                SD = sqrt(2*4.114*gamma1/dt); % sqrt(2*(pN*nm)*(pN*S/nm)/S)
                Fx = SD*randn(1);
                Fy = SD*randn(1);
                
                for f = idx1' % Loop through all the filaments attached to this structure
                    yprevious = Filaments.XYCoords{f}(end,2);
                    
                    % START Calculate New Filament Position --------------------------------------------------------------------------
                    Filaments.XYCoords{f}(:,1) = Filaments.XYCoords{f}(:,1) + (Fx/gamma1)*dt + (Fc/gamma1)*dt + (FAx/gamma1)*dt; % Compute new position in X direction
                    Filaments.XYCoords{f}(:,2) = Filaments.XYCoords{f}(:,2) + (Fy/gamma1)*dt + (Fm/gamma1)*dt + (FAy/gamma1)*dt; % Compute new position in Y direction
                    % END Calculate New Filament Position ----------------------------------------------------------------------------
                    
                    % Record Filament Speed, Position, Region, etc of each tip .................................................
                    R = FindWhichRegionTheFilamentIsIn(Filaments.XYCoords{f}(end,1),Membrane,ModelParameters);
                    Data.FilamentTips.YSpeed       = [Data.FilamentTips.YSpeed; (Filaments.XYCoords{f}(end,2) - yprevious)/ModelParameters.TimeStep];
                    Data.FilamentTips.Region       = [Data.FilamentTips.Region; R];
                    Data.FilamentTips.FilamentName = [Data.FilamentTips.FilamentName; Filaments.Name(f,1)];
                    Data.FilamentTips.XYPosition   = [Data.FilamentTips.XYPosition; Filaments.XYCoords{f}(end,:)];
                    
                    % Mirror Filaments horizontally if their tip moves outside the boundary edge ------------------
                    MemWidth = Membrane.Nodes(end,1) - Membrane.Nodes(1,1) - 40;
                    BreakIdx = [];
                    
                    if  Filaments.XYCoords{f}(end,1) < Membrane.Nodes(1,1)
                        Offset = +MemWidth; 
                        BreakIdx = find(FAConnections.FilamentName == Filaments.Name(f,1)); % Break adhesions connection if it exists
                    elseif Filaments.XYCoords{f}(end,1) > Membrane.Nodes(end,1)
                        Offset = -MemWidth;
                        BreakIdx = find(FAConnections.FilamentName == Filaments.Name(f,1)); % Break adhesions connection if it exists
                    end
                    
                    if ~isempty(BreakIdx)
                        FAConnections.AdhesionIndex(BreakIdx,:) = [];
                        FAConnections.FilamentName(BreakIdx,:)  = [];
                        FAConnections.MonomerIndex(BreakIdx,:)  = [];
                    end
                    %----------------------------------------------------------------------------------------------    
                end
                
                % Apply horizontal offset (if it applies) to all filaments in structure
                if Offset ~= 0 
                    for f = idx1'
                        Filaments.XYCoords{f}(:,1) = Filaments.XYCoords{f}(:,1) + Offset;    
                    end
                end
                
        end
    end
    
    % Set limits for membrane and calculate new positions -------------------------------
        if ModelParameters.BoundaryFixed
            MemStart = 2;
            MemEnd = nS-1;
        else
            MemStart = 1;
            MemEnd = nS;
        end

        gamma3 = ModelParameters.MembraneGamma; % gamma1 % membrane gamma
        % MEMBRANE SECTION STEP 1 (move membrane based on filament positions)
        if ~isempty(Filaments.Name)
            FilamentTips = GetFilamentsTipLocations(Filaments);
            MembraneOld = Membrane;
            for s = MemStart:MemEnd
                FF = CalculateForceFromFilamentsHittingMembraneSegment(FilamentForce,Filaments,FilamentTips,MembraneOld,s,ModelParameters);
                Membrane.Nodes(Membrane.Segments(s,:),2) = MembraneOld.Nodes(MembraneOld.Segments(s,:),2) + (FF/gamma3)*dt; % Move y-position of segment nodes
            end
        end
        % MEMBRANE SECTION STEP 2 (move membrane based on spring tensions of neighboring segments)
        MembraneOld = Membrane;
        for s = MemStart:MemEnd
            FN = CalculateForceFromNeighboringMembraneSegments(MembraneOld,s,nS,k,D);
            Membrane.Nodes(Membrane.Segments(s,:),2) = MembraneOld.Nodes(MembraneOld.Segments(s,:),2) + (FN/gamma3)*dt; % Move y-position of segment nodes
        end
    %-------------------------------------------------------------------------------
    
end

%======================================================================================================
%======================================================================================================

function Fm = TotalForceExertedByMembraneOnFilamentStructure(FilInd,MemInd,Filaments,Membrane,k,D,nS,ModelParameters)

        Fm = 0;
        % Calculate membrane segment spring force
        if ~isempty(MemInd)
            for n = 1:length(MemInd)
                    m = MemInd(n,1);
                    f = FilInd(n,1);
                    % Initalize spring force on left and right side of segment as zero
                    Fy1 = 0;
                    Fy2 = 0;
                    % Calculate spring force from spring on left side of membrane segment
                    if isequal(m,1) % Don't calcualte if this is the far left membrane segment
                        Fdir1 = sign(Membrane.Nodes(Membrane.Segments(nS,2),2) - Membrane.Nodes(Membrane.Segments(m,1),2));
                        delX1 =  abs(Membrane.Nodes(Membrane.Segments(nS,2),1) - Membrane.Nodes(Membrane.Segments(m,1),1));
                        delY1 =  abs(Membrane.Nodes(Membrane.Segments(nS,2),2) - Membrane.Nodes(Membrane.Segments(m,1),2));
                        Fs1 = k*(sqrt(delX1^2 + delY1^2)-D); % Total spring force
                        Fy1 = Fdir1*Fs1*sin(atan(delY1/delX1));  % Y component of spring force
                    else
                        Fdir1 = sign(Membrane.Nodes(Membrane.Segments(m-1,2),2) - Membrane.Nodes(Membrane.Segments(m,1),2));
                        delX1 =  abs(Membrane.Nodes(Membrane.Segments(m-1,2),1) - Membrane.Nodes(Membrane.Segments(m,1),1));
                        delY1 =  abs(Membrane.Nodes(Membrane.Segments(m-1,2),2) - Membrane.Nodes(Membrane.Segments(m,1),2));
                        Fs1 = k*(sqrt(delX1^2 + delY1^2)-D); % Total spring force
                        Fy1 = Fdir1*Fs1*sin(atan(delY1/delX1));  % Y component of spring force
                    end
                    % Calculate spring force from spring on right side of membrane segment
                    if isequal(m,nS) % Don't calcualte if this is the far left membrane segment
                        Fdir2 = sign(Membrane.Nodes(Membrane.Segments(1,1),2) - Membrane.Nodes(Membrane.Segments(m,2),2));
                        delX2 =  abs(Membrane.Nodes(Membrane.Segments(m,2),1) - Membrane.Nodes(Membrane.Segments(1,1),1));
                        delY2 =  abs(Membrane.Nodes(Membrane.Segments(m,2),2) - Membrane.Nodes(Membrane.Segments(1,1),2));
                        Fs2 = k*(sqrt(delX2^2+delY2^2)-D); % Total spring force
                        Fy2 = Fdir2*Fs2*sin(atan(delY2/delX2));  % Y component of spring force
                    else
                        Fdir2 = sign(Membrane.Nodes(Membrane.Segments(m+1,1),2) - Membrane.Nodes(Membrane.Segments(m,2),2));
                        delX2 =  abs(Membrane.Nodes(Membrane.Segments(m,2),1) - Membrane.Nodes(Membrane.Segments(m+1,1),1));
                        delY2 =  abs(Membrane.Nodes(Membrane.Segments(m,2),2) - Membrane.Nodes(Membrane.Segments(m+1,1),2));
                        Fs2 = k*(sqrt(delX2^2+delY2^2)-D); % Total spring force
                        Fy2 = Fdir2*Fs2*sin(atan(delY2/delX2));  % Y component of spring force
                    end
                    % Calculate Boundary Force from filament crossing membrane
                    CrossDelta = Filaments.XYCoords{f}(end,2) - Membrane.Nodes(Membrane.Segments(m,1),2); % How much filament crosses membrane segment
                    CrossDelta(CrossDelta < 0) = 0;
                    Fb = -ModelParameters.BoundaryForceSpringConstant*CrossDelta;

                    if (Fy1+Fy2) < 0 % Only add the contribution if the total force of the two springs points points towards the filament
                        Fm = Fm + Fy1 + Fy2 + Fb;
                    end
            end
        end
end

%======================================================================================================
%======================================================================================================

function [FilInd,MemInd] = FindWhichFilamentsAreHittingMembraneAndWhere(idx,Filaments,Membrane,ModelParameters)

        FilInd = [];
        MemInd = [];
        SpringWidth = ModelParameters.SpringWidth;
        nS = size(Membrane.Segments,1);
         
        for f = idx' % For each filament of a filament structure
            Xend = Filaments.XYCoords{f}(end,1); % X-coordinate of filament tip
            Yend = Filaments.XYCoords{f}(end,2); % Y-coordinate of filament tip
            for m = 1:nS % For each membrane segment
                   Yseg = Membrane.Nodes(Membrane.Segments(m,1),2); %  y-coordinate current Membrane segment
                   % Test if filament tip is between membrane segment endpoints and is above membrane's y-position
                   if Xend >= Membrane.Nodes(Membrane.Segments(m,1),1)-SpringWidth/2 &&... % Filament X-end >= Left segment endpoint
                      Xend <= Membrane.Nodes(Membrane.Segments(m,2),1)+SpringWidth/2   % Filament X-end <= Right segment endpoint
                          if Yend > (Yseg-ModelParameters.ContactThresholdForMembrane) % Filament tip > Membrane segment     
                                % Record this filament and the membrane it is hitting
                                FilInd = [FilInd; f];
                                MemInd = [MemInd; m];
                          end
                   end   
                end
        end

end

%======================================================================================================
%======================================================================================================

function Region = FindWhichRegionTheFilamentIsIn(Xpt,Membrane,ModelParameters)

        % Xend = X coordinate of filament tip
        SpringWidth = ModelParameters.SpringWidth;
        Region = find((Membrane.Nodes(Membrane.Segments(:,1),1)-SpringWidth/2) <= Xpt &... % Filament X-end >= Left segment endpoint
                      (Membrane.Nodes(Membrane.Segments(:,2),1)+SpringWidth/2) >= Xpt);   
                     
        if isempty(Region)
            Region = NaN;
        end

end

%======================================================================================================
%======================================================================================================

% function FilamentTips = GetFilamentsTipLocations(Filaments)
%     
%         nF = size(Filaments.XYCoords,1);
%         FilamentTips = zeros(nF,2);
%         for f = 1:nF
%             FilamentTips(f,:) = Filaments.XYCoords{f}(end,:);
%         end
% 
% end

%======================================================================================================
%======================================================================================================

function FF = CalculateForceFromFilamentsHittingMembraneSegment(FilamentForce,Filaments,FilamentTips,Membrane,s,ModelParameters)
        
        FF = 0;
        NodeLeft    = Membrane.Nodes(Membrane.Segments(s,1),:);
        NodeRight   = Membrane.Nodes(Membrane.Segments(s,2),:);
        SegmentYPos = Membrane.Nodes(Membrane.Segments(s,1),2);
        SpringWidth = ModelParameters.SpringWidth;
        
        idx = find( FilamentTips(:,1) >= NodeLeft(:,1)  - SpringWidth/2 & ...
                    FilamentTips(:,1) <= NodeRight(:,1) + SpringWidth/2 & ...                                     
                    FilamentTips(:,2) >  SegmentYPos );
                
        % Add up all the y-components of all the filament forces on this membrane segment        
        for f = idx'
            FF = FF + abs(FilamentForce*Filaments.UnitVector(f,2)); % Get Y component of filament force using y component of filament unit vector
        end   
        
end

%======================================================================================================
%======================================================================================================

function FN = CalculateForceFromNeighboringMembraneSegments(Membrane,s,nS,k,D)

        FN = 0;
        Fy1 = 0;
        Fy2 = 0;
        % Calcualte spring force from spring on left side of membrane segment
        if isequal(s,1) % Don't calcualte if this is the far left membrane segment
            Fdir1 = -sign(Membrane.Nodes(Membrane.Segments(nS,2),2) - Membrane.Nodes(Membrane.Segments(s,1),2));
            delX1 =   abs(Membrane.Nodes(Membrane.Segments(nS,2),1) - Membrane.Nodes(Membrane.Segments(s,1),1));
            delY1 =   abs(Membrane.Nodes(Membrane.Segments(nS,2),2) - Membrane.Nodes(Membrane.Segments(s,1),2));
            Fs1 = k*(D-sqrt(delX1^2 + delY1^2)); % Total spring force
            Fy1 = Fdir1*Fs1*sin( atan(delY1/delX1) );  % Y component of spring force
        else
            Fdir1 = -sign(Membrane.Nodes(Membrane.Segments(s-1,2),2) - Membrane.Nodes(Membrane.Segments(s,1),2));
            delX1 =   abs(Membrane.Nodes(Membrane.Segments(s-1,2),1) - Membrane.Nodes(Membrane.Segments(s,1),1));
            delY1 =   abs(Membrane.Nodes(Membrane.Segments(s-1,2),2) - Membrane.Nodes(Membrane.Segments(s,1),2));
            Fs1 = k*(D-sqrt(delX1^2 + delY1^2)); % Total spring force
            Fy1 = Fdir1*Fs1*sin( atan(delY1/delX1) );  % Y component of spring force
        end
        % Calcualte spring force from spring on right side of membrane segment
        if isequal(s,nS) % Don't calcualte if this is the far left membrane segment
            Fdir2 = -sign(Membrane.Nodes(Membrane.Segments(1,1),2) - Membrane.Nodes(Membrane.Segments(s,2),2));
            delX2 =   abs(Membrane.Nodes(Membrane.Segments(s,2),1) - Membrane.Nodes(Membrane.Segments(1,1),1));
            delY2 =   abs(Membrane.Nodes(Membrane.Segments(s,2),2) - Membrane.Nodes(Membrane.Segments(1,1),2));
            Fs2 = k*(D-sqrt(delX2^2 + delY2^2)); % Total spring force
            Fy2 = Fdir2*Fs2*sin( atan(delY2/delX2) );  % Y component of spring force
        else
            Fdir2 = -sign(Membrane.Nodes(Membrane.Segments(s+1,1),2) - Membrane.Nodes(Membrane.Segments(s,2),2));
            delX2 =   abs(Membrane.Nodes(Membrane.Segments(s,2),1) - Membrane.Nodes(Membrane.Segments(s+1,1),1));
            delY2 =   abs(Membrane.Nodes(Membrane.Segments(s,2),2) - Membrane.Nodes(Membrane.Segments(s+1,1),2));
            Fs2 = k*(D-sqrt(delX2^2 + delY2^2)); % Total spring force
            Fy2 = Fdir2*Fs2*sin( atan(delY2/delX2) );  % Y component of spring force
        end
        FN = FN + Fy1 + Fy2;
end

%======================================================================================================
%======================================================================================================

function Fc = CreateContainmentForce(Filaments,idx1,ModelParameters)
    
    % NOT USED ------
    nF = length(idx1);
    Xtip = zeros(length(idx1),1);
    xthresh = 900;
    
    for n = 1:nF
        Xtip(n,1) = Filaments.XYCoords{idx1(n,1)}(end,1);
    end
    
    X = ModelParameters.MemleadingEdgeLengthInNanometers/2;
    [~,idx2] = max(abs(Xtip));
    Fc = -X./(abs(Xtip(idx2)-xthresh).^2) + X./(abs(Xtip(idx2)+xthresh).^2);
    Fc(Xtip(idx2) < -xthresh + 1) =  1000;
    Fc(Xtip(idx2) >  xthresh - 1) = -1000;
    
%      figure(1)
%      x  = linspace(-1000,1000,20001)';
%     Fc = -10./(abs(x-xthresh).^0.5) + 10./(abs(x+xthresh).^0.5);
%     Fc(x <= -xthresh + 1) =  10;
%     Fc(x >=  xthresh - 1) = -10;
%      plot(x,Fc,'.-')
%     
end

%======================================================================================================
%======================================================================================================

function [FAx,FAy,Adhesions,FAConnections,Tensions,kBreaks,Regions,AttachedFilament,AdhesionIndex] = ...
                                   CalculateForceDueToFAConnections(idx,Filaments,Adhesions,FAConnections,ModelParameters)

        slope  = ModelParameters.DoubleExp_Slope;
        base   = ModelParameters.DoubleExp_Base;
        slope2 = ModelParameters.DoubleExp_Slope2;
        A = ModelParameters.DoubleExp_A;
        B = ModelParameters.DoubleExp_B;
        C = ModelParameters.DoubleExp_C;
        Tensions = [];
        kBreaks = [];
        AttachedFilament = [];
        AdhesionIndex    = [];
        FAx = 0;
        FAy = 0;
        BreakIdx = [];
        Regions = [];
        
        if ~isempty(FAConnections.AdhesionIndex)
            for f = idx'
                 con = find(FAConnections.FilamentName == Filaments.Name(f,1));
                 if ~isempty(con)
                         for c = con'
                             midx = find(Filaments.MonomerIndices{f} == FAConnections.MonomerIndex(c));
                             aidx = FAConnections.AdhesionIndex(c,1);
                             xDist = Adhesions.XYPoints(aidx,1) - Filaments.XYCoords{f}(midx,1); 
                             yDist = Adhesions.XYPoints(aidx,2) - Filaments.XYCoords{f}(midx,2); 
                             StretchDist = sqrt(xDist^2 + yDist^2) - ModelParameters.AdhesionSpringEqLength;
                             StretchDist(StretchDist < 0) = 0;
                             ConnectionTension = StretchDist*ModelParameters.AdhesionSpringConstant;
                             Tensions = [Tensions;ConnectionTension];
                             AttachedFilament = [AttachedFilament; f];
                             AdhesionIndex = [AdhesionIndex; aidx];
                             ForceX = ConnectionTension*(xDist/sqrt(xDist^2+yDist^2));
                             ForceY = ConnectionTension*(yDist/sqrt(xDist^2+yDist^2));
                             FAx = FAx + ForceX;
                             FAy = FAy + ForceY;
                             kBreakRate =  base*exp(A*slope*ConnectionTension) + C*exp(B*slope2*ConnectionTension);
                             %kBreakRate = 0.1;
                             kBreaks = [kBreaks; kBreakRate];
                             % Test if Adhesion-Filament connection breaks (Molecular Clutch)
                             if rand(1) < (kBreakRate)*ModelParameters.TimeStep && ModelParameters.AllowAdhesionFilamentBondBreak
                                 BreakIdx = [BreakIdx; c];
                             end
                             Regions = [Regions; Adhesions.RegionLocation(aidx,1)];
                         end
                 end
            end
            
            % Remove Adhesion-Filament broken connections -----------------
            if ~isempty(BreakIdx)
                 a = FAConnections.AdhesionIndex(BreakIdx,1);
                 Adhesions.AttachedFilamentName(a,1) = NaN;
                 FAConnections.AdhesionIndex(BreakIdx,:) = [];
                 FAConnections.FilamentName(BreakIdx,:) = [];
                 FAConnections.MonomerIndex(BreakIdx,:) = [];
                 if ModelParameters.Adhesion_OFF_DependentOnAdhesionFilamentTension
                        Adhesions.XYPoints(a,1) = NaN;
                        Adhesions.XYPoints(a,2) = NaN;
                        Adhesions.RegionLocation(a,:) = NaN;
                        Adhesions.ActiveStatus(a,1) = false;
                        Adhesions.AttachedFilamentName(a,1) = NaN;
                 end
            end 
            %--------------------------------------------------------------
        end


end

