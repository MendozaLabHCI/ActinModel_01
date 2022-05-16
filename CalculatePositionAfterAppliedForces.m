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
    
    Data.FAConnections.Tensions = [];
    Data.FAConnections.Regions  = [];
    Data.FAConnections.FilamentName = [];
    Data.FAConnections.AdhesionIdx = [];
    Data.FilamentTips.YSpeed    = [];
    Data.FilamentTips.Region    = [];
    Data.FilamentTips.FilamentName = [];
    Data.FilamentTips.XYPosition = [];
    Data.FilamentTips.StructName = [];
    
    % FILAMENT SECTION
    if ~isempty(Filaments.Name)
        for MF = 1:nMF % For each filament structure
            L = 0; % Total length of Filament
            idx1 = find( Filaments.MainIndex == Filaments.Name(idxMF(MF)) ); % Find all filaments attached to same structure
            nF = length(idx1); % number of filaments in this structure
            % Calculate total length of filament structure
                for f = idx1'  
                    L = L + length(Filaments.MonomerIndices{f})*ModelParameters.MonomerLength/1000; % L is total length of filament structure in microns
                end 
            % CALCULATE FORCES DUE TO FILAMENT-ADHESION CONNECTIONS and break connections where force is greater than threshold
               [FAx,FAy,Adhesions,FAConnections,Ts,kBs,Regions,AttachedFilament,AdhesionIndex] = CalculateForceDueToFAConnections(idx1,Filaments,Adhesions,FAConnections,ModelParameters);
                % Concatentate all tensions and data about adhesions attached to filaments in this structure.
                Tensions = [Tensions;Ts];
                kBreaks  = [kBreaks; kBs];
                Data.FAConnections.Tensions = [Data.FAConnections.Tensions; single(Ts)];
                Data.FAConnections.Regions  = [Data.FAConnections.Regions; single(Regions)];
                Data.FAConnections.FilamentName = [Data.FAConnections.FilamentName; single(AttachedFilament)];
                Data.FAConnections.AdhesionIdx = [Data.FAConnections.AdhesionIdx; single(AdhesionIndex)];
            % CALCULATE FORCES ACTING ON FILAMENT FROM MEMBRANE
            % first check if and which filaments of the current structure are hitting the membrane ---------------
                [FilInd,MemInd] = FindWhichFilamentsAreHittingMembraneAndWhere(idx1,Filaments,Membrane,ModelParameters); %Filaments.XYCoords{f}(end,:)
                % Calculate Membrane Spring and Boundary Force on Filamentary structure
                Fm = TotalForceExertedByMembraneOnFilamentStructure(FilInd,MemInd,Filaments,Membrane,k,D,nS,ModelParameters);
            % Create Containment Force pointing at [0,0] from both sides in the X-dimension
                Fc = 0; % Fc = CreateContainmentForce(Filaments,idx1,ModelParameters);
            % CALCULATE NEW POSITION FOR FILAMENT STRUCTURE
                Nu = ModelParameters.CytoplasmViscosity;
                gamma1 = ( (4*pi*Nu*L) / (0.84+log(L/0.007)) )/1000; %  gamma for random fluid motion (pN*s/nm) 
                SD = sqrt(2*4.114*gamma1/dt); % sqrt(2*(pN*nm)*(pN*S/nm)/S)
                Fx = SD*randn(1);
                Fy = SD*randn(1);
                
                tipXpositions = NaN(nF,1);
                for n = 1:nF % Loop through all the filaments attached to this structure
                    f = idx1(n);
                    yprevious = Filaments.XYCoords{f}(end,2);
                    
                    % START Calculate New Filament Position --------------------------------------------------------------------------
                    Filaments.XYCoords{f}(:,1) = Filaments.XYCoords{f}(:,1) + (Fx/gamma1)*dt + (Fc/gamma1)*dt + (FAx/gamma1)*dt; % Compute new position in X direction
                    Filaments.XYCoords{f}(:,2) = Filaments.XYCoords{f}(:,2) + (Fy/gamma1)*dt + (Fm/gamma1)*dt + (FAy/gamma1)*dt; % Compute new position in Y direction
                    % END Calculate New Filament Position ----------------------------------------------------------------------------
                    
                    % Record Filament Speed, Position, Region, etc of each tip .................................................
                    R = FindWhichRegionTheFilamentIsIn(Filaments.XYCoords{f}(end,1),Membrane,ModelParameters);
                    Data.FilamentTips.YSpeed       = [Data.FilamentTips.YSpeed; single((Filaments.XYCoords{f}(end,2) - yprevious)/ModelParameters.TimeStep)];
                    Data.FilamentTips.Region       = [Data.FilamentTips.Region; uint16(R)];
                    Data.FilamentTips.FilamentName = [Data.FilamentTips.FilamentName; single(Filaments.Name(f,1))];
                    Data.FilamentTips.XYPosition   = [Data.FilamentTips.XYPosition;   single(Filaments.XYCoords{f}(end,:))];
                    Data.FilamentTips.StructName   = [Data.FilamentTips.StructName;   single(Filaments.MainIndex(f,1))];
                    
                    tipXpositions(n,1) = Filaments.XYCoords{f}(end,1); % Record for mirroring
                end
                
                % Mirror Filaments horizontally if their tip moves outside the boundary edge ------------------
                BreakIdx = [];
                Offset = 0;

                % Check if any of the filament tips in the structure are crossing the left or right membrane edge
                if  min(tipXpositions) < Membrane.Nodes(1,1) 
                    Offset = Membrane.Nodes(end,1) - max(tipXpositions); 
                elseif max(tipXpositions) > Membrane.Nodes(end,1)
                    Offset = Membrane.Nodes(1,1) - min(tipXpositions);
                end
                
                % If there are any filaments out of bounds, move to the otherside just within the bounds
                if ~isequal(Offset,0)
                    % Apply offset to each filament in the structure
                    for n = 1:nF
                        f = idx1(n);
                        Filaments.XYCoords{f}(:,1) = Filaments.XYCoords{f}(:,1) + Offset; % Apply Offset
                        BreakIdx = [BreakIdx; find(FAConnections.FilamentName == Filaments.Name(f,1))]; % Find adhesion connections if they exists
                    end
                    % Remove any adhesion connections. 
                    if ~isempty(BreakIdx)
                         a = FAConnections.AdhesionIndex(BreakIdx,1);
                         Adhesions.AttachedFilamentName(a,1) = NaN;
                         FAConnections.AdhesionIndex(BreakIdx,:) = [];
                         FAConnections.FilamentName(BreakIdx,:) = [];
                         FAConnections.MonomerIndex(BreakIdx,:) = [];
                         % If Molecular clutch is ON deactivate adhesions
                         if ModelParameters.Adhesion_MolecularClutchOn 
                             Adhesions.XYPoints(a,1) = NaN;
                             Adhesions.XYPoints(a,2) = NaN;
                             Adhesions.RegionLocation(a,:) = NaN;
                             Adhesions.ActiveStatus(a,1) = false;
                             Adhesions.AttachedFilamentName(a,1) = NaN;
                         end
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
                    
                    Fsprings = Fy1 + Fy2;
                    Fsprings(Fsprings > 0) = 0; 
                    % Calculate Boundary Force from filament crossing membrane
                    CrossDelta = Filaments.XYCoords{f}(end,2) - (Membrane.Nodes(Membrane.Segments(m,1),2)); % How much filament crosses membrane segment
                    
                    Fb = -ModelParameters.BoundaryForceSpringConstant*CrossDelta;
                    Fb(CrossDelta < 0) = 0;
                    %if (Fy1+Fy2) < 0 % Only add the contribution if the total force of the two springs points points towards the filament
                    Fm = Fm + Fsprings + Fb;
                   % end
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
% NOT USED
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


        [slope1,base,slope2,A,B,C] = MolecularClutchPeakParameters(ModelParameters.MolecularClutch_PeakNumber);
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
                 con = find(FAConnections.FilamentName == Filaments.Name(f,1)); % Find adhesions attached to this filament
                 if ~isempty(con)
                         for c = con'
                             midx = find( Filaments.MonomerIndices{f} == FAConnections.MonomerIndex(c) ); % Find index of attached monomer (used to get XY coords of monomer)
                             aidx = FAConnections.AdhesionIndex(c,1);                                     % Find index of adhesion in Adhesions (used to get XY coords of adhesion)
                             xDist = Adhesions.XYPoints(aidx,1) - Filaments.XYCoords{f}(midx,1); % X Distance between adhesion and attached filament monomer
                             yDist = Adhesions.XYPoints(aidx,2) - Filaments.XYCoords{f}(midx,2); % Y Distance between adhesion and attached filament monomer
                             SeparationDist = sqrt( xDist^2 + yDist^2 );                         % Distance between adhesion and attached filament monomer
                             StretchDist = SeparationDist - ModelParameters.AdhesionSpringEqLength;  % Stretch distance
                             StretchDist(StretchDist < 0) = 0;                                       % If stretch distance is less than Equilibrium length, StretchDist = 0;
                             ConnectionTension = StretchDist*ModelParameters.AdhesionSpringConstant; % Fa = ka*x (Calcualte spring force between adhesion and filament
                             
                             Tensions = [Tensions; ConnectionTension];              % concatenate Tensions, Filament Name, Adhesion Index, and Adhesion region information
                             AttachedFilament = [AttachedFilament; Filaments.Name(f)];
                             AdhesionIndex = [AdhesionIndex; aidx];
                             Regions = [Regions; Adhesions.RegionLocation(aidx,1)];
                             
                             ForceX = ConnectionTension*(xDist/SeparationDist);   % Fx = F*(x component of separation distance unit vector)
                             ForceY = ConnectionTension*(yDist/SeparationDist);   % Fy = F*(y component of separation distance unit vector)
                             
                             FAx = FAx + ForceX;                                  % Accumulate all the forces from adhesions being exerted on this filament structure
                             FAy = FAy + ForceY;
                             
                             % If Molecular Clutch is turned on, deactivation of adhesion rate is dependent on tension
                             kBreakRate =  base*exp(A*slope1*ConnectionTension) + C*exp(B*slope2*ConnectionTension);
                             kBreaks = [kBreaks; kBreakRate]; % Accumulate calculated deactivation rates (if needed for data recording)
                             
                             % IF rand < dt*kBreakRate AND Molecular Clutch is ON 
                             if rand(1) < (kBreakRate)*ModelParameters.TimeStep && ModelParameters.Adhesion_MolecularClutchOn
                                 BreakIdx = [BreakIdx; c]; % Acculumate indices of adhesions that need to be deactivated
                             end
                             
                         end
                 end
            end
            
            % Remove connections of deactivated Adhesion-Filament connections (only executes if MC is on. Otheriwse BreakIdx is always empty)-----------------
            if ~isempty(BreakIdx)
                 a = FAConnections.AdhesionIndex(BreakIdx,1);
                 Adhesions.AttachedFilamentName(a,1) = NaN;
                 FAConnections.AdhesionIndex(BreakIdx,:) = [];
                 FAConnections.FilamentName(BreakIdx,:) = [];
                 FAConnections.MonomerIndex(BreakIdx,:) = [];
                 %if ModelParameters.Adhesion_MolecularClutchOn
                 Adhesions.XYPoints(a,1) = NaN;
                 Adhesions.XYPoints(a,2) = NaN;
                 Adhesions.RegionLocation(a,:) = NaN;
                 Adhesions.ActiveStatus(a,1) = false;
                 Adhesions.AttachedFilamentName(a,1) = NaN;
                 %end
            end 
            %--------------------------------------------------------------
        end


end

