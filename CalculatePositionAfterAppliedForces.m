function [Filaments, Membrane, Adhesions, FAConnections, Tensions, kBreaks, Data] = ...
                        CalculatePositionAfterAppliedForces(Filaments,Membrane,Adhesions,FAConnections,ModelParameters)

   
    
    % Find all the main filaments of connected filament structures (group of attached filaments)
    idxMF = find(Filaments.Parent == 0);
    nMF = length(idxMF); % Total number of filament structures
    dt = ModelParameters.TimeStep;
    nS = size(Membrane.Segments,1); % Number of membrane segments
    %FilamentForce = ModelParameters.FilamentForce; %1; % Generic filament force along filament direction (pN)
    k  = ModelParameters.MembraneSpringConstant;
    D  = ModelParameters.SpringWidth;
    kT = ModelParameters.kT;
    FilamentTips = GetFilamentsTipLocations(Filaments);
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
    Data.FilamentTips.FilForce = [];
    Data.FilamentTips.FilLength = [];
    Data.FilamentTips.FilAngle = [];
    Data.FilamentTips.FilDist = [];
    Data.Forces = []; % Tx, Ty, Fm, Fb, Ax, Ay
    
    % Create membrane bins:
    MemBins =  (Membrane.Nodes(Membrane.Segments(1,1),1) - ModelParameters.SpringWidth/2:...
                ModelParameters.SegmentWidth+ModelParameters.SpringWidth:...
                Membrane.Nodes(Membrane.Segments(end,2),1) + ModelParameters.SpringWidth/2)';

    % FILAMENT SECTION
    if ~isempty(Filaments.Name)
        for MF = 1:nMF % For each filament structure
                L = 0; % Total length of Filament
                idx1 = find( Filaments.MainIndex == Filaments.Name(idxMF(MF)) ); % Find all filaments attached to same structure
                nF = length(idx1); % number of filaments in this structure
            % Calculate total length of filament structure
                for f = idx1'  
                    L = L + length(Filaments.MonomerIndices{f})*ModelParameters.MonomerLength; % L is total length of filament structure in nm
                end 
            % CALCULATE FORCES DUE TO FILAMENT-ADHESION CONNECTIONS and break connections where force is greater than threshold
               [FAx,FAy,Adhesions,FAConnections,Ts,kBs,Regions,AttachedFilament,AdhesionIndex] = CalculateForceDueToFAConnections(idx1,Filaments,Adhesions,FAConnections,ModelParameters);
                % Concatentate all tensions and data about adhesions attached to filaments in this structure.
                Tensions = [Tensions;Ts];
                kBreaks  = [kBreaks; kBs];
                Data.FAConnections.Tensions     = [Data.FAConnections.Tensions; single(Ts)];
                Data.FAConnections.Regions      = [Data.FAConnections.Regions; single(Regions)];
                Data.FAConnections.FilamentName = [Data.FAConnections.FilamentName; single(AttachedFilament)];
                Data.FAConnections.AdhesionIdx  = [Data.FAConnections.AdhesionIdx; single(AdhesionIndex)];

            % CALCULATE FORCES ACTING ON FILAMENT FROM MEMBRANE
                % Calculate sum of the normal forces from filements pushing on membrane (all from the same structure)
                [Fm,FilForces,FilLengths,FilAngles,FilDists] = TotalForceExertedByMembraneOnFilamentStructure(idx1,MemBins,Filaments,FilamentTips,Membrane,k,D,nS,ModelParameters);
                
            % CALCULATE GAMMA AND NEW POSITION FOR FILAMENT STRUCTURE
                Nu = ModelParameters.CytoplasmViscosity*(10^-6); % convert viscosity to pN/nm^2
                r = ModelParameters.FilamentDiameter; % nm;

                L(L<r) = r; % 
                gamma1 =  (4*pi*Nu*L) / ( 0.84+log(L/r) ); % gamma for random fluid motion (pN*s/nm) 2
                SD = sqrt(2*kT*gamma1/dt);                   
                Fx = SD*randn(1);
                Fy = SD*randn(1);
                
                tipXpositions = NaN(nF,1);
                for n = 1:nF % Loop through all the filaments attached to this structure
                    f = idx1(n);
                    yprevious = Filaments.XYCoords{f}(end,2);
                    
                    % START Calculate New Filament Position --------------------------------------------------------------------------
                    Filaments.XYCoords{f}(:,1) = Filaments.XYCoords{f}(:,1) + (Fx/gamma1)*dt + (FAx/gamma1)*dt;                  % Compute new position in X direction
                    Filaments.XYCoords{f}(:,2) = Filaments.XYCoords{f}(:,2) + (Fy/gamma1)*dt + (FAy/gamma1)*dt + (Fm/gamma1)*dt; % Compute new position in Y direction
                    % END Calculate New Filament Position ----------------------------------------------------------------------------
                   
                    % Record Filament Speed, Position, Region, etc of each tip .................................................
                    R = FindWhichRegionTheFilamentIsIn(Filaments.XYCoords{f}(end,1),Membrane,ModelParameters);
                    YDist = (Filaments.XYCoords{f}(end,2) - yprevious);
                    Data.FilamentTips.YSpeed       = [ Data.FilamentTips.YSpeed;       single(YDist/ModelParameters.TimeStep) ];
                    Data.FilamentTips.Region       = [ Data.FilamentTips.Region;       uint16(R) ];
                    Data.FilamentTips.FilamentName = [ Data.FilamentTips.FilamentName; single(Filaments.Name(f,1)) ];
                    Data.FilamentTips.XYPosition   = [ Data.FilamentTips.XYPosition;   single(Filaments.XYCoords{f}(end,:)) ];
                    Data.FilamentTips.StructName   = [ Data.FilamentTips.StructName;   single(Filaments.MainIndex(f,1)) ];
                    Data.FilamentTips.FilForce     = [ Data.FilamentTips.FilForce;     single(FilForces(n))];
                    Data.FilamentTips.FilLength    = [ Data.FilamentTips.FilLength;    single(FilLengths(n))];
                    Data.FilamentTips.FilAngle     = [ Data.FilamentTips.FilAngle;     single(FilAngles(n))];
                    Data.FilamentTips.FilDist      = [ Data.FilamentTips.FilDist;      single(FilDists(n))];
                    Data.Forces = [ Data.Forces; single([Fx, Fy, Fm, FAx, FAy]) ]; % Fx, Fy, Fm, Ax, Ay
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
                         FAConnections.FilamentName(BreakIdx,:)  = [];
                         FAConnections.MonomerIndex(BreakIdx,:)  = [];
                         for nn = 1:length(a)
                             aidx = find(Data.FAConnections.AdhesionIdx == a(nn));
                             Data.FAConnections.AdhesionIdx(aidx,:) = [];
                             Data.FAConnections.Tensions(aidx,:)    = [];
                             Data.FAConnections.Regions(aidx,:)     = [];
                             Data.FAConnections.FilamentName(aidx,:)= [];
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
                FF = CalculateForceFromFilamentsHittingMembraneSegment(Filaments,FilamentTips,MembraneOld,s,ModelParameters);
                if isnan(FF) || isinf(FF)
                    disp('NaN in FF')
                end
                Membrane.Nodes(Membrane.Segments(s,:),2) = MembraneOld.Nodes(MembraneOld.Segments(s,:),2) + (FF/gamma3)*dt; % Move y-position of segment nodes
            end
        end
        % MEMBRANE SECTION STEP 2 (move membrane based on spring tensions of neighboring segments)
        MembraneOld = Membrane;
        for s = MemStart:MemEnd
            FN = CalculateForceFromNeighboringMembraneSegments(MembraneOld,s,nS,k,D);
            Membrane.Nodes(Membrane.Segments(s,:),2) = MembraneOld.Nodes(MembraneOld.Segments(s,:),2) + (FN/gamma3)*dt; % Move y-position of segment nodes
            if any(isnan( Membrane.Nodes(Membrane.Segments(s,:),2) ) )
                disp('NaN in FN')
            end
        end
    %-------------------------------------------------------------------------------
    
end

%======================================================================================================
%======================================================================================================

function [Fm,FilForces,FilLengths,FilAngle,FilDist] = TotalForceExertedByMembraneOnFilamentStructure(FilInd,MemBins,Filaments,FilamentTips,Membrane,k,D,nS,ModelParameters)
            
            FilForces = [];  % Record all the individual filament forces
            FilLengths = []; % Record all the individual filament lengths
            FilAngle = [];   % Record all the individual filament orientation angles
            FilDist = [];

            Fm = 0; % Fm = Total force exerted on the filament structure.

            % Calculate the normal force each filament agains the membrane OR the force the membrane exerts on the filament structure.
            for f = FilInd'
%                     m = find( Membrane.Nodes(Membrane.Segments(:,1),1) - ModelParameters.SpringWidth/2  <=  FilamentTips(f,1) &... % Left segment endpoint  <= Filament X-end
%                               Membrane.Nodes(Membrane.Segments(:,2),1) + ModelParameters.SpringWidth/2  >=  FilamentTips(f,1) );   % Right segment endpoint >= Filament X-end

                    m = discretize(FilamentTips(f,1),MemBins); % Find interacting membrane segment for this filament
                    if ModelParameters.BrownianRatchetOn
                            if ~isnan(m)
                                y = Membrane.Nodes(Membrane.Segments(m,1),2) - FilamentTips(f,2); % Distance of filament tip from membrane
                                %|||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||
                                % Evaluate normal force of filament against membrane velocity eq.5 (Mogliner and Oser 1996);
                                theta = abs(atand( Filaments.UnitVector(f,1)/Filaments.UnitVector(f,2) ));
                                theta(theta < 1 ) = 1; % For model stability, keep theta greater than 1 degree
                                lambda = ModelParameters.PersistenceLength;
                                kT = ModelParameters.kT;
                                L  = length(Filaments.MonomerIndices{f})*ModelParameters.MonomerLength;
                                Lactual = L;
                                L(L< 30) =  30;
                                L(L>150) = 150;
                                delta = ModelParameters.MonomerLength*cosd(theta);
                                K = 4*lambda*kT/(L.^3*sind(theta).^2); % eq. B.1
                                if y >= 0
                                    Fn = -sqrt(2*kT*K/pi).*exp(-K*y.^2/(2*kT))./( erf(y*sqrt(K/(2*kT))) + 1 ); % Normal force calculation
                                else
                                    Fn = -sqrt(2*kT*K/pi).*exp(-K*y.^2/(2*kT))./( erfc( abs( y*sqrt(K/(2*kT)) ) ) );
                                end
                                %|||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||
                            else % If filment tip is outside either end of the membrane set F = 0
                                y = NaN;
                                theta = abs(atand( Filaments.UnitVector(f,1)/Filaments.UnitVector(f,2) ));
                                theta(theta < 1 ) = 1; % For model stability, keep theta greater than 1 degree
                                L  = length(Filaments.MonomerIndices{f})*ModelParameters.MonomerLength;
                                Lactual = L;
                                Fn = -10;
                            end
        
                            % for the cases where equation Fn fails
                            if isnan(Fn); Fn = -50; end
                            if isinf(Fn); Fn = -50; end
        
                            Fm = Fm + Fn; % Accumulate the forces applied to this filament structure
        
                            FilForces  = [FilForces;      -Fn]; 
                            FilLengths = [FilLengths; Lactual]; 
                            FilAngle   = [FilAngle;     theta];
                            FilDist    = [FilDist;      y];

                    else % Case where BR is turned off                           

                            if isnan(m) % If filament is outside membrane edge
                                y = NaN;
                                Fn = -1;
                            elseif FilamentTips(f,2) > (Membrane.Nodes(Membrane.Segments(m,1),2) - ModelParameters.MCmodelContactThreshold)
                                y = Membrane.Nodes(Membrane.Segments(m,1),2) - FilamentTips(f,2); % Distance of filament tip from membrane
                                Fn = -ModelParameters.MaxFilamentForceWithoutBR*Filaments.UnitVector(f,2);
                            else
                                y = Membrane.Nodes(Membrane.Segments(m,1),2) - FilamentTips(f,2); % Distance of filament tip from membrane
                                Fn = 0;
                            end
                            
                            theta = abs(atand( Filaments.UnitVector(f,1)/Filaments.UnitVector(f,2) ));
                            L  = length(Filaments.MonomerIndices{f})*ModelParameters.MonomerLength;
                            Lactual = L;

                            Fm = Fm + Fn; % Accumulate the forces applied to this filament structure

                            FilForces  = [FilForces;      -Fn]; 
                            FilLengths = [FilLengths; Lactual]; 
                            FilAngle   = [FilAngle;     theta];
                            FilDist    = [FilDist;          y];
                    end

            end

end

%======================================================================================================
%======================================================================================================

% function [FilInd,MemInd] = FindWhichFilamentsAreHittingMembraneAndWhere(idx,Filaments,Membrane,ModelParameters)
% 
%         FilInd = [];
%         MemInd = [];
%         SpringWidth = ModelParameters.SpringWidth;
%         nS = size(Membrane.Segments,1);
%          
%         for f = idx' % For each filament of a filament structure
%             Xend = Filaments.XYCoords{f}(end,1); % X-coordinate of filament tip
%             Yend = Filaments.XYCoords{f}(end,2); % Y-coordinate of filament tip
%             for m = 1:nS % For each membrane segment
%                    Yseg = Membrane.Nodes(Membrane.Segments(m,1),2); %  y-coordinate current Membrane segment
%                    % Test if filament tip is between membrane segment endpoints and is above membrane's y-position
%                    if Xend >= Membrane.Nodes(Membrane.Segments(m,1),1)-SpringWidth/2 &&... % Filament X-end >= Left segment endpoint
%                       Xend <= Membrane.Nodes(Membrane.Segments(m,2),1)+SpringWidth/2   % Filament X-end <= Right segment endpoint
%                           %if Yend > (Yseg-ModelParameters.ContactThresholdForMembrane) % Filament tip > Membrane segment     
%                                 % Record this filament and the membrane it is hitting
%                                 FilInd = [FilInd; f];
%                                 MemInd = [MemInd; m];
%                           %end
%                    end   
%                 end
%         end
% 
% end

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

function Ff = CalculateForceFromFilamentsHittingMembraneSegment(Filaments,FilamentTips,Membrane,m,ModelParameters)
        
        %global FF
        NodeLeft    = Membrane.Nodes(Membrane.Segments(m,1),:);
        NodeRight   = Membrane.Nodes(Membrane.Segments(m,2),:);
        SegmentYPos = Membrane.Nodes(Membrane.Segments(m,1),2);
        
        idx = find( FilamentTips(:,1) >= NodeLeft(:,1)  - ModelParameters.SpringWidth/2 & ...
                    FilamentTips(:,1) <= NodeRight(:,1) + ModelParameters.SpringWidth/2);
                
        % Add up all the y-components of all the filament forces on this membrane segment  
        Ff = 0;
        for f = idx'
             if ModelParameters.BrownianRatchetOn
                    y = SegmentYPos - FilamentTips(f,2); % Distance of filament tip from membrane
                    %|||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||
                    % Evaluate force of filament against membrane velocity eq.5 (Mogliner and Oser 1996);
                    theta = abs(atand( Filaments.UnitVector(f,1)/Filaments.UnitVector(f,2) ));
                    theta(theta < 1 ) = 1; % For model stability, keep theta greater than 1 degree 
                    lambda = ModelParameters.PersistenceLength;
                    kT = ModelParameters.kT;
                    L = length(Filaments.MonomerIndices{f})*ModelParameters.MonomerLength;
                    L(L<30) = 30;
                    delta = ModelParameters.MonomerLength*cosd(theta);
                    K = 4*lambda*kT/(L.^3*sind(theta).^2); % eq. B.1
                    %Fm = (kT*exp(-K*y.^2/(2*kT)))./( sqrt(pi*kT/(2*K))*( erf(y*sqrt(K/(2*kT))) + 1 ) );
                    if y >= 0
                        Fm = sqrt(2*kT*K/pi).*exp(-K*y.^2/(2*kT))./( erf(y*sqrt(K/(2*kT))) + 1 ); % Normal force calculation
                    else
                        Fm = sqrt(2*kT*K/pi).*exp(-K*y.^2/(2*kT))./( erfc( abs( y*sqrt(K/(2*kT)) ) ) );
                    end
        
                    Fm( isnan(Fm) ) = 0;
                    Fm( isinf(Fm) ) = 50;
                    %|||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||
                    Ff = Ff + Fm;

            else % Case where BR is turned off
                    if ~isnan(m) && FilamentTips(f,2) >= (Membrane.Nodes(Membrane.Segments(m,1),2) - ModelParameters.MCmodelContactThreshold)
                        Fm = ModelParameters.MaxFilamentForceWithoutBR*Filaments.UnitVector(f,2);
                    else
                        Fm = 0;
                    end
                    Fm( isnan(Fm) ) = 0;
                    Fm( isinf(Fm) ) = 0;
                    Ff = Ff + Fm; 
            end
        end   
       
        %Ff(abs(Ff) > 1) = 1;
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
        FAx = [];
        FAy = [];
        BreakIdx = [];
        BreakIdx2 = [];
        Regions = [];
        index = 0;
        if ~isempty(FAConnections.AdhesionIndex)
            for f = idx'
                 con = find(FAConnections.FilamentName == Filaments.Name(f,1)); % Find adhesions attached to this filament
                 if ~isempty(con)
                         for c = con'
                             index = index + 1;
                             midx  = find( Filaments.MonomerIndices{f} == FAConnections.MonomerIndex(c) ); % Find index of attached monomer (used to get XY coords of monomer)
                             aidx  = FAConnections.AdhesionIndex(c,1);                                     % Find index of adhesion in Adhesions (used to get XY coords of adhesion)
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
                             
                             FAx = [FAx;ForceX];                                  % Accumulate all the forces from adhesions being exerted on this filament structure
                             FAy = [FAy;ForceY];
                             
                             % If Molecular Clutch is turned on, deactivation of adhesion rate is dependent on tension
                             kBreakRate =  base*exp(A*slope1*ConnectionTension) + C*exp(B*slope2*ConnectionTension);
                             kBreaks    = [kBreaks; kBreakRate]; % Accumulate calculated deactivation rates (if needed for data recording)
                            
                             % IF rand < dt*kBreakRate AND Molecular Clutch is ON 
                             if rand(1) < (kBreakRate)*ModelParameters.TimeStep && ModelParameters.Adhesion_MolecularClutchOn
                                 BreakIdx  = [BreakIdx; c]; % Acculumate indices of adhesions that need to be deactivated
                                 BreakIdx2 = [BreakIdx2; index];
                             end
                             
                         end
                 end
            end
            
            % Remove connections of deactivated Adhesion-Filament connections (only executes if MC is on. Otheriwse BreakIdx is always empty)-----------------
            if ~isempty(BreakIdx)
                 a = FAConnections.AdhesionIndex(BreakIdx,1);
                 FAConnections.AdhesionIndex(BreakIdx,:) = [];
                 FAConnections.FilamentName(BreakIdx,:)  = [];
                 FAConnections.MonomerIndex(BreakIdx,:)  = [];
                 
                 Adhesions.XYPoints(a,1) = NaN;
                 Adhesions.XYPoints(a,2) = NaN;
                 Adhesions.RegionLocation(a,1) = NaN;
                 Adhesions.ActiveStatus(a,1) = false;
                 Adhesions.AttachedFilamentName(a,1) = NaN;
                 
                 Tensions(BreakIdx2) = [];
                 kBreaks(BreakIdx2) = [];
                 AttachedFilament(BreakIdx2) = [];
                 AdhesionIndex(BreakIdx2)    = [];
                 FAx(BreakIdx2) = [];
                 FAy(BreakIdx2) = [];
                 Regions(BreakIdx2) = [];
            end 
            %--------------------------------------------------------------
        end
        FAx = sum(FAx);  if isempty(FAx); FAx = 0; end                                % Accumulate all the forces from adhesions being exerted on this filament structure
        FAy = sum(FAy);  if isempty(FAy); FAy = 0; end             

end

