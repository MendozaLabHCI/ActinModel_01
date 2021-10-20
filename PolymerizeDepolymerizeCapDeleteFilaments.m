function [Filaments,Adhesions,FAConnections,PolymCoeffData] = PolymerizeDepolymerizeCapDeleteFilaments(Filaments,Adhesions,Membrane,FAConnections,ModelParameters)
    %
    % Filaments = PolymerizeDepolymerizeCapDeleteFilaments(Filaments,ModelParamters)
    %
    % KRC 12/07/2020
   
        nF = length(Filaments.XYCoords); % number of filaments
        D = ModelParameters.MonomerLength; %(nm)
        FilamentTips = GetFilamentsTipLocations(Filaments);
        PolymCoeffData.PolymCoeff = [];
        PolymCoeffData.FilamentName = [];
        
        % FOR BARBED END OF FILAMENT
        for f = 1:nF
            % Cap selected filaments first ----------------------------------------
            cap_test = rand(1) < ModelParameters.k_cap*ModelParameters.TimeStep;
            if cap_test
                Filaments.IsCapped(f,1) = true;
            end

            % Check if it's capped first before polymerizing or depolymerizing
            if ~Filaments.IsCapped(f,1) 
                    [~,~,PolymCoeff] = TestIfFilamentIsHittingMembraneAndWhere(f,Filaments,FilamentTips,Membrane,ModelParameters);
                    % Run random selection test ----------------------------------------------------------------------------
                    barbed_on_test  =  rand(1) < PolymCoeff*ModelParameters.k_on_barbed*ModelParameters.TimeStep;
                    barbed_off_test =  rand(1) < PolymCoeff*ModelParameters.k_off_barbed*ModelParameters.TimeStep;

                    PolymCoeffData.PolymCoeff   = [PolymCoeffData.PolymCoeff; PolymCoeff];
                    PolymCoeffData.FilamentName = [PolymCoeffData.FilamentName; Filaments.Name(f,1)];

                    % Polymerize barbed end for selected filaments (as long as depolymerized doesn't cancel it out) ------------------------------------   
                    if barbed_on_test && ~barbed_off_test   
                        % Add next indices and next XY point on barbed end of filament
                        Filaments.MonomerIndices{f} = [Filaments.MonomerIndices{f}; Filaments.MonomerIndices{f}(end,1)+1]; 
                        Filaments.XYCoords{f} = [Filaments.XYCoords{f}; [Filaments.XYCoords{f}(end,1) + D*Filaments.UnitVector(f,1),...
                                                                         Filaments.XYCoords{f}(end,2) + D*Filaments.UnitVector(f,2)] ];
                    % De-Polymerize barbed end for selected filaments (as long as polymerized doesn't cancel it out) ----------------------------------------
                    elseif barbed_off_test && ~barbed_on_test 
                        L = length(Filaments.MonomerIndices{f}); % Get current length of filament
                        if L > 1 && isempty(find(Filaments.Parent == Filaments.Name(f),1)) % If length > 1, remove the end monomer. 
                            midx = Filaments.MonomerIndices{f}(end,1); % Grab index of monomer about to be deleted
                            Filaments.MonomerIndices{f} = Filaments.MonomerIndices{f}(1:end-1,1);
                            Filaments.XYCoords{f} = Filaments.XYCoords{f}(1:end-1,:);
                            % Check if Adhesion is attached --------------------------------------------------------
                            idx = find( FAConnections.FilamentName == Filaments.Name(f,1) & FAConnections.MonomerIndex == midx);
                            % If so, reset adhesion, and remove connection  
                            if ~isempty(idx) 
                                Adhesions.AttachedFilamentName(FAConnections.AdhesionIndex(idx,1),1) = NaN;
                                FAConnections.AdhesionIndex(idx,:) = [];
                                FAConnections.FilamentName(idx,:)  = [];
                                FAConnections.MonomerIndex(idx,:)  = [];
                            end
                            %---------------------------------------------------------------------------------------
                        end
                    end
            else
                PolymCoeffData.PolymCoeff   = [PolymCoeffData.PolymCoeff; 0];
                PolymCoeffData.FilamentName = [PolymCoeffData.FilamentName; Filaments.Name(f,1)];
            end
        end 
      
        % FOR POINTED END OF FILAMENT
        idxD = [];
        for f = 1:nF
            % De-Polymerize pointed end for selected filaments ----------------------------------------
            pointed_off_test = rand(1) < ModelParameters.k_off_pointed*ModelParameters.TimeStep;
            L = length(Filaments.MonomerIndices{f}); % Get current length of filament
            if pointed_off_test && Filaments.Parent(f,1) == 0 
                % Record Name and Current Monomer index in case this filament is deleted 
                % and we can still find its daughers (BranchIdx) afterwards
                Name = Filaments.Name(f,1);
                MonomerIndice = Filaments.MonomerIndices{f}(1,1);
                
                % If length is 1 unit record filament as needing to be deleted
                if L == 1  
                    idxD = [idxD;f];
                end
                
                % Remove monomer
                midx = Filaments.MonomerIndices{f}(1,1); % First, grab index of monomer about to be deleted
                Filaments.MonomerIndices{f} = Filaments.MonomerIndices{f}(2:end,1);
                Filaments.XYCoords{f} = Filaments.XYCoords{f}(2:end,:);
                % Check if adhesion is attached --------------------------------------------------------
                idx = find( FAConnections.FilamentName == Filaments.Name(f,1) & ...
                            FAConnections.MonomerIndex == midx);
                % If so, reset adhesion, and remove connection  
                if ~isempty(idx) 
                    Adhesions.AttachedFilamentName(FAConnections.AdhesionIndex(idx,1),1) = NaN;
                    FAConnections.AdhesionIndex(idx,:) = [];
                    FAConnections.FilamentName(idx,:)  = [];
                    FAConnections.MonomerIndex(idx,:)  = [];
                end
                %---------------------------------------------------------------------------------------
                % Check if there was a branch point there. 
                % If so make the branch a main filament (i.e. set Filaments.Parent = 0;
                BranchIdx = find(Filaments.Parent==Name);
                idx1 = [];
                for k = 1:length(BranchIdx)
                    if Filaments.ParentIndex(BranchIdx(k),1) == MonomerIndice
                        idx1 = BranchIdx(k); break;
                    end
                end
                
                % Set branch as Main filament, and find all its daughter
                % filaments and rename their MainIndex (Main Filament) to
                % this branch.
                if ~isempty(idx1)
                   % Setting first branch as Main Filament
                   Filaments.Parent(idx1,1) = 0; 
                   Filaments.ParentIndex(idx1,1) = 0;
                   Filaments.MainIndex(idx1,1) = Filaments.Name(idx1,1);
                   %------------------------------------------------------------------------------
                   % Find all sub-daughters of the new Main Filament and
                   % change their MainIndex to the new Filament name
                   idx2 = find(Filaments.Parent == Filaments.Name(idx1,1));
                   if ~isempty(idx2)
                       idx3 = idx2;
                       while true
                            idx4 = find( ismember(Filaments.Parent,Filaments.Name(idx2,1)) );
                            if isempty(idx4); break; end
                            idx3 = [idx3; idx4];
                            idx2 = idx4;
                       end
                       Filaments.MainIndex(idx3,1) = Filaments.Name(idx1,1);
                   end
                   %------------------------------------------------------------------------------
                end
                
                
            end
        end
        % Remove deleted filaments before continuing on
        Filaments.Name(idxD) = [];
        Filaments.MonomerIndices(idxD) = [];
        Filaments.XYCoords(idxD) = [];
        Filaments.UnitVector(idxD,:) = [];
        Filaments.IsCapped(idxD) = [];
        Filaments.MainIndex(idxD) = [];
        Filaments.Parent(idxD) = [];
        Filaments.ParentIndex(idxD) = [];
        
end

%======================================================================================================
%======================================================================================================

% function FilamentTips = GetFilamentsTipLocations(Filaments)
%     
%     nF = size(Filaments.XYCoords,1);
%     nF(isempty(Filaments.XYCoords)) = 0;
%     FilamentTips = zeros(nF,2);
%     for f = 1:nF
%         FilamentTips(f,:) = Filaments.XYCoords{f}(end,:);
%     end
%     
% end

%======================================================================================================
%======================================================================================================

function [MemInd,DistFromMembrane,PolymCoeff] = TestIfFilamentIsHittingMembraneAndWhere(f,Filaments,FilamentTips,Membrane,ModelParameters)
       
        DistFromMembrane = [];
        PolymCoeff = 1;
        k  = ModelParameters.MembraneSpringConstant;
        D  = ModelParameters.SpringWidth;
        nS = size(Membrane.Segments,1); 
        SpringWidth = ModelParameters.SpringWidth;
        % Find which membrane segment's x-coordinates are surrounding the endponit of the filament
        MemInd = find( Membrane.Nodes(Membrane.Segments(:,1),1)-SpringWidth/2 <= FilamentTips(f,1) &... % Left segment endpoint  <= Filament X-end
                       Membrane.Nodes(Membrane.Segments(:,2),1)+SpringWidth/2 >= FilamentTips(f,1) );   % Right segment endpoint >= Filament X-end 
                   
        if ~isempty(MemInd)
            DistFromMembrane = Membrane.Nodes(Membrane.Segments(MemInd,1),2) - FilamentTips(f,2);
            if DistFromMembrane > ModelParameters.ContactThresholdForMembrane
                PolymCoeff = 1;
            else
%               mSlope = 1/ModelParameters.ContactThresholdForMembrane;
                % Initalize spring force on left and right side of segment as zero
                Fy1 = 0;
                Fy2 = 0;
                % Calcualte spring force from spring on left side of membrane segment
                if ~isequal(MemInd,1) % Don't calcualte if this is the far left membrane segment
                    Fdir1 = sign(Membrane.Nodes(Membrane.Segments(MemInd-1,2),2) - Membrane.Nodes(Membrane.Segments(MemInd,1),2));
                    delX1 = abs(Membrane.Nodes(Membrane.Segments(MemInd-1,2),1) - Membrane.Nodes(Membrane.Segments(MemInd,1),1));
                    delY1 = abs(Membrane.Nodes(Membrane.Segments(MemInd-1,2),2) - Membrane.Nodes(Membrane.Segments(MemInd,1),2));
                    Fs1 = k*(sqrt(delX1^2 + delY1^2)-D); % Total spring force
                    Fy1 = Fdir1*Fs1*sin(atan(delY1/delX1));  % Y component of spring force
                end
                % Calcualte spring force from spring on right side of membrane segment
                if ~isequal(MemInd,nS) % Don't calcualte if this is the far left membrane segment
                    Fdir2 = sign(Membrane.Nodes(Membrane.Segments(MemInd+1,1),2) - Membrane.Nodes(Membrane.Segments(MemInd,2),2));
                    delX2 = abs(Membrane.Nodes(Membrane.Segments(MemInd,2),1) - Membrane.Nodes(Membrane.Segments(MemInd+1,1),1));
                    delY2 = abs(Membrane.Nodes(Membrane.Segments(MemInd,2),2) - Membrane.Nodes(Membrane.Segments(MemInd+1,1),2));
                    Fs2 = k*(sqrt(delX2^2+delY2^2)-D); % Total spring force
                    Fy2 = Fdir2*Fs2*sin(atan(delY2/delX2));  % Y component of spring force
                end
                % Calculate Boundary Force of Filament crossing membrane
                CrossDelta = Filaments.XYCoords{f}(end,2) - Membrane.Nodes(Membrane.Segments(MemInd,1),2); % How much filament crosses membrane segment
                CrossDelta(CrossDelta < 0) = 0;
                Fb = -ModelParameters.BoundaryForceSpringConstant*CrossDelta;
        
                %display(Fb)
                %-------------------------------------------------------------
                if (Fy1+Fy2) < 0 % Only add the contribution if the total force of the two springs points 'down' (-y direction) towards the filament
                    Fm = abs( Fy1 + Fy2 + Fb );% + ModelParameters.Fmin;
                else
                    Fm = abs( Fb );
                end
                % Create a polymerization coefficient that varies from 0-1 in the 15nm window below membrane.
                % It varies with how much force the membrane is pushing on it.
                PolymCoeff = exp(-1*Fm/4.11); 
                
                % Create a polymerization coefficient that varies from 0-1 in the contact threshold region of membrane
                % and below that it is 1 and above membrane it is 0;
%                 PolymCoeff(DistFromMembrane < 0) = 0; % For the case the filament is above the membrane
%                 PolymCoeff(DistFromMembrane > ModelParameters.ContactThresholdForMembrane) = 1;
            end
        end

end

