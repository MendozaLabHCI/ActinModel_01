function [Adhesions,FAConnections] = ManageAdhesions(Filaments,Adhesions,FAConnections,Membrane,ModelParameters)

    SpringWidth = ModelParameters.SpringWidth;
    nR = size(Membrane.Segments,1);
    nA = size(Adhesions.ActiveStatus,1);     
    Adhesions.RegionNodes = zeros(4,2,nR);   
    
            % Update adhesion regions and decide if adhesion points should
            % stay or be deleted
            for r = 1:nR % 
               XY_LT = [ Membrane.Nodes(Membrane.Segments(r,1),1) - SpringWidth/2,... %  Left-Top xy-coordinates of current adhesion region
                         Membrane.Nodes(Membrane.Segments(r,1),2) ];        
               XY_RT = [ Membrane.Nodes(Membrane.Segments(r,2),1) + SpringWidth/2,...  %  Right-Top xy-coordinates of current adhesion region
                         Membrane.Nodes(Membrane.Segments(r,2),2) ]; 
               XY_RB = [ Membrane.Nodes(Membrane.Segments(r,2),1) + SpringWidth/2,...  %  Right-Bottom xy-coordinates of current adhesion region
                         Membrane.Nodes(Membrane.Segments(r,2),2) - ModelParameters.AdhesionRegionDepth ];
               XY_LB = [ Membrane.Nodes(Membrane.Segments(r,1),1) - SpringWidth/2,... %  Left-Bottom xy-coordinates of current adhesion region
                         Membrane.Nodes(Membrane.Segments(r,1),2) - ModelParameters.AdhesionRegionDepth ];   
               
               % Update Current Regions
               Adhesions.RegionNodes(1,:,r) = XY_LT;           
               Adhesions.RegionNodes(2,:,r) = XY_RT;            
               Adhesions.RegionNodes(3,:,r) = XY_RB; 
               Adhesions.RegionNodes(4,:,r) = XY_LB; 
            end  
            
            
            % START Mirror Points below lower boundary ---------------------------------------------------------------------------------
                idx3 = find(~isnan(Adhesions.RegionLocation));
                R = unique(Adhesions.RegionLocation(idx3))';
                for r = R % Find regions with adhesions in them
                    % For the adhesions in this regionn, find which ones are active and below the bottom edge of the region
                    idx4 = find( (Adhesions.ActiveStatus == true) & ...
                                 (Adhesions.XYPoints(:,2) < Adhesions.RegionNodes(4,2,r)) & ...
                                 (Adhesions.RegionLocation == r));
                    if ~isempty(idx4)  
                        for idx5 = idx4'
                            Adhesions.XYPoints(idx5,2) = Adhesions.XYPoints(idx5,2) + ModelParameters.AdhesionRegionDepth; 
                            % Find out if adhesion is connected to a filament
                            idx6 = find(FAConnections.AdhesionIndex == idx5);
                                % If so, reset adhesion, and remove connection  
                                if ~isempty(idx6) 
                                    Adhesions.AttachedFilamentName(idx5,1) = NaN;
                                    FAConnections.AdhesionIndex(idx6,:) = [];
                                    FAConnections.FilamentName(idx6,:)  = [];
                                    FAConnections.MonomerIndex(idx6,:)  = [];
                                end
                        end
                    end
                end
            % END Mirror Points below lower boundary -----------------------------------------------------------------------------------
           
            
            % ACTIVATE Adhesions randomly selected for activation ----------
            idx1 = find(~Adhesions.ActiveStatus);
            
            for a = idx1' 
                r = randi(nR,1);
                if rand(1) < ModelParameters.Adhesion_ActivationRate * ModelParameters.TimeStep
                    delX = Adhesions.RegionNodes(2,1,r) - Adhesions.RegionNodes(1,1,r);
                    delY = Adhesions.RegionNodes(1,2,r) - Adhesions.RegionNodes(4,2,r);
                    
                    Adhesions.XYPoints(a,:) = [delX*rand(1) + Adhesions.RegionNodes(1,1,r), ...
                                               delY*rand(1) + Adhesions.RegionNodes(4,2,r)];
                    Adhesions.RegionLocation(a,1) = r;
                    Adhesions.ActiveStatus(a,1) = true;
                    Adhesions.AttachedFilamentName(a,1) = NaN;
                end
            end
            
            % DE-ACTIVATE Adhesions randomly elected for de-activation ----------
            if ModelParameters.Adhesion_MolecularClutchOn
                ModelParameters.Adhesion_DeActivationRate = 0.0;
            end
            
            idx2 = find(Adhesions.ActiveStatus); % Find all the active adhesions
            for a = idx2' % For each adhesion, test if it is deactivated
                if rand(1) < ModelParameters.Adhesion_DeActivationRate * ModelParameters.TimeStep
                    % Deactivate selected adhesion by resetting its parameters
                    Adhesions.XYPoints(a,:) = [NaN,NaN];
                    Adhesions.RegionLocation(a,:) = NaN;
                    Adhesions.ActiveStatus(a,1) = false;
                    Adhesions.AttachedFilamentName(a,1) = NaN;
                    % Check if the adhesion is attached to a filament monomer ------------------------------
                    idx6 = find(FAConnections.AdhesionIndex == a);
                    % If so, reset adhesion, and remove connection  
                    if ~isempty(idx6) 
                        FAConnections.AdhesionIndex(idx6,:) = [];
                        FAConnections.FilamentName(idx6,:)  = [];
                        FAConnections.MonomerIndex(idx6,:)  = [];
                    end
                    %---------------------------------------------------------------------------------------
                end
            end

end