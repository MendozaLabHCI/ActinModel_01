function [FAConnections, Adhesions] = FilamentAndAdhesionConnect(FAConnections,Filaments,Adhesions,ModelParameters)

    nF = size(Filaments.XYCoords,1);
    
    MinDist = ModelParameters.MonomerLength; % Minimum distance for connection. 
    % Test if Adhesions are within a distance D of a filament/Monomer
    if ~isempty(Filaments.XYCoords)
        for f = 1:nF
           % Grab start and end points of filament
           FS = Filaments.XYCoords{f}(1,:); % Filament starting point
           FE = Filaments.XYCoords{f}(end,:); % Filament end point
           FilamentLength = sqrt( (FE(1)-FS(1))^2 + (FE(2)-FS(2))^2 );
           % Calculate distance from Filament starting point to all adhesions
                distFStoAD = sqrt( (FS(1)-Adhesions.XYPoints(:,1)).^2 + (FS(2)-Adhesions.XYPoints(:,2)).^2 );
           % Calculate distance from Filament end point to all adhesions
                distFEtoAD = sqrt( (FE(1)-Adhesions.XYPoints(:,1)).^2 + (FE(2)-Adhesions.XYPoints(:,2)).^2 );
           % Calculate Triangle height using Heron's formula using FilamentLength as the base
                s = (FilamentLength + distFStoAD + distFEtoAD)./2;
                h = sqrt( 4.*s.*(s-distFStoAD).*(s-distFEtoAD).*(s-FilamentLength)./(FilamentLength^2) );
           % Now find Adhesions within a distance D of Filament line
                idx = find( h <= MinDist | distFStoAD <= MinDist | distFEtoAD <= MinDist );
           % If there are any adhesions within a distance D of the filament line, 
           % find out if these adhesions are within a distance D of any of the monomers. (this may be overdoing it).
                for a = idx'
                   D = sqrt((Filaments.XYCoords{f}(:,1)-Adhesions.XYPoints(a,1)).^2 + ...
                            (Filaments.XYCoords{f}(:,2)-Adhesions.XYPoints(a,2)).^2);
                   [Dmin,idx2] = min(D); 
                   if Dmin <= MinDist
                        % Check if this adhesion and monomer have already been connected -----------
                        check1 = isempty(find( FAConnections.FilamentName  == Filaments.Name(f,1) & ...
                                               FAConnections.MonomerIndex  == Filaments.MonomerIndices{f}(idx2,1) ,1));
                        % ------------------------------------------------
                        if check1 && isnan(Adhesions.AttachedFilamentName(a,1)) 
                            FAConnections.AdhesionIndex = [FAConnections.AdhesionIndex; a];
                            FAConnections.FilamentName  = [FAConnections.FilamentName; Filaments.Name(f,1)];
                            FAConnections.MonomerIndex  = [FAConnections.MonomerIndex; Filaments.MonomerIndices{f}(idx2,1)];
                            Adhesions.AttachedFilamentName(a,1) = Filaments.Name(f,1);
                        end
                   end
                end
                
        end
    end
    
end