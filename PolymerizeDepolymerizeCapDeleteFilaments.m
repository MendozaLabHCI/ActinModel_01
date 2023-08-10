function [Filaments,Adhesions,FAConnections,PolymCoeffData] = PolymerizeDepolymerizeCapDeleteFilaments(Filaments,Adhesions,Membrane,FAConnections,ModelParameters)
    %
    % Filaments = PolymerizeDepolymerizeCapDeleteFilaments(Filaments,ModelParamters)
    %
    % KRC 12/07/2020
   
        nF = length(Filaments.XYCoords); % number of filaments
        D = ModelParameters.MonomerLength; %(nm)
        FilamentTips = GetFilamentsTipLocations(Filaments);
        FilamentEnds = GetFilamentsEndLocations(Filaments);
        PolymCoeffData.PolymCoeff = [];
        PolymCoeffData.FilamentName = [];
        
        %% FOR BARBED END OF FILAMENT
        for f = 1:nF
                    % Cap selected filaments first ----------------------------------------
                    cap_test = rand(1) < ModelParameters.k_cap*ModelParameters.TimeStep;
                    if cap_test
                        Filaments.IsCapped(f,1) = true;
                    end

                    % Check if it's capped before polymerizing barbed end
                    if ~Filaments.IsCapped(f,1) 
                            PolymCoeff = CalculatePolymerizationCoefficient(f,Filaments,FilamentTips,Membrane,ModelParameters);
                            
                            % Run random selection test ----------------------------------------------------------------------------
                            barbed_on_test  =  rand(1) < PolymCoeff*ModelParameters.k_on_barbed*ModelParameters.TimeStep;
                            barbed_off_test =  rand(1) < PolymCoeff*ModelParameters.k_off_barbed*ModelParameters.TimeStep;

                            PolymCoeffData.PolymCoeff   = [PolymCoeffData.PolymCoeff; single(PolymCoeff)];
                            PolymCoeffData.FilamentName = [PolymCoeffData.FilamentName; single(Filaments.Name(f,1))];

                            % Polymerize barbed end for selected filaments (as long as depolymerized doesn't cancel it out) ------------------------------------   
                            if barbed_on_test && ~barbed_off_test   
                                % Add next indices and next XY point on barbed end of filament
                                Filaments.MonomerIndices{f} = [Filaments.MonomerIndices{f}; Filaments.MonomerIndices{f}(end,1)+1]; 
                                Filaments.XYCoords{f} = [Filaments.XYCoords{f}; [Filaments.XYCoords{f}(end,1) + D*Filaments.UnitVector(f,1),...
                                                                                 Filaments.XYCoords{f}(end,2) + D*Filaments.UnitVector(f,2)] ];
                            % De-Polymerize barbed end for selected filaments (as long as polymerized doesn't cancel it out) ----------------------------------------
                            %elseif barbed_off_test && ~barbed_on_test 
%                                 L = length(Filaments.MonomerIndices{f}); % Get current length of filament
%                                 if L > 1 && isempty(find(Filaments.Parent == Filaments.Name(f),1)) % If length > 1, remove the end monomer. 
%                                     midx = Filaments.MonomerIndices{f}(end,1); % Grab index of monomer about to be deleted
%                                     Filaments.MonomerIndices{f} = Filaments.MonomerIndices{f}(1:end-1,1);
%                                     Filaments.XYCoords{f} = Filaments.XYCoords{f}(1:end-1,:);
%                                     % Check if Adhesion is attached --------------------------------------------------------
%                                     idx = find( FAConnections.FilamentName == Filaments.Name(f,1) & FAConnections.MonomerIndex == midx);
%                                     % If so, reset adhesion, and remove connection  
%                                     if ~isempty(idx) 
%                                         Adhesions.AttachedFilamentName( FAConnections.AdhesionIndex(idx,1),1 ) = NaN;
%                                         FAConnections.FilamentName(idx,:)  = [];
%                                     end
%                                     %---------------------------------------------------------------------------------------
%                                 end
                            end
                    else
                        PolymCoeffData.PolymCoeff   = [PolymCoeffData.PolymCoeff; -1]; % Set PC for capped filaments to -1
                        PolymCoeffData.FilamentName = [PolymCoeffData.FilamentName; Filaments.Name(f,1)];
                    end
        end 
      
        %% FOR POINTED END OF FILAMENT
        idxD = [];
        for f = 1:nF
                    % De-Polymerize pointed end for selected filaments ----------------------------------------
                    pointed_off_test = rand(1) < ModelParameters.k_off_pointed*ModelParameters.TimeStep;
                    L = length(Filaments.MonomerIndices{f}); % Get current length of filament

                    if pointed_off_test && Filaments.Parent(f,1) == 0 % Parent = 0 means filament is a main filament
                        % Record Name and Current Monomer index in case this filament is deleted 
                        % and we can still find its daughers (BranchIdx) afterwards
                        Name = Filaments.Name(f,1);
                        MonomerIndice = Filaments.MonomerIndices{f}(1,1);

                        % If parent filement length is 2 monomers long remove it 
                        if L <= 2  
                            dName = Filaments.Name(f,1);
                            idxD = [idxD; dName];
                            % Check if there are branches (and branches on branches...)
                            ind = find(  Filaments.MainIndex == dName );
                            if ~isempty(ind)
                                BranchName = Filaments.Name(ind,1);
                                idxD = [idxD;BranchName];
                            end
                        end

                        % Remove monomer
                        midx = Filaments.MonomerIndices{f}(1,1); % First, grab index of monomer about to be deleted
                        Filaments.MonomerIndices{f} = Filaments.MonomerIndices{f}(2:end,1);
                        Filaments.XYCoords{f} = Filaments.XYCoords{f}(2:end,:);
                        % Check if adhesion is attached --------------------------------------------------------
                        idx = find( FAConnections.FilamentName == Filaments.Name(f,1) & ...
                                    FAConnections.MonomerIndex == midx);
                        % If so, keep adhesion active, but remove connection  
                        if ~isempty(idx) 
                            Adhesions.AttachedFilamentName( FAConnections.AdhesionIndex(idx,1),1 ) = NaN;
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
        
        
        %|||||||||||||||||||||||||||||||||||||||||||||||||||||||||||
        % Remove deleted filaments before continuing on
        idxD = unique(idxD);
        for n = 1:length(idxD)
            name = idxD(n);
                idx = find( Filaments.Name == name );
                % Delete filament
                Filaments.Name(idx) = [];
                Filaments.MonomerIndices(idx) = [];
                Filaments.XYCoords(idx) = [];
                Filaments.UnitVector(idx,:) = [];
                Filaments.IsCapped(idx) = [];
                Filaments.MainIndex(idx) = [];
                Filaments.Parent(idx) = [];
                Filaments.ParentIndex(idx) = [];
                
                % Remove any Filament-Adhesion connections
                idx2 = find( FAConnections.FilamentName == name );
                if ~isempty(idx2)
                    FAConnections.AdhesionIndex(idx2) = [];   
                    FAConnections.FilamentName(idx2) = []; 
                    FAConnections.MonomerIndex(idx2) = []; 
                end
                
                % Update attached Adhesions connections   
                idx3 = find( Adhesions.AttachedFilamentName == name );
                if ~isempty(idx3)
                    Adhesions.AttachedFilamentName(idx3) = NaN;
                end
        end
        
        
        
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

function PolymCoeff = CalculatePolymerizationCoefficient(f,Filaments,FilamentTips,Membrane,ModelParameters)
       
        DistFromMembrane = [];
        PolymCoeff = 0;
        % k  = ModelParameters.MembraneSpringConstant;
        % D  = ModelParameters.SpringWidth;
        nS = size(Membrane.Segments,1); 
        
        % Find which membrane segment's x-coordinates are surrounding the endpoint of the filament
        m = find( Membrane.Nodes(Membrane.Segments(:,1),1) - ModelParameters.SpringWidth/2  <=  FilamentTips(f,1) &... % Left segment endpoint  <= Filament X-end
                  Membrane.Nodes(Membrane.Segments(:,2),1) + ModelParameters.SpringWidth/2  >=  FilamentTips(f,1) );   % Right segment endpoint >= Filament X-end 
                   
        if ~isempty(m)
                if  ModelParameters.BrownianRatchetOn
                        y = Membrane.Nodes(Membrane.Segments(m,1),2) - FilamentTips(f,2); % Distance of filament tip from membrane
                        %|||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||
                        % Evaluate normalized polymerization velocity eq.3 (Mogliner and Oser 1996) for k_off = 0;
                        theta = abs(atand( Filaments.UnitVector(f,1)/Filaments.UnitVector(f,2) ));
                        theta(theta < 1 ) = 1; % theta can't be zero otherwise K has 0 in the denominator
                        %theta(theta > 85) = 85;
                        lambda = ModelParameters.PersistenceLength;
                        kT = ModelParameters.kT;
                        L = length(Filaments.MonomerIndices{f})*ModelParameters.MonomerLength;
                        L(L< 30) =  30;
                        L(L>150) = 150;
                        delta = ModelParameters.MonomerLength*cosd(theta);
                        K = 4*lambda*kT/(L.^3*sind(theta).^2); % eq. B.1
                        
                        if y <= delta && y >= 0
                            P_top = sqrt(pi*kT/(2*K))*( erfc(abs((y-delta)*sqrt(K/(2*kT)))) ); 
                            P_bot = sqrt(pi*kT/(2*K))*( erf (y*sqrt(K/(2*kT))) + 1 );
                        elseif y < 0 %-10
                            P_top = sqrt(pi*kT/(2*K))*( erfc(abs((y-delta)*sqrt(K/(2*kT)))) ); 
                            P_bot = sqrt(pi*kT/(2*K))*( erfc(abs(y*sqrt(K/(2*kT)))) );
                        else % y > delta
                            P_top = sqrt(pi*kT/(2*K))*( erf((y-delta)*sqrt(K/(2*kT))) + 1 ); 
                            P_bot = sqrt(pi*kT/(2*K))*( erf(y*sqrt(K/(2*kT))) + 1 );
                        end
                        PolymCoeff1 = P_top./P_bot;
                        
                        % Now incorporate a second polymerization coefficient
                        % From "Cell Protrusion and Retraction Driven by Fluctuations in Actin Polymerization: A Two-Dimensinoal Model"
                        % Gillian L. Ryan et al:  https://pubmed.ncbi.nlm.nih.gov/28752950/
                        PolymCoeff2 = 0.84*exp(-1e-3*y/0.5) + 0.16*exp(-1e-3*y/4); % y input is in nm
                        PolymCoeff2( PolymCoeff2 > 1) = 1; 
                        
                        % Final PolymCoefficient is the product of the two
                        PolymCoeff = PolymCoeff1 * PolymCoeff2;
    

                %|||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||
                else  
                        y = Membrane.Nodes(Membrane.Segments(m,1),2) - FilamentTips(f,2); % Distance of filament tip from membrane
                        if y < ModelParameters.MonomerLength
                            PolymCoeff = 0;
                        else
                            PolymCoeff = 1;
                        end
                end

        else % If filment tip is outside either end of the membrane set coefficient = 0
            PolymCoeff = 0;
        end
end

