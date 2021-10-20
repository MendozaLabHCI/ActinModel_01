function Filaments = InitializeActinFilaments(ModelParameters)
    % This function creates all the initial actin filaments
    N = ModelParameters.StartingNumberOfFilaments;
    theta = 180*rand(N,1);

    %----------------------------------------------------------
    InitialLengths =  max([ones(N,1), abs(round(ModelParameters.FilamentsInitialLength+2*randn(N,1))) ],[],2);
    Filaments.Name = (1:N)';
    Filaments.MonomerIndices = cell(N,1);
    Filaments.XYCoords = cell(N,1);
    Filaments.UnitVector = zeros(N,2);
    D = ModelParameters.MonomerLength; %(nm)
    
    for f = 1:N
        % Create initial two point of filament as a horizontal line (origin is implied as second point).
            x = ModelParameters.MonomerLength;
            y = 0;
        % Now randomly rotate between 0 and 180 degrees
            R = [ [cosd(theta(f,1)), -sind(theta(f,1))];... % create rotation matrix
                  [sind(theta(f,1)),  cosd(theta(f,1))] ];
             
            M = R*[x;y]; 
            xr = M(1); % Where second point is the end closest to the membrane
            yr = M(2); % Where second point is the end closest to the membrane
            
            Filaments.MonomerIndices{f,1}  = 1; %% Index values for the first point of the filament
            Filaments.XYCoords{f,1} = [0,0]; % coordinate of first point
            Filaments.UnitVector(f,:) = [xr,yr]./vecnorm([xr;yr]);
            % plot(x,y,'-b',xr,yr,'-r',x(2),y(2),'.b',xr(2),yr(2),'.r'); axis equal; axis([-3,3,-3,3]); pause
            
            % Make Filaments the desired length
            for k = 1:InitialLengths(f)
                Filaments.MonomerIndices{f} = [Filaments.MonomerIndices{f}; Filaments.MonomerIndices{f}(end,1)+1;]; 
                Filaments.XYCoords{f} = [Filaments.XYCoords{f}; [Filaments.XYCoords{f}(end,1) + D*Filaments.UnitVector(f,1),...
                                                                 Filaments.XYCoords{f}(end,2) + D*Filaments.UnitVector(f,2)] ];
            end
            
            % Make the filament's center of mass at [0,0] 
            Filaments.XYCoords{f} = Filaments.XYCoords{f} - Filaments.XYCoords{f}(end,:)/2;
    end
    
    % Setup starting/default parameters
    Filaments.IsCapped = false(N,1);          % Has the filament been capped?
    Filaments.MainIndex = (1:N)';             % Main Filament group that this filament is a part of. (based on parent filament name)
    Filaments.Parent = zeros(N,1);            % The filament parent name that each filament branched off of (zero if it's a main filament)
    Filaments.ParentIndex  = zeros(N,1);      % The Monomer Index on the parent filament where this filament branched from.
  
    
    % Spread out Filaments in the x direction and apply y-offset --------------------------------------
    Xrand = ModelParameters.SpreadWidth*(rand(N,1)-0.5);
    for f = 1:N
         Filaments.XYCoords{f}(:,1) = Filaments.XYCoords{f}(:,1) + Xrand(f,1);
         Filaments.XYCoords{f}(:,2) = Filaments.XYCoords{f}(:,2) + ModelParameters.VerticalOffSet;
    end
    %------------------------------------------------------------
end 


