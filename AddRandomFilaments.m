function [Filaments] = AddRandomFilaments(Filaments,Membrane,ModelParameters,nMonomers)
    

    if nMonomers < ModelParameters.FilamentMassThreshold
        
        InitialLength = abs(round(ModelParameters.FilamentsInitialLength+2*randn(1)));
        theta = 45+90*rand;%%180*rand(1);
        %disp(theta)
        D = ModelParameters.MonomerLength; %(nm)
        x = ModelParameters.MonomerLength;
        y = 0;

        % Now randomly rotate between 0 and 180 degrees
        R = [ [cosd(theta), -sind(theta)];... % create rotation matrix
              [sind(theta),  cosd(theta)] ];

        M = R*[x;y]; 
        xr = M(1); % Where second point is the end closest to the membrane
        yr = M(2); % Where second point is the end closest to the membrane

        MonomerIndices  = 1; % Index values for the first point of the filament
        XYCoords = [0,0];    % coordinate of first point
        UnitVector = [xr,yr]./vecnorm([xr;yr]);
        %plot(x,y,'-b',xr,yr,'-r',x(2),y(2),'.b',xr(2),yr(2),'.r'); axis equal; axis([-3,3,-3,3]); pause

        % Make Filaments the desired length
        for k = 1:InitialLength
            MonomerIndices = [MonomerIndices; MonomerIndices(end,1)+1;]; 
            XYCoords = [XYCoords; [XYCoords(end,1) + D*UnitVector(1,1),...
                                   XYCoords(end,2) + D*UnitVector(1,2)] ];
        end

        % Make the filament's center of mass at [0,0]
        XYCoords = XYCoords - XYCoords(end,:)/2;
        SpreadWidth = Membrane.Nodes(end,1) - Membrane.Nodes(1,1);
        Xrand = (SpreadWidth)*(rand(1)-0.5); % Calculate random horizonal position
        XYCoords(:,1) = XYCoords(:,1) + Xrand;
        % Find membrane segment closest to this new horizontal position and use it's y coordinate to reference where to move new filament to.
        MembraneSegmentCenters = mean([ Membrane.Nodes(Membrane.Segments(:,1),1) ,...
                                        Membrane.Nodes(Membrane.Segments(:,2),1) ],2);
        [~,idx] = min( abs(MembraneSegmentCenters - XYCoords(end,1)) );
        
        %XYCoords(:,1) = XYCoords(:,1) + Xrand;
        m = find( Membrane.Nodes(Membrane.Segments(:,1),1) - ModelParameters.SpringWidth/2  <=  XYCoords(end,1) &... % Left segment endpoint  <= Filament X-end
                  Membrane.Nodes(Membrane.Segments(:,2),1) + ModelParameters.SpringWidth/2  >=  XYCoords(end,1) );   % Right segment endpoint >= Filament X-end 
        
        if isempty(m); m = 1; end
            
        XYCoords(:,2) = XYCoords(:,2) - max(XYCoords(:,2)) + Membrane.Nodes(m,2) - 20; % Set filament tip y-position -20 below closest membrane segment
    
        NewFilamentName = max(Filaments.Name) + 1; 
        
        Filaments.Name           = [ Filaments.Name; NewFilamentName];
        Filaments.MonomerIndices = [ Filaments.MonomerIndices; {MonomerIndices} ];
        Filaments.XYCoords       = [ Filaments.XYCoords; {XYCoords}];
        Filaments.UnitVector     = [ Filaments.UnitVector; UnitVector ];
        Filaments.IsCapped       = [ Filaments.IsCapped; false ];
        Filaments.MainIndex      = [ Filaments.MainIndex; NewFilamentName ];
        Filaments.Parent         = [ Filaments.Parent; 0 ];
        Filaments.ParentIndex    = [ Filaments.ParentIndex; 0 ];

    end
end

