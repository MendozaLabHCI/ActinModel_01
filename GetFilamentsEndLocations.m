function FilamentEnds = GetFilamentsEndLocations(Filaments)
    
        nF = size(Filaments.XYCoords,1);
        nF(isempty(Filaments.XYCoords)) = 0;
        FilamentEnds = zeros(nF,2);
        for f = 1:nF
            FilamentEnds(f,:) = Filaments.XYCoords{f}(1,:);
        end

end