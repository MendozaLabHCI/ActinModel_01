function FilamentTips = GetFilamentsTipLocations(Filaments)
    
        nF = size(Filaments.XYCoords,1);
        nF(isempty(Filaments.XYCoords)) = 0;
        FilamentTips = zeros(nF,2);
        for f = 1:nF
            FilamentTips(f,:) = Filaments.XYCoords{f}(end,:);
        end

end