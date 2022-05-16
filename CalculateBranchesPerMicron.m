function BranchesPerMicron = CalculateBranchesPerMicron(Filaments)

    BranchesPerMicron = [];

    idx1 = find( Filaments.Parent == 0 )';
    for f = idx1
        Name = Filaments.Name(f,1);
        idx2 = find(Filaments.Parent == Name); % Find all the filaments that branched off Main Filament: Name
        FilamentLength = sqrt( (Filaments.XYCoords{f,1}(1,1) - Filaments.XYCoords{f,1}(end,1))^2 + (Filaments.XYCoords{f,1}(1,2) - Filaments.XYCoords{f,1}(end,2))^2 );
        if FilamentLength ~= 0
            BranchesPerMicron = [BranchesPerMicron; length(idx2)/(FilamentLength/1000)];
        end
    end

    idx3 = find(Filaments.Parent ~= 0)';
    Parents = unique( Filaments.Parent(idx3) );
    for n = 1:length(Parents)
        f = find(Filaments.Name == Parents(n,1));
        FilamentLength = sqrt( (Filaments.XYCoords{f,1}(1,1) - Filaments.XYCoords{f,1}(end,1))^2 + (Filaments.XYCoords{f,1}(1,2) - Filaments.XYCoords{f,1}(end,2))^2 );
        idx4 = find( Filaments.Parent == Parents(n,1) );
        if FilamentLength ~= 0
            BranchesPerMicron = [BranchesPerMicron; length(idx4)/(FilamentLength/1000)];
        end
    end