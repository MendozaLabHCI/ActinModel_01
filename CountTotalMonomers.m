function nMonomers = CountTotalMonomers(Filaments)

    nF = length(Filaments.MonomerIndices);
    nMonomers = 0;
    for f = 1:nF
        nMonomers = nMonomers + length(Filaments.MonomerIndices{f});
    end

end