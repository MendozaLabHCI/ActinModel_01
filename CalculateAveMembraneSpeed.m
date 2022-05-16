function MemVel = CalculateAveMembraneSpeed(Membrane,MembranePrevious,ModelParameters)

    if isempty(MembranePrevious)
        MemVel = 0;
    else
       AvePreMemY = mean(MembranePrevious.Nodes(MembranePrevious.Segments(:,1),2));
       AveCurMemY = mean(Membrane.Nodes(Membrane.Segments(:,1),2));
       
       MemVel = (AveCurMemY - AvePreMemY)/ModelParameters.TimeStep;
    end
end