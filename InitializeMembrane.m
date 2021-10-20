function Membrane = InitializeMembrane(MP)
   
    
    N = floor( (MP.MemleadingEdgeLengthInNanometers-MP.SegmentWidth)/(MP.SegmentWidth + MP.SpringWidth));
    
    Xs = [0; MP.SegmentWidth; (MP.SegmentWidth+MP.SpringWidth)];
    X = Xs;
    Segments = [1,2];
    Springs  = [2,3];
   
    for n = 1:N
        X = [X; X(end)+Xs(2); X(end)+Xs(3)]; 
        Segments = [Segments; [Segments(end,1)+2, Segments(end,2)+2]];
        Springs  = [Springs;  [ Springs(end,1)+2,  Springs(end,2)+2]];
    end
    X = [X; X(end)+Xs(2)]; 
    X = X - max(X)/2;
    Segments = [Segments; [Segments(end,1)+2, Segments(end,2)+2]];
    
    Membrane.Nodes = [X,MP.MemleadingEdgeHeightInNanometers*ones(length(X),1)];
    Membrane.Segments = Segments;
    Membrane.Springs = Springs;
  
    
end