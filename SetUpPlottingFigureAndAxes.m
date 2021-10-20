function [FH,AH1,AH2,AH3] = SetUpPlottingFigureAndAxes

    FH = figure(1); 
    clf(FH)
    set(FH,'Position',[500,56,1500,1100],'Color',[1,1,1])
    AH1 = axes('Parent',FH,'Position',[0.07,0.30,0.88,0.7]);
    %AH2 = axes('Parent',FH,'Position',[0.07,0.08,0.88,0.20]);
    AH2 = axes('Parent',FH,'Position',[0.07,0.08,0.39,0.20]);
    AH3 = axes('Parent',FH,'Position',[0.56,0.08,0.39,0.20]);
end
