function [FH,AH1,AH2,AH3] = SetUpPlottingFigureAndAxes

    FH = figure(20); 
    clf(FH)
    set(FH,'Position',[ 543,61,1471,1100],'Color',[1,1,1])
    AH1 = axes('Parent',FH,'Position',[0.07,0.30,0.88,0.7]);
    %AH2 = axes('Parent',FH,'Position',[0.07,0.08,0.88,0.20]);
    AH2 = axes('Parent',FH,'Position',[0.07,0.08,0.39,0.20]);
    AH3 = axes('Parent',FH,'Position',[0.56,0.08,0.39,0.20]);
end
