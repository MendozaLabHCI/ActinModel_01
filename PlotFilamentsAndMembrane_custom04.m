function [FAConnections,count] = PlotFilamentsAndMembrane_custom04(Ymax,area1,nth,count,Filaments,Membrane,Adhesions,FAConnections,FH,AH1,AH2,AH3,t,nMono,nAdhes,MemVel,index,TimeVec,ModelParameters)

    if (count >= nth) || isequal(t,0)% Plot every 10th, 100th, or nth frame etc.
        count = 0;
        cla(AH1)

        idxMF = find(Filaments.Parent == 0); % Find Main Filaments (Filaments that are at the top of the parent/daughter structure)
        nC = length(idxMF);
        C = lines(nC); % Create colors for each of the attached groups of filaments
        
        % For all the attached groups of filaments, plot each group with its own color
        TopLayerIndices = [];
        
        for MF = 1:nC 
            if rem(Filaments.Name(idxMF(MF)),10) ~= 0
                FColor = [0.75,0.75,0.75];
                LWidth = 1;
%                 if isequal(MF,1)
%                     hold(AH1,'on')
%                 end
                % Find all the filaments that are attached to the main parent
                % filament, and the main parent filament as well.
                idx = find(Filaments.MainIndex == Filaments.Name(idxMF(MF)));
                for f = idx'
                    if isequal(Filaments.Parent(f,1),0)
                        % Plotting for Main Filaments
    %                     plot(AH1, Filaments.XYCoords{f}(:,1),...
    %                               Filaments.XYCoords{f}(:,2),'-','Color',C(MF,:),'LineWidth',1)
                          plot(AH1, Filaments.XYCoords{f}(:,1),...
                                    Filaments.XYCoords{f}(:,2),'-','Color',FColor,'LineWidth',LWidth)
                           hold(AH1,'on')
                    else
                        % Plotting for daughter filaments
                        idxP  = find(Filaments.Name == Filaments.Parent(f,1));
                        idxXY = find(Filaments.MonomerIndices{idxP} == Filaments.ParentIndex(f,1));
                        XYstart = Filaments.XYCoords{idxP}(idxXY,:);
    %                     plot(AH1, [XYstart(1,1); Filaments.XYCoords{f}(:,1)],...
    %                               [XYstart(1,2); Filaments.XYCoords{f}(:,2)],'-','Color',C(MF,:),'LineWidth',1)
                        plot(AH1, [XYstart(1,1); Filaments.XYCoords{f}(:,1)],...
                                  [XYstart(1,2); Filaments.XYCoords{f}(:,2)],'-','Color',FColor,'LineWidth',LWidth)
                         hold(AH1,'on')
                    end    
                end
                    
            else
                TopLayerIndices = [TopLayerIndices,MF];
            end
        end


        for MF = TopLayerIndices 
                FColor = [0,0,1];
                LWidth = 2;
%                 if isequal(MF,1)
%                     hold(AH1,'on')
%                 end
                % Find all the filaments that are attached to the main parent
                % filament, and the main parent filament as well.
                idx = find(Filaments.MainIndex == Filaments.Name(idxMF(MF)));
                for f = idx'
                    if isequal(Filaments.Parent(f,1),0)
                        % Plotting for Main Filaments
    %                     plot(AH1, Filaments.XYCoords{f}(:,1),...
    %                               Filaments.XYCoords{f}(:,2),'-','Color',C(MF,:),'LineWidth',1)
                          plot(AH1, Filaments.XYCoords{f}(:,1),...
                                    Filaments.XYCoords{f}(:,2),'-','Color',FColor,'LineWidth',LWidth)
                    else
                        % Plotting for daughter filaments
                        idxP  = find(Filaments.Name == Filaments.Parent(f,1));
                        idxXY = find(Filaments.MonomerIndices{idxP} == Filaments.ParentIndex(f,1));
                        XYstart = Filaments.XYCoords{idxP}(idxXY,:);
    %                     plot(AH1, [XYstart(1,1); Filaments.XYCoords{f}(:,1)],...
    %                               [XYstart(1,2); Filaments.XYCoords{f}(:,2)],'-','Color',C(MF,:),'LineWidth',1)
                        plot(AH1, [XYstart(1,1); Filaments.XYCoords{f}(:,1)],...
                                  [XYstart(1,2); Filaments.XYCoords{f}(:,2)],'-','Color',FColor,'LineWidth',LWidth)
                    end    
                end
            
        end
        
        
        set(gca,'LineWidth',1.5)
        
        % Plot Red dot for capped filaments
        if any(Filaments.IsCapped)
            idxC = find(Filaments.IsCapped);
            for f = 1:length(idxC)
                Xc = Filaments.XYCoords{idxC(f)}(end,1);
                Yc = Filaments.XYCoords{idxC(f)}(end,2);
                plot(AH1, Xc,Yc,'.','MarkerSize',9,'Color',[0,0.6,0])
            end
        end 
        
       
       % Plot adhesions regions -------------------------------
%        for m = 1:Nseg
%             XY_LT = Adhesions.RegionNodes(1,:,m);
%             XY_RT = Adhesions.RegionNodes(2,:,m);
%             XY_RB = Adhesions.RegionNodes(3,:,m);
%             XY_LB = Adhesions.RegionNodes(4,:,m);
%             plot(AH1,[XY_LT(1); XY_LB(1); XY_RB(1); XY_RT(1)],...
%                      [XY_LT(2); XY_LB(2); XY_RB(2); XY_RT(2)], '-', 'Color', [0.9,0.9,0.9])
%        end

        % Calculate adhesion regions area for density calculation --------------------------
       delX = Membrane.Nodes(end,1) -  Membrane.Nodes(1,1);
       delY = ModelParameters.AdhesionRegionDepth;
       AreaOfAllRegions = delX*delY/1000^2; % um^2

       % Plot active adhesion indices and plot
       idx = find(Adhesions.ActiveStatus);
       %[~,~,val] = find(idx ~= FAConnections.AdhesionIndex);
       plot(AH1, Adhesions.XYPoints(idx,1), Adhesions.XYPoints(idx,2), '*k', 'MarkerSize', 5)
       
       % plot attached adhesions as red
       idx = FAConnections.AdhesionIndex;
       plot(AH1, Adhesions.XYPoints(idx,1),Adhesions.XYPoints(idx,2),  '*r', 'MarkerSize', 5)
       AD = round( (nAdhes( find(~isnan(nAdhes),1,'last'))/area1) );

       if logical(exist('TH'))
            if isgraphics(TH); delete(TH); end
       end
       TH = text(AH1, 0.01, 1.13,...
                          [{'\bfBR-MC model    \fontsize{16}\rm\it\tau_{max}\rm = 3.0 s    \itk_{m}\rm = 0.3 pN/nm    \itk_{B}\rm = 10 pN/nm    dt = 0.2 ms'},...
                           {['\fontsize{20}\bfAdhesion density = ',char(pad(string(sprintf('%0.0f', AD )),4,'left')),...
                             ' \mum^{-2}                    Time = ', char(pad(string(sprintf('%#0.4f',t)),7,'left')), ' s   ']}],...
                          'FontName','monospaced', 'FontSize', 20, 'FontWeight','bold','Units','normalized');

%        title(AH1,[{'\bfBR-MC     \rm\it\tau_{max}\rm = 3.0 s    \itk_{B}\rm = 10 pN/nm     dt = 0.2 ms'},...
%                   {['\bfAdhesion density = ',char(pad(string(sprintf('%0.0f', AD )),4,'left')),...
%                   ' \mum^{-2}            Time = ', char(pad(string(sprintf('%#0.3f',t)),6,'left')), ' s   ']}],...
%                   'FontName','monospaced', 'FontSize', 24, 'FontWeight','bold')
%        AH1.TitleHorizontalAlignment = 'left';

              
      
      

       % Plot Adhesions connections --------------------------
       idx = find(~isnan(FAConnections.AdhesionIndex));
       for n = idx'
                Axy = Adhesions.XYPoints(FAConnections.AdhesionIndex(n,1),:); 
                f = find(Filaments.Name == FAConnections.FilamentName(n,1)); 
                Midx = find(Filaments.MonomerIndices{f} == FAConnections.MonomerIndex(n,1));
                Mxy = Filaments.XYCoords{f}(Midx,:);
                plot(AH1, [Axy(1);Mxy(1)], [Axy(2);Mxy(2)],':k')
       end
       
       
       % Plot Membrane ------------------------------------
            Nseg = size(Membrane.Segments,1);
            Nspr = size(Membrane.Springs,1);
            % Plot Segments ---------------------------------------------------
            for n = 1:Nseg
                plot( AH1, [Membrane.Nodes(Membrane.Segments(n,1),1);Membrane.Nodes(Membrane.Segments(n,2),1)],...
                           [Membrane.Nodes(Membrane.Segments(n,1),2); Membrane.Nodes(Membrane.Segments(n,2),2)], 'k-','LineWidth',3)
                if isequal(n,1)
                    hold(AH1,'on'); 
                end
            end
            % Plots Springs ---------------------------------------------------
            for n = 1:Nspr
                plot( AH1, [Membrane.Nodes(Membrane.Springs(n,1),1);Membrane.Nodes(Membrane.Springs(n,2),1)],...
                           [Membrane.Nodes(Membrane.Springs(n,1),2); Membrane.Nodes(Membrane.Springs(n,2),2)], 'k-')
            end
      % Plot Boundary Force region ------------------------------------
            Nseg = size(Membrane.Segments,1);
            Nspr = size(Membrane.Springs,1);
            offset = 10;% nm
            % Plot Segments ---------------------------------------------------
            for n = 1:Nseg
                plot( AH1, [Membrane.Nodes(Membrane.Segments(n,1),1);Membrane.Nodes(Membrane.Segments(n,2),1)],...
                           [Membrane.Nodes(Membrane.Segments(n,1),2) + offset; Membrane.Nodes(Membrane.Segments(n,2),2) + offset], 'k:','LineWidth',2)
                if isequal(n,1)
                    hold(AH1,'on'); 
                end
            end
            % Plots Springs ---------------------------------------------------
            for n = 1:Nspr
                plot( AH1, [Membrane.Nodes(Membrane.Springs(n,1),1);Membrane.Nodes(Membrane.Springs(n,2),1)],...
                           [Membrane.Nodes(Membrane.Springs(n,1),2) + offset; Membrane.Nodes(Membrane.Springs(n,2),2) + offset], 'k-','LineWidth',2)
            end
      %---------------------------------------------------
      % Plot Contact Region ------------------------------------
%             Nseg = size(Membrane.Segments,1);
%             Nspr = size(Membrane.Springs,1);
%             offset = -15;% nm
%             % Plot Segments ---------------------------------------------------
%             for n = 1:Nseg
%                 plot( AH1, [Membrane.Nodes(Membrane.Segments(n,1),1);Membrane.Nodes(Membrane.Segments(n,2),1)],...
%                            [Membrane.Nodes(Membrane.Segments(n,1),2) + offset; Membrane.Nodes(Membrane.Segments(n,2),2) + offset], 'k:','LineWidth',2)
%                 if isequal(n,1)
%                     hold(AH1,'on'); 
%                 end
%             end
%             % Plots Springs ---------------------------------------------------
%             for n = 1:Nspr
%                 plot( AH1, [Membrane.Nodes(Membrane.Springs(n,1),1);Membrane.Nodes(Membrane.Springs(n,2),1)],...
%                            [Membrane.Nodes(Membrane.Springs(n,1),2) + offset; Membrane.Nodes(Membrane.Springs(n,2),2) + offset], 'k-','LineWidth',2)
%             end
      %---------------------------------------------------
      
      
       hold(AH1,'off')
       axis(AH1,'equal') 
       box(AH1, 'on')
       grid(AH1,'on')%'Color',[0.8,0.8,0.8])
       padding = 100; % 0.1*( Membrane.Nodes(end,1) -  Membrane.Nodes(1,1));
       PositionVector = [( Membrane.Nodes(1,1) - padding) ( Membrane.Nodes(end,1) + padding) Ymax-500 Ymax];
       axis(AH1,PositionVector);
       yticks(AH1, round(Ymax,-2)-800:100:round(Ymax,-2)+100)
       yticklabels(compose('%0.1f',yticks'./1000)')
       xticks(gca,-500:200:500)
       xticklabels(compose('%0.1f',(0:0.2:1)'))
       %set(AH1,'Position',[0.10124   0.20846   0.803755   0.61806])
       set(AH1,'FontSize',18)
       xlabel(AH1,'X(\mum)','FontSize',22,'FontWeight','bold')
       ylabel(AH1,'Y(\mum)','FontSize',22,'FontWeight','bold')
       %disp([diff(ylim),diff(xlim)])
       %disp(PositionVector)
    end
end



