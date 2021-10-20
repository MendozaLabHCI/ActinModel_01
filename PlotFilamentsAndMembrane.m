function [FAConnections,count] = PlotFilamentsAndMembrane(nth,count,Filaments,Membrane,Adhesions,FAConnections,FH,AH1,AH2,AH3,t,nMono,nAdhes,MemVel,index,TimeVec,ModelParameters)

    if (count >= nth) || isequal(t,0)% Plot every 10th, 100th, or nth frame etc.
        count = 0;
        figure(FH)
        cla(AH1)
        cla(AH2)
        % Plot Membrane ------------------------------------
        Nseg = size(Membrane.Segments,1);
        Nspr = size(Membrane.Springs,1);
        %---------------------------------------------------
        for n = 1:Nseg
            plot( AH1, [Membrane.Nodes(Membrane.Segments(n,1),1);Membrane.Nodes(Membrane.Segments(n,2),1)],...
                       [Membrane.Nodes(Membrane.Segments(n,1),2);Membrane.Nodes(Membrane.Segments(n,2),2)], 'k-','LineWidth',3)
            if isequal(n,1)
                hold(AH1,'on'); 
            end
        end
        %---------------------------------------------------
        for n = 1:Nspr
            plot( AH1, [Membrane.Nodes(Membrane.Springs(n,1),1);Membrane.Nodes(Membrane.Springs(n,2),1)],...
                       [Membrane.Nodes(Membrane.Springs(n,1),2);Membrane.Nodes(Membrane.Springs(n,2),2)], 'k-')
        end
        %---------------------------------------------------
        idxMF = find(Filaments.Parent == 0); % Find Main Filaments (Filaments that are at the top of the parent/daughter structure)
        nC = length(idxMF);
        C = lines(nC); % Create colors for each of the attached groups of filaments

        % For all the attached groups of filaments, plot each group with its own color
        for MF = 1:nC 
            if isequal(MF,1)
                hold(AH1,'on')
            end
            % Find all the filaments that are attached to the main parent
            % filament, and the main parent filament as well.
            idx = find(Filaments.MainIndex == Filaments.Name(idxMF(MF)));
            for f = idx'
                if isequal(Filaments.Parent(f,1),0)
                    % Plotting for Main Filaments
                    plot(AH1, Filaments.XYCoords{f}(:,1),...
                              Filaments.XYCoords{f}(:,2),'-','Color',C(MF,:),'LineWidth',2)
                else
                    % Plotting for daughter filaments
                    idxP  = find(Filaments.Name == Filaments.Parent(f,1));
                    idxXY = find(Filaments.MonomerIndices{idxP} == Filaments.ParentIndex(f,1));
                    XYstart = Filaments.XYCoords{idxP}(idxXY,:);
                    plot(AH1, [XYstart(1,1); Filaments.XYCoords{f}(:,1)],...
                              [XYstart(1,2); Filaments.XYCoords{f}(:,2)],'-','Color',C(MF,:),'LineWidth',2)
                end    
            end
        end

        % Plot Red dot for capped filaments
        if any(Filaments.IsCapped)
            idxC = find(Filaments.IsCapped);
            for f = 1:length(idxC)
                Xc = Filaments.XYCoords{idxC(f)}(end,1);
                Yc = Filaments.XYCoords{idxC(f)}(end,2);
                plot(AH1, Xc,Yc,'.r','MarkerSize',8)
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

        % Calculate adhesion density --------------------------
       delX = Adhesions.RegionNodes(2,1,1) - Adhesions.RegionNodes(1,1,1);
       delY = Adhesions.RegionNodes(1,2,1) - Adhesions.RegionNodes(4,2,1);
       AreaOfAllRegions = delX*delY*size(Adhesions.RegionNodes,3)/1000^2; % um^2

       % Plot active adhesion indices and plot
       idx = find(Adhesions.ActiveStatus);
       %[~,~,val] = find(idx ~= FAConnections.AdhesionIndex);
       plot(AH1, Adhesions.XYPoints(idx,1), Adhesions.XYPoints(idx,2), '*k', 'MarkerSize', 6)
       
       % plot attached adhesions as red
       idx = FAConnections.AdhesionIndex;
       plot(AH1, Adhesions.XYPoints(idx,1),Adhesions.XYPoints(idx,2),  '*r', 'MarkerSize', 6)
       
       title(AH1,['Time = ', sprintf('%#0.3f',t), ' sec     Adhesion Density (\mum^{-2}) = ',... 
                  num2str(nAdhes(index,1)/AreaOfAllRegions,'%#3.1f'),...
                  '     Adhesion-Filament Connections = ', num2str(size(FAConnections.AdhesionIndex,1)),... 
                  '     dt = ', num2str(ModelParameters.TimeStep),' sec'],...
                  'FontName','monospaced','FontSize', 14, 'FontWeight','bold')
              
       MemYmax = ceil(max(Membrane.Nodes(:,2)));
       axis(AH1,'equal') 
       %axis(AH1,[-ModelParameters.MemleadingEdgeLengthInNanometers/2 ModelParameters.MemleadingEdgeLengthInNanometers/2 (MemYmax-800) (MemYmax+200)])
       axis(AH1,[-1000 1000 -300 700])
       box(AH1, 'on')
       xlabel(AH1,'X(nm)','FontSize',14)
       ylabel(AH1,'Y(nm)','FontSize',14)
       

       % Plot Adhesions connections --------------------------
       idx = find(~isnan(FAConnections.AdhesionIndex));
       for n = idx' %n = 1:nC
           try
                Axy = Adhesions.XYPoints(FAConnections.AdhesionIndex(n,1),:); 
                f = find(Filaments.Name == FAConnections.FilamentName(n,1)); 
                Midx = find(Filaments.MonomerIndices{f} == FAConnections.MonomerIndex(n,1));
                Mxy = Filaments.XYCoords{f}(Midx,:);
                plot(AH1, [Axy(1);Mxy(1)], [Axy(2);Mxy(2)],':k')
           catch
               FAConnections.AdhesionIndex(n) = NaN;
               FAConnections.FilamentName(n) = NaN;
               FAConnections.MonomerIndex(n) = NaN;
           end
       end
       
       hold(AH1,'off')
       
       
       
       % Plot Total Monomers
       yyaxis(AH2,'left')
           plot(AH2, TimeVec, nMono,'-b')
           xlim(AH2,[0,max(TimeVec)])
           ylim(AH2,[0,max(nMono)+100])
           xlabel(AH2,'Time (s)','FontSize',14)
           ylabel(AH2,[{'Filament Mass'};{'(Monomers)'}],'FontSize',14)
       yyaxis(AH2,'right')
           plot(AH2, TimeVec, nAdhes,'-r')
           ylim(AH2,[0,max(nAdhes)+10])
           ylabel(AH2,'nAdhesions','FontSize',14)
       drawnow
       
       % Plot Membrane Speed
       idx = find(~isnan(MemVel));
       VelSmooth = NaN(size(MemVel));
       VelSmooth(idx,1) = smooth(MemVel(idx),1/ModelParameters.TimeStep,'sgolay');
       plot(AH3,TimeVec,MemVel,'-','Color',[0.8,0.8,0.8]); hold(AH3,'on')
       plot(AH3,TimeVec,VelSmooth,'r-','LineWidth',2);     hold(AH3,'off')
       xlim(AH3,[0,max(TimeVec)])
       %ylim(AH3,[0,max(nMono)+10])
       xlabel(AH3,'Time (s)','FontSize',14)
       ylabel(AH3,'Membrane Speed (nm/s)','FontSize',14)
    end
end



