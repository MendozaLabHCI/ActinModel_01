
            clc
            % Initialize Model paramaeters, filaments, membrane, and adhesions
            ModelParameters = InitializeModelParameters;
            %==========================================================================
            ModelParameters.TimeStep = 0.0002;
            ModelParameters.TotalSimulationTime = 30;
            ModelParameters.branchWindowSize = 50; % nm
            ModelParameters.VariableBranchRate = false;
            ModelParameters.BoundaryFixed = false;
            ModelParameters.SegmentWidth = 18; % 18;
            ModelParameters.SpringWidth  = 2; % 2;
            ModelParameters.MonomerLength = 2.7; %15; % 2.7; (nm) 
            ModelParameters.k_branch = 1.0; % 0.5;
            ModelParameters.FilamentsInitialLength = round(1*(15/2.7)*6); % monomers
            ModelParameters.k_off_pointed = 12;   % 12(s^-1)
            ModelParameters.k_on_barbed   = 165;   % 13 s(-1) 
            ModelParameters.k_cap = 3; %0.6; %0.6;
            ModelParameters.MembraneSpringConstant = 0.3; %0.3;
            ModelParameters.AdhesionSpringConstant = 10;
            ModelParameters.CytoplasmViscosity = 10; % Pa*s
            ModelParameters.MembraneViscosity = 100;
            NUmem = ModelParameters.MembraneViscosity*(10^-6); % Convert Pa*s to pN*s/nm^2
            ModelParameters.MembraneGamma = 6*pi*NUmem*ModelParameters.SegmentWidth; % 100 Pa*s.
            %==========================================================================
                
            ModelParameters.Adhesion_ActivationRate = 1; %0.20; % events/sec
            ModelParameters.Adhesion_DeActivationRate = 1/3; %0.1; % events/sec
            ModelParameters.Adhesion_MolecularClutchOn = true;
            ModelParameters.MolecularClutch_PeakNumber = 1;
            ModelParameters.BrownianRatchetOn = true;
            
            ModelParameters.MemleadingEdgeLengthInNanometers = 1000; %1000; %2000
            
            % Filament and adhesion density fine adjustment Scaling Factor ---------------------------
                Membrane = InitializeMembrane(ModelParameters);
                width = (Membrane.Nodes(end,1)- Membrane.Nodes(1,1))/1000;
                area1 = width*ModelParameters.AdhesionRegionDepth/1000;
            %-----------------------------------------------------------------------------------------
                
            ModelParameters.FilamentMassThreshold = round(0.5*(15/2.7)*4500*width); 
            ModelParameters.AdhesionTotalStartNumber = round(area1*2000);
            %==========================================================================
            
            ShowPlot = true;
            Membrane  = InitializeMembrane(ModelParameters);
            Filaments = InitializeActinFilaments(ModelParameters,Membrane); % (#Filaments, InitialLength, SpreadWidth, ModelParameters)
            Adhesions = InitializeAdhesions(Membrane,ModelParameters);
            FAConnections = InitializeFAConnections;
            if ShowPlot  
                [FH,AH1,AH2,AH3] = SetUpPlottingFigureAndAxes; 
            end
            delete(AH2)
            delete(AH3)
            clf(FH)
            TL = tiledlayout(FH,1,1,'TileSpacing','normal');
            AH1 = nexttile(TL,1,[1,1]);
            set(AH1,'LineWidth',2,'TickDir','in','FontSize',20)

            %set(gcf,'Position',[445,95,1459,1090])
             set(gcf,'Position',[  555          96        1317         701],'Resize','off')
            % Pre-allocate space for recording certain parameters to be saved after model completes (main saved variable is SimData)
            TimeVec = 0:ModelParameters.TimeStep:ModelParameters.TotalSimulationTime;
            nMono   = NaN(length(TimeVec),1);
            nAdhes  = NaN(length(TimeVec),1);
            MemVel  = NaN(length(TimeVec),1);
            SimData = SetUpRecordMembrane(TimeVec,Membrane);
            SimData.ModelParameters = ModelParameters;
            DATA = cell(length(TimeVec),1);
            MembranePrevious = [];
            
            % Initialize parameters to control the "nth" frame to plot
            index = 0;  
            count = 0;
            frame = 0;
            nth = 10; % nth time frame to plot
            tic
            RecordFrames = false;
            BranchesPerMicron = [];
            % Get adhesions up to "stable" level before starting model --------
            disp('Creating stable adhesion level.... please wait')
            for k = 1:round(10/ModelParameters.TimeStep)
                [Adhesions,FAConnections] = ManageAdhesions(Filaments,Adhesions,FAConnections,Membrane,ModelParameters);
            end
            
            Ymax = 0; %NaN(length(TimeVec),1);
            disp('Starting model....')
            
            % START Model ----------------------------------------------------------------------------------------------- 
            for t = TimeVec
                
                index = index + 1;
                [Filaments, Adhesions, FAConnections, PolymCoeffData] = PolymerizeDepolymerizeCapDeleteFilaments(Filaments,Adhesions,Membrane,FAConnections,ModelParameters);
                Filaments = BranchFilamentsInBranchWindowIfSelected(t,Filaments,ModelParameters,Membrane);
                [Adhesions, FAConnections] = ManageAdhesions(Filaments,Adhesions,FAConnections,Membrane,ModelParameters);
                [FAConnections, Adhesions] = FilamentAndAdhesionConnect(FAConnections,Filaments,Adhesions,ModelParameters);
                [Filaments, Membrane, Adhesions, FAConnections, Tensions, kBreaks, Data] = CalculatePositionAfterAppliedForces(Filaments,Membrane,Adhesions,FAConnections,ModelParameters);
           
                DATA{index,1}  = Data;
                nMono(index,1) = CountTotalMonomers(Filaments);
                [Filaments] = AddRandomFilaments(Filaments,Membrane,ModelParameters,nMono(index,1));

                nAdhes(index,1) = length(find(Adhesions.ActiveStatus));
                MemVel(index,1) = CalculateAveMembraneSpeed(Membrane,MembranePrevious,ModelParameters);

                %Ymax(index,1) = max(Membrane.Nodes(:,2))+100;
                Ymax = max([Ymax,mean(Membrane.Nodes(:,2))+150]);

                if t == 0.000
                    %nth = 20; 
                    RecordFrames = true;
                end
                if ShowPlot
                    %[FAConnections,count] = PlotFilamentsAndMembrane_custom01(nth,count,Filaments,Membrane,Adhesions,FAConnections,FH,AH1,AH2,AH3,t,nMono,nAdhes,MemVel,index,TimeVec,ModelParameters);
                    [FAConnections,count] = PlotFilamentsAndMembrane_custom04(Ymax,area1,nth,count,Filaments,Membrane,Adhesions,FAConnections,FH,AH1,AH2,AH3,t,nMono,nAdhes,MemVel,index,TimeVec,ModelParameters);
                    %BranchesPerMicron = [BranchesPerMicron; CalculateBranchesPerMicron(Filaments)];
                    drawnow
                    %Write snapshots of figure as tif images
                    if RecordFrames && rem(count,nth) == 0 % nth = 1; RecordFrames = true;
                        frame = frame + 1;
                        F = getframe(FH);
                        imwrite(F.cdata,...
                        ['X:\Mendoza Lab\MATLAB\Actin Growth Network Modeling\NEW_Model_Runs3\MovieData\BR-MC_500nm_tau3s_ver01\Frame_',sprintf('%04d',frame),'.tif'])
                    end
                end

                SimData.Nodes(:,:,index) = single(Membrane.Nodes); % Record membrane segment positions
                count = count + 1;
                MembranePrevious = Membrane;
                
                
          
            end
            % END Model --------------------------------------------------------------------------------------------------
toc




   % SimData.nMonomers = nMono;
   % SimData.nAdhesions = nAdhes;
   % SimData.MembraneVelocity = MemVel;
%     SimData.DATA = DATA;
%     SaveName = ['Run00',num2str(r)];
%     savefig(FH,fullfile('D:\Keith\TwoMinuteRun_001',SaveName))
%     save(fullfile('D:\Keith\TwoMinuteRun_001',[SaveName,'.mat']),'SimData','-mat','-v7.3')
    

%     try
%         save(fullfile(SaveDir,SaveName),'SimData','-mat','-v7.3')
%     catch
%         disp(['File not saved: ',SaveName])
%     end
%end

