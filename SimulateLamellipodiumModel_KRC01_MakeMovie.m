
            clc
            % Initialize Model paramaeters, filaments, membrane, and adhesions
            ModelParameters = InitializeModelParameters;
            %==========================================================================
                ModelParameters.TotalSimulationTime = 61;
                ModelParameters.TimeStep = 0.001;
               
                ModelParameters.MemleadingEdgeLengthInNanometers = 1000; %2000
                ModelParameters.Adhesion_ActivationRate = 1.0; % events/sec
                ModelParameters.Adhesion_DeActivationRate = 0.1; % events/sec
                ModelParameters.Adhesion_MolecularClutchOn = true;
                ModelParameters.MolecularClutch_PeakNumber = 1;
                ModelParameters.BrownianRatchetOn = true;
                ModelParameters.BoundaryFixed = false;
                
                % Filament and adhesion density fine adjustment Scaling Factor ---------------------------
                    Membrane = InitializeMembrane(ModelParameters);
                        width = (Membrane.Nodes(end,1)- Membrane.Nodes(1,1))/1000;
                        area1 = width*ModelParameters.AdhesionRegionDepth/1000;
                %-----------------------------------------------------------------------------------------
                
                ModelParameters.FilamentMassThreshold = round(4500*width); %9000
                
                if ModelParameters.Adhesion_MolecularClutchOn
                    ModelParameters.AdhesionTotalStartNumber = round(area1*2000); % 800 MC = on, 400 MC = off; (for LeadingEdgeLength = 2000nm)
                else
                    ModelParameters.AdhesionTotalStartNumber = round(area1*1000);
                end
                
                
                
                ModelParameters.SegmentWidth = 18; % 18;
                ModelParameters.MembraneGamma =  0.05;%0.05;
                ModelParameters.SpringWidth  = 2; % 2;
                ModelParameters.MembraneSpringConstant = 0.3;
                ModelParameters.AdhesionSpringConstant = 10;
                ModelParameters.BoundaryForceSpringConstant = 10;
            %==========================================================================
            
            ShowPlot = true;
            Membrane  = InitializeMembrane(ModelParameters);
            Filaments = InitializeActinFilaments(ModelParameters,Membrane); % (#Filaments, InitialLength, SpreadWidth, ModelParameters)
            Adhesions = InitializeAdhesions(Membrane,ModelParameters);
            FAConnections = InitializeFAConnections;
            if ShowPlot  
                [FH,AH1,AH2,AH3] = SetUpPlottingFigureAndAxes; 
            end
            set(gcf,'Position',[445,95,1459,1090])
            
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
            nth = 50; % nth time frame to plot
            tic
            RecordFrames = false;
            BranchesPerMicron = [];
            % Get adhesions up to "stable" level before starting model --------
            disp('Creating stable adhesion level.... please wait')
            for k = 1:round(10/ModelParameters.TimeStep)
                [Adhesions,FAConnections] = ManageAdhesions(Filaments,Adhesions,FAConnections,Membrane,ModelParameters);
            end
            
            disp('Starting model....')
            
            % START Model ----------------------------------------------------------------------------------------------- 
            for t = TimeVec
                
                index = index + 1;
                [Filaments, Adhesions, FAConnections, PolymCoeffData] = PolymerizeDepolymerizeCapDeleteFilaments(Filaments,Adhesions,Membrane,FAConnections,ModelParameters);
                Filaments = BranchFilamentsInBranchWindowIfSelected(Filaments,ModelParameters,Membrane);
                [Adhesions, FAConnections] = ManageAdhesions(Filaments,Adhesions,FAConnections,Membrane,ModelParameters);
                [FAConnections, Adhesions] = FilamentAndAdhesionConnect(FAConnections,Filaments,Adhesions,ModelParameters);
                [Filaments, Membrane, Adhesions, FAConnections, Tensions, kBreaks, Data] = CalculatePositionAfterAppliedForces(Filaments,Membrane,Adhesions,FAConnections,ModelParameters);
        
                DATA{index,1}  = Data;
                nMono(index,1) = CountTotalMonomers(Filaments);
                [Filaments] = AddRandomFilaments(Filaments,Membrane,ModelParameters,nMono(index,1));

                nAdhes(index,1) = length(find(Adhesions.ActiveStatus));
                MemVel(index,1) = CalculateAveMembraneSpeed(Membrane,MembranePrevious,ModelParameters);
                if t == 0.000
                    nth = 20; 
                    RecordFrames = true;
                end
                if ShowPlot
                    %[FAConnections,count] = PlotFilamentsAndMembrane_custom01(nth,count,Filaments,Membrane,Adhesions,FAConnections,FH,AH1,AH2,AH3,t,nMono,nAdhes,MemVel,index,TimeVec,ModelParameters);
                    [FAConnections,count] = PlotFilamentsAndMembrane_custom04(area1,nth,count,Filaments,Membrane,Adhesions,FAConnections,FH,AH1,AH2,AH3,t,nMono,nAdhes,MemVel,index,TimeVec,ModelParameters);
                    %BranchesPerMicron = [BranchesPerMicron; CalculateBranchesPerMicron(Filaments)];
                    drawnow
                    %Write snapshots of figure as tif images
                    if RecordFrames && rem(count,nth) == 0 % nth = 1; RecordFrames = true;
                        frame = frame + 1;
                        F = getframe(FH);
                        imwrite(F.cdata,...
                        ['X:\Mendoza Lab\MATLAB\Actin Growth Network Modeling\NEW_Model_Runs2\Movies\BR-MC_Tau3_Kb10_60sec_Movie01\Frame_',sprintf('%04d',frame),'.tif'])
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

