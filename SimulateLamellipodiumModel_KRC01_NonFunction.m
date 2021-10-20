%function SimData = SimulateLamelipodiumModel_KRC01(ModelParameters,SaveDir,SaveName,ShowPlot)
            clc
            % Initialize Model paramaeters, filaments, membrane, and adhesions
            ModelParameters = InitializeModelParameters;
            ModelParameters.FilamentMassThreshold = 2000; %9000;
            ModelParameters.TotalSimulationTime = 100;
            ModelParameters.Adhesion_OFF_DependentOnAdhesionFilamentTension = true;
            ModelParameters.AllowAdhesionFilamentBondBreak = true;
            ModelParameters.BoundaryFixed = false;
            ModelParameters.MembraneSpringConstant = 0.3;
            
            ShowPlot = true;
            Filaments = InitializeActinFilaments(ModelParameters); % (#Filaments, InitialLength, SpreadWidth, ModelParameters)
            Membrane  = InitializeMembrane(ModelParameters);
            Adhesions = InitializeAdhesions(Membrane,ModelParameters);
            FAConnections = InitializeFAConnections;
            if ShowPlot  
                [FH,AH1,AH2,AH3] = SetUpPlottingFigureAndAxes; 
            end
            
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
            nth = 100; % nth time frame to plot
            tic
            
            % Get adhesions up to "stable" level before starting model --------
            disp('Creating stable adhesion level.... please wait')
            for k = 1:round(10/ModelParameters.TimeStep)
                [Adhesions,FAConnections] = ManageAdhesions(Filaments,Adhesions,FAConnections,Membrane,ModelParameters);
            end
            
            frame = 0;
            
            disp('Starting model....')
            % START Model ----------------------------------------------------------------------------------------------- 
            for t = TimeVec
                
                index = index + 1;
                [Filaments,Adhesions,FAConnections] = PolymerizeDepolymerizeCapDeleteFilaments(Filaments,Adhesions,Membrane,FAConnections,ModelParameters);
                Filaments = BranchFilamentsInBranchWindowIfSelected(Filaments,ModelParameters,Membrane);
                [Adhesions,FAConnections]  = ManageAdhesions(Filaments,Adhesions,FAConnections,Membrane,ModelParameters);
                [FAConnections, Adhesions] = FilamentAndAdhesionConnect(FAConnections,Filaments,Adhesions,ModelParameters);
                [Filaments, Membrane, Adhesions, FAConnections, Tensions, kBreaks, Data] = ...
                                                  CalculatePositionAfterAppliedForces(Filaments,Membrane,Adhesions,FAConnections,ModelParameters);

                DATA{index,1}  = Data;
                nMono(index,1) = CountTotalMonomers(Filaments);
                [Filaments] = AddRandomFilaments(Filaments,Membrane,ModelParameters,nMono(index,1));

                nAdhes(index,1) = length(find(Adhesions.ActiveStatus));
                MemVel(index,1) = CalculateAveMembraneSpeed(Membrane,MembranePrevious,ModelParameters);
                if ShowPlot
                    [FAConnections,count] = PlotFilamentsAndMembrane(nth,count,Filaments,Membrane,Adhesions,FAConnections,FH,AH1,AH2,AH3,t,nMono,nAdhes,MemVel,index,TimeVec,ModelParameters);
%                   Write snapshots of figure as tif images
%                     if count == 0
%                     frame = frame + 1;
%                     F = getframe(1);
%                     imwrite(F.cdata,['X:\Mendoza Lab\MATLAB\Actin Growth Network Modeling\ModelImages\Frame_',sprintf('%04d',frame),'.tif'])
%                     end
                end

                SimData.Nodes(:,:,index) = single(Membrane.Nodes); % Record membrane segment positions
                count = count + 1;
                MembranePrevious = Membrane;
                
                drawnow
          
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

