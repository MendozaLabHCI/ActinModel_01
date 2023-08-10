function SimData = SimulateLamellipodiumModel_KRC01(ModelParameters,SaveDir,SaveName,ShowPlot,RecordLargeData,RecordPCDdata)

    % ModelParameters = InitializeModelParameters;
    % ShowPlot = true;
    Membrane  = InitializeMembrane(ModelParameters);
    Filaments = InitializeActinFilaments(ModelParameters,Membrane); % (#Filaments, InitialLength, SpreadWidth, ModelParameters)
    Adhesions = InitializeAdhesions(Membrane,ModelParameters);
    FAConnections = InitializeFAConnections;
    if ShowPlot  
        [FH,AH1,AH2,AH3] = SetUpPlottingFigureAndAxes; 
    end

    TimeVec = 0:ModelParameters.TimeStep:ModelParameters.TotalSimulationTime;
    nMono = NaN(length(TimeVec),1);
    nAdhes = NaN(length(TimeVec),1);
    MemVel = NaN(length(TimeVec),1);
    SimData = SetUpRecordMembrane(TimeVec,Membrane);
    SimData.ModelParameters = ModelParameters;
    
    DATA = cell(0,1);
    
    %PCD  = cell(length(TimeVec),1);
    MembranePrevious = [];
    index1 = 0;  
    index2 = 0;
    count = 0;
    nth = 1; % nth time frame to plot
    
    
    % Get adhesions up to "stable" level before starting model ----------------------------------------------------------------------------------------------------------
    for k = 1:10000 %round(1/ModelParameters.TimeStep)
        [Adhesions,FAConnections] = ManageAdhesions(Filaments,Adhesions,FAConnections,Membrane,ModelParameters);
    end
    
    tic
    
    for t = TimeVec
            index1 = index1 + 1;
            
            [Filaments, Adhesions, FAConnections, PolymCoeffData] = PolymerizeDepolymerizeCapDeleteFilaments(Filaments,Adhesions,Membrane,FAConnections,ModelParameters);
            Filaments = BranchFilamentsInBranchWindowIfSelected(t,Filaments,ModelParameters,Membrane);
            [Adhesions,FAConnections] = ManageAdhesions(Filaments,Adhesions,FAConnections,Membrane,ModelParameters);
            [FAConnections, Adhesions] = FilamentAndAdhesionConnect(FAConnections,Filaments,Adhesions,ModelParameters);
            [Filaments, Membrane, Adhesions, FAConnections, Tensions, kBreaks, Data] = CalculatePositionAfterAppliedForces(Filaments,Membrane,Adhesions,FAConnections,ModelParameters);
            
            Data.Time = t;
            PolymCoeffData.Time = t;

            % Recording options 
            if RecordLargeData &&  RecordPCDdata && rem(t,0.001) == 0
                Data.PolymCoeff = PolymCoeffData ;
                DATA = [DATA;{Data}];
            elseif ~RecordLargeData &&  RecordPCDdata && rem(t,0.001) == 0
                Data = [];
                Data.PolymCoeff = PolymCoeffData ;
                DATA = [DATA;{Data}];
            elseif  RecordLargeData && ~RecordPCDdata && rem(t,0.001) == 0
                DATA = [DATA;{Data}];
            end
            
            S = whos('DATA');
            if (S.bytes/10^9 >= 4 || t == TimeVec(end)) && (RecordLargeData || RecordPCDdata) 
                try
                    index2 = index2 + 1;
                    DATASaveName = [SaveName(1:end-4),'__DATA',sprintf('%03d',index2),SaveName(end-3:end)];
                    save(fullfile(SaveDir,DATASaveName),'DATA','-mat','-v7.3')
                    DATA = cell(0,1);
                catch
                    if index2 > 10
                        index2 = index2 - 1;
                    end
                    disp(['File not saved: ',DATASaveName])
                    pause(60)
                end
            end

            nMono(index1,1)  = CountTotalMonomers(Filaments);
            nAdhes(index1,1) = length(find(Adhesions.ActiveStatus));
            MemVel(index1,1) = CalculateAveMembraneSpeed(Membrane,MembranePrevious,ModelParameters);
            
            if ShowPlot
                [FAConnections,count] = PlotFilamentsAndMembrane(nth,count,Filaments,Membrane,Adhesions,FAConnections,FH,AH1,AH2,AH3,t,nMono,nAdhes,MemVel,index1,TimeVec,ModelParameters);
            end
            
            [Filaments] = AddRandomFilaments(Filaments,Membrane,ModelParameters,nMono(index1,1));
            % Record Small Data by default (membrane position)
            SimData.Nodes(:,:,index1) = single(Membrane.Nodes); % Record membrane segment positions
            count = count + 1;
            MembranePrevious = Membrane;
        
    end
  
    SimData.nMonomers = nMono;
    SimData.nAdhesions = nAdhes;
    SimData.MembraneVelocity = MemVel;

%     if RecordLargeData || RecordPCDdata
%         % Remove empty cells in DATA
%         idx = cellfun(@isempty,DATA); 
%         DATA(idx) = [];
%         SimData.DATA = DATA;
%     end
    
TryWriting = true;
while TryWriting
    try
        save(fullfile(SaveDir,SaveName),'SimData','-mat','-v7.3')
        TryWriting = false;
    catch
        disp(['File not saved: ',SaveName])
    end
end
toc
