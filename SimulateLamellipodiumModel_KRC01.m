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
    if RecordLargeData
        DATA = cell(length(TimeVec),1);
    end
    %PCD  = cell(length(TimeVec),1);
    MembranePrevious = [];
    index = 0;  
    count = 0;
    nth = 1; % nth time frame to plot
    
    
    % Get adhesions up to "stable" level before starting model ----------------------------------------------------------------------------------------------------------
    for k = 1:10000 %round(1/ModelParameters.TimeStep)
        [Adhesions,FAConnections] = ManageAdhesions(Filaments,Adhesions,FAConnections,Membrane,ModelParameters);
    end
    
    tic
    for t = TimeVec
        
        index = index + 1;
        
        [Filaments, Adhesions, FAConnections, PolymCoeffData] = PolymerizeDepolymerizeCapDeleteFilaments(Filaments,Adhesions,Membrane,FAConnections,ModelParameters);
        Filaments = BranchFilamentsInBranchWindowIfSelected(Filaments,ModelParameters,Membrane);
        [Adhesions,FAConnections] = ManageAdhesions(Filaments,Adhesions,FAConnections,Membrane,ModelParameters);
        [FAConnections, Adhesions] = FilamentAndAdhesionConnect(FAConnections,Filaments,Adhesions,ModelParameters);
        [Filaments, Membrane, Adhesions, FAConnections, Tensions, kBreaks, Data] = CalculatePositionAfterAppliedForces(Filaments,Membrane,Adhesions,FAConnections,ModelParameters);
        
        % Recording options 
        if RecordLargeData &&  RecordPCDdata
            Data.PolymCoeff = PolymCoeffData ;
            DATA{index,1} = Data;
        elseif ~RecordLargeData &&  RecordPCDdata
            Data = [];
            Data.PolymCoeff = PolymCoeffData ;
            DATA{index,1} = Data;
        elseif  RecordLargeData && ~RecordPCDdata 
            DATA{index,1} = Data;
        end
        
        nMono(index,1)  = CountTotalMonomers(Filaments);
        nAdhes(index,1) = length(find(Adhesions.ActiveStatus));
        MemVel(index,1) = CalculateAveMembraneSpeed(Membrane,MembranePrevious,ModelParameters);
        
        if ShowPlot
            [FAConnections,count] = PlotFilamentsAndMembrane(nth,count,Filaments,Membrane,Adhesions,FAConnections,FH,AH1,AH2,AH3,t,nMono,nAdhes,MemVel,index,TimeVec,ModelParameters);
        end
        
        [Filaments] = AddRandomFilaments(Filaments,Membrane,ModelParameters,nMono(index,1));
        % Record Small Data by default (membrane position)
        SimData.Nodes(:,:,index) = single(Membrane.Nodes); % Record membrane segment positions
        count = count + 1;
        MembranePrevious = Membrane;
        
    end
  
    SimData.nMonomers = nMono;
    SimData.nAdhesions = nAdhes;
    SimData.MembraneVelocity = MemVel;
    if RecordLargeData || RecordPCDdata
        SimData.DATA = DATA;
    end
    
    try
        save(fullfile(SaveDir,SaveName),'SimData','-mat','-v7.3')
    catch
        disp(['File not saved: ',SaveName])
    end
    toc
end

