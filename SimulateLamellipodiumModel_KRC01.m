function SimData = SimulateLamellipodiumModel_KRC01(ModelParameters,SaveDir,SaveName,ShowPlot)

    % ModelParameters = InitializeModelParameters;
    % ShowPlot = true;
    Filaments = InitializeActinFilaments(ModelParameters); % (#Filaments, InitialLength, SpreadWidth, ModelParameters)
    Membrane  = InitializeMembrane(ModelParameters);
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
    DATA = cell(length(TimeVec),1);
    PCD  = cell(length(TimeVec),1);
    MembranePrevious = [];
    index = 0;  
    count = 0;
    nth = 100; % nth time frame to plot
    
    
    % Get adhesions up to "stable" level before starting model ----------------------------------------------------------------------------------------------------------
    for k = 1:round(30/ModelParameters.TimeStep)
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
        
        PCD{index,1} = PolymCoeffData;                                        
        DATA{index,1} = Data;
        nMono(index,1) = CountTotalMonomers(Filaments);
        [Filaments] = AddRandomFilaments(Filaments,Membrane,ModelParameters,nMono(index,1));
        
        nAdhes(index,1) = length(find(Adhesions.ActiveStatus));
        MemVel(index,1) = CalculateAveMembraneSpeed(Membrane,MembranePrevious,ModelParameters);
        if ShowPlot
            [FAConnections,count] = PlotFilamentsAndMembrane(nth,count,Filaments,Membrane,Adhesions,FAConnections,FH,AH1,AH2,AH3,t,nMono,nAdhes,MemVel,index,TimeVec,ModelParameters);
        end
        
        SimData.Nodes(:,:,index) = single(Membrane.Nodes); % Record membrane segment positions
        count = count + 1;
        MembranePrevious = Membrane;
        
    end
  
    SimData.nMonomers = nMono;
    SimData.nAdhesions = nAdhes;
    SimData.MembraneVelocity = MemVel;
    SimData.DATA = DATA;
    SimData.PCD = PCD;
    
    try
        save(fullfile(SaveDir,SaveName),'SimData','-mat','-v7.3')
    catch
        disp(['File not saved: ',SaveName])
    end
    toc
end

