% These sample BR model simulations vary r_on from 0.1 to 0.5 (0.1 steps) with a total of 10 runs each
clear; 
clc;
SaveDir = 'D:\KEITH\Simulation_10secRuns_007b';    
    
   %Simulation setup -------------------------------------------------------------------------------------
    values = [];
    nRuns = 10;
    r_on = 0.1:0.1:0.5;
    %tau_peak  = [1,2,3];
   
    for m = 1:length(r_on)
        for r = 1:nRuns
            values = [values; [r_on(m),r] ];
        end
    end
    
    nProcessors = 4;
    parfor (k = 1:size(values,1),nProcessors)
        
             ModelParameters = InitializeModelParameters;
            % Edit default ModelParameterrs==========================================================================
                ModelParameters.TotalSimulationTime = 10;
                ModelParameters.TimeStep = 0.001;
                
                ModelParameters.MemleadingEdgeLengthInNanometers = 2000; %2000
                ModelParameters.SegmentWidth = 18;
                SW = ModelParameters.SegmentWidth;
                ModelParameters.MembraneGamma = 4.3092e-7*SW^2 + 2.4082e-3*SW + 6.0469e-3; % 0.05
                ModelParameters.MembraneSpringConstant = 0.3;
                 
                ModelParameters.k_branch = 0.5;
                
                ModelParameters.Adhesion_MolecularClutchOn = false;
                ModelParameters.MolecularClutch_PeakNumber = 1;
                ModelParameters.Adhesion_ActivationRate   = values(k,1); % events/sec
                ModelParameters.Adhesion_DeActivationRate = 0.3; % events/sec
                
                % Filament and adhesion density fine adjustment Scaling Factor ---------------------------
                    Membrane = InitializeMembrane(ModelParameters);
                        width = (Membrane.Nodes(end,1)- Membrane.Nodes(1,1))/1000;
                        area1 = width*ModelParameters.AdhesionRegionDepth/1000;
                %-----------------------------------------------------------------------------------------
                % Set maximum allowed filaments
                ModelParameters.FilamentMassThreshold = round(4500*width); %9000
                % Set maximum allowed adhesions
                if ModelParameters.Adhesion_MolecularClutchOn
                    ModelParameters.AdhesionTotalStartNumber = round(area1*2000); % 800 MC = on, 400 MC = off; (for LeadingEdgeLength = 2000nm)
                else
                    ModelParameters.AdhesionTotalStartNumber = round(area1*1000);
                end
                
                ModelParameters.BrownianRatchetOn = true;
                ModelParameters.BoundaryFixed     = false;
            %==========================================================================
            SaveName = ['SIMULATION-007b__','Ron_',SimFormat(values(k,1)), '_run_',sprintf('%02.0f',values(k,2)), '.mat'];
        
        if ~isfile(fullfile(SaveDir,SaveName))  % If code is re-run it will only run values for files not created/completedsimulation yet
            disp(SaveName);
            [~] = SimulateLamellipodiumModel_KRC01(ModelParameters,SaveDir,SaveName,false,false,false); % Run simulation
        end
    end
    disp('SIMULATION-007b complete') 
    
    