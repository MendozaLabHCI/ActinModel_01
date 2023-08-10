% Sample BR-MC model runs: 5 runs of Tau_peak = 3,7.5,12 [1,2,3]
clear; 
clc;
SaveDir = 'C:\DATA\Simulation_Output';   % Set a save location for Data files
    
   % Simulation setup -------------------------------------------------------------------------------------
    values = [];
    Runs = 1:5; 
    tau_peak  = [1,2,3]; 

    for n = 1:length(tau_peak)
            for r = 1:length(Runs)
                values = [values; [tau_peak(n),Runs(r)] ];
            end
    end

% Uncomment the following three lines for parallel processing -----------------------------------   
%    nProcessors = 4;
%    PoolObj = parpool('local',nProcessors);  
%    parfor (k = 1:size(values,1),nProcessors)
%------------------------------------------------------------------------------------------------

    for k = 1:size(values,1) % comment this line if parallel processing is used

               ModelParameters = InitializeModelParameters;
               % Filament and adhesion density fine adjustment Scaling Factor ---------------------------
                    Membrane = InitializeMembrane(ModelParameters);
                        width = (Membrane.Nodes(end,1)- Membrane.Nodes(1,1))/1000;
                        area1 = width*ModelParameters.AdhesionRegionDepth/1000;
                %-----------------------------------------------------------------------------------------
                 ModelParameters.AdhesionTotalStartNumber = round(area1*2000);
            %==========================================================================
                ModelParameters.TotalSimulationTime = 0.1;
                ModelParameters.TimeStep = 0.0002;
                ModelParameters.k_branch = 1.0;
                ModelParameters.branchWindowSize = 50; % nm
                ModelParameters.MembraneSpringConstant = 0.3;
                ModelParameters.CytoplasmViscosity = 10; % Pa*s
                % -------------------------------------------------------------------------------------
                ModelParameters.FilamentMassThreshold  = round(0.5*(15/2.7)*4500*width); % 
                ModelParameters.FilamentsInitialLength = round(1*(15/2.7)*6);
                ModelParameters.k_off_pointed =  12; % 12(s^-1)
                ModelParameters.k_on_barbed   = 165; % 13 s(-1) 
                ModelParameters.k_cap         =   3; % (s^-1)
                %--------------------------------------------------------------------------------------
                ModelParameters.MembraneViscosity = 100;
                NUmem = ModelParameters.MembraneViscosity*(10^-6); % Convert Pa*s to pN*s/nm^2
                ModelParameters.MembraneGamma = 6*pi*NUmem*ModelParameters.SegmentWidth; % 100 Pa*s.... previous value: 0.05;

                ModelParameters.Adhesion_MolecularClutchOn = true;
                ModelParameters.MolecularClutch_PeakNumber = values(k,1);
                ModelParameters.Adhesion_ActivationRate    = 1.0; % events/sec
                ModelParameters.Adhesion_DeActivationRate  = 0.1; % events/sec
                
                ModelParameters.BrownianRatchetOn = true;
                ModelParameters.BoundaryFixed     = false;
            %==========================================================================
            SaveName = ['SIMULATION__','TauPeak_',sprintf('%02.0f',values(k,1)),'_run_',sprintf('%02.0f',values(k,2)), '.mat'];
        
        if ~isfile(fullfile(SaveDir,SaveName))  % If code is re-run it will only run values for files not created/completedsimulation yet
                disp(SaveName);
                [~] = SimulateLamellipodiumModel_KRC01(ModelParameters,SaveDir,SaveName,false,false,false); % Run simulation
        end
   end
    disp('SIMULATION complete') 
    
    