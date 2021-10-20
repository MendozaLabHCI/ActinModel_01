clear; 
clc;
% Update SaveDir to a directory on your computer
SaveDir = 'X:\Mendoza Lab\MATLAB\Actin Growth Network Modeling\Test_ModelRun';    

% Create all value combinations first so that I can feed them all into parallel processing loops --------------------------------------------------------------------------------------
% In this simulation I'm verying adhesion lifetime curves (Pslope2, PA, Base) (See Test_DoubleExp.m)    
    Pslope2 = [  1.040500,  0.922800,  0.86060 ];
    PA      = [ -0.009727, -0.007286, -0.00612 ];
    Base    = [ 2, 0.5, 0.25];
    Pnum    = [1,2,3];
    nRuns = 10;
    values = zeros(0,nRuns);
    
    for p = Pnum
        for r = 1:nRuns
            values = [values; [p,Base(p),Pslope2(p),PA(p),r]];
        end
    end
    
    % I set this up for parallel processing, but I commented that line out (line 24) and just have running one simulation
    % Notice I set the simulation time length as 120 seconds (line 26)
    for k = 1:1
    %parfor (k = 1:size(values,1),4)
        ModelParameters = InitializeModelParameters;   % Create default values
        ModelParameters.TotalSimulationTime = 120; 
        ModelParameters.DoubleExp_Base   = values(k,2);
        ModelParameters.DoubleExp_Slope2 = values(k,3);
        ModelParameters.DoubleExp_A      = values(k,4);
        ModelParameters.Adhesion_OFF_DependentOnAdhesionFilamentTension = true;
        ModelParameters.AllowAdhesionFilamentBondBreak = true;
        ModelParameters.BoundaryFixed = true;
        ModelParameters.MembraneSpringConstant = 0.3; %3.0;
        ModelParameters.Adhesion_ON_ActivationRate = 1;
        ModelParameters.AdhesionTotalStartNumber = 2000;
        SaveName = ['Simulation-001__','Peak', num2str(values(k,1)), '_Run', sprintf('%02.0f',values(k,5)), '.mat'];
        if ~isfile(fullfile(SaveDir,SaveName))  % If code is re-run it will only run values not simulated yet
            disp(SaveName)
            [~] = SimulateLamellipodiumModel_KRC01(ModelParameters,SaveDir,SaveName,false);
        end
    end



