function ModelParameters = InitializeModelParameters

%% Time Properties
  
        ModelParameters.TotalSimulationTime = 5*60; %30; % seconds
        ModelParameters.TimeStep = 0.0002; % seconds 

%% General Properties
        
        ModelParameters.kT = 4.28; % k_B * T = 4.28 pN*nm (for T = 310 K)
        ModelParameters.CytoplasmViscosity = 10; % (Pa*s) https://summerschool.tugraz.at/images/phocadownload/Wirtz-Annu_Rev_Biophys-2009.pdf
        ModelParameters.MembraneViscosity = 100; % (Pa*s)

%% Adhesion Properties

        ModelParameters.AdhesionTotalStartNumber = 2000; % 800 MC = on, 400 MC = off; (for LeadingEdgeLength = 2000 nm and RegionsDepth = 200 nm)
        ModelParameters.AdhesionRegionDepth = 200;       % 500; % nm
        ModelParameters.Adhesion_ActivationRate = 1;     % events/sec
        ModelParameters.Adhesion_DeActivationRate = 0.1; % events/sec
        ModelParameters.Adhesion_MolecularClutchOn = true;
        ModelParameters.MolecularClutch_PeakNumber = 3;
        ModelParameters.AdhesionSpringConstant = 10; % 10 (picoNewtons/nanometer)
        ModelParameters.AdhesionSpringEqLength = 2;  % (nm) Adhesion spring equilibrium length
        
        %ModelParameters.AllowAdhesionFilamentBondBreak = true;
        % Peak 3 is default settings -----------------------------
%         ModelParameters.DoubleExp_Slope  =  6.8;
%         ModelParameters.DoubleExp_Base   =  0.25;
%         ModelParameters.DoubleExp_Slope2 =  0.8606;
%         ModelParameters.DoubleExp_A      = -0.00612;
%         ModelParameters.DoubleExp_B      =  0.3;
%         ModelParameters.DoubleExp_C      =  0.5E-5;
        
        
%% Membrane properties

        ModelParameters.BoundaryFixed = false;
        ModelParameters.BrownianRatchetOn = true;
        %ModelParameters.ContactThresholdForMembrane = 15; % (nm)
        ModelParameters.MemleadingEdgeLengthInNanometers = 2000; % 20000;
        ModelParameters.MemleadingEdgeHeightInNanometers = 0;
        ModelParameters.SegmentWidth = 18; % 18;   (nm)
        ModelParameters.SpringWidth = 2;   % 2; %10;  (nm)
        ModelParameters.MembraneSpringConstant = 0.3; % 0.03; pN/nm
            NUmem = ModelParameters.MembraneViscosity*(10^-6); % Convert Pa*s to pN*s/nm^2
        ModelParameters.MembraneGamma = 6*pi*NUmem*ModelParameters.SegmentWidth; % 100 Pa*s.

%% Filament properties
 
        ModelParameters.StartingNumberOfFilaments = 1;
        ModelParameters.FilamentsInitialLength = round(1*(15/2.7)*6); % Ave Number of monomers in an initialized filament
        ModelParameters.VerticalOffSet = -150; % (nm)
        ModelParameters.MonomerLength = 2.7; % 2.7; (nm)    
        ModelParameters.FilamentMassThreshold = round(0.5*(15/2.7)*4500*(ModelParameters.MemleadingEdgeLengthInNanometers/1000)); % (units of nMonomers) Threshold for adding new filaments to model
        ModelParameters.k_off_barbed = 0; % (s^-1)
        ModelParameters.k_off_pointed = 12; %6.5; %13; (s^-1)
        ModelParameters.k_on_barbed = 165; %12 s(-1)   http://cytomorpholab.com/wp-content/uploads/2021/11/publication_1-23.pdf (table 2)
        ModelParameters.k_cap = 3; % (s^-1)
        ModelParameters.k_branch = 1.0; %0.2  (s^-1)
        ModelParameters.branchWindowSize = 50; %15; % (nm)
        ModelParameters.branchAngle = 70;
        ModelParameters.branchAngleSTD = 10; % Standard deviation of branching angle in degrees
        ModelParameters.PersistenceLength = 1000; % (nm) Persistence length of actin filament (17.7 um)
        ModelParameters.VariableBranchRate = false;
        ModelParameters.VariableBranchRatePeriod = 120; % s
        ModelParameters.MinimumBranchSeparation = 50; % nm
        ModelParameters.MaxFilamentForceWithoutBR = 1; % pN
        ModelParameters.MCmodelContactThreshold = 2.7; % nm
        ModelParameters.FilamentDiameter = 3.5; % nm

        
        
        
        
           