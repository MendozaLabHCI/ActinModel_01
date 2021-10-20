function ModelParameters = InitializeModelParameters

%% Time Properties
  
        ModelParameters.TotalSimulationTime = 10; %30; % seconds
        ModelParameters.TimeStep = 0.001; % seconds 
        
%% Adhesion Properties

        ModelParameters.AdhesionTotalStartNumber = 1000; % 1000 adhesions/um^2 or 1e-4 adhesions/nm^2
        ModelParameters.AdhesionRegionDepth = 200; % 500; % nm
        ModelParameters.Adhesion_ON_ActivationRate = 0.1; % events/sec
        ModelParameters.Adhesion_OFF_ActivationRate = 0.1; % events/sec
        ModelParameters.Adhesion_OFF_DependentOnAdhesionFilamentTension = false;
        ModelParameters.AdhesionSpringConstant = 10; % 10 (picoNewtons/nanometer)
        ModelParameters.AdhesionSpringEqLength = 2; % (nm) Adhesion spring equilibrium length
        
        ModelParameters.AllowAdhesionFilamentBondBreak = false;
        % Peak 3 is default settings
        ModelParameters.DoubleExp_Slope  =  6.8;
        ModelParameters.DoubleExp_Base   =  0.25;
        ModelParameters.DoubleExp_Slope2 =  0.8606;
        ModelParameters.DoubleExp_A      = -0.00612;
        ModelParameters.DoubleExp_B      =  0.3;
        ModelParameters.DoubleExp_C      =  0.5E-5;
        
        
%% Membrane properties
    
        ModelParameters.BoundaryFixed = false;
        ModelParameters.ContactThresholdForMembrane = 15; % (nm)
        ModelParameters.MemleadingEdgeLengthInNanometers = 2000; % 20000;
        ModelParameters.MemleadingEdgeHeightInNanometers = 0;
        ModelParameters.SegmentWidth = 18;% 18; %90;   (nm)
        ModelParameters.SpringWidth = 2;% 2; %10;  (nm)
        ModelParameters.MembraneSpringConstant = 0.3; % 0.03; pN/nm
        ModelParameters.BoundaryForceSpringConstant = 10;% 0.1; pN/nm
        ModelParameters.MembraneGamma = 0.1; % (pN*S/nm)
        

%% Filament properties

        ModelParameters.StartingNumberOfFilaments = 1;
        ModelParameters.FilamentsInitialLength = 15; % Ave Number of monomers in an initialized filament
        ModelParameters.SpreadWidth = 1800; %1400; % (nm) % Range in the x-direction where filaments are randomly created
        ModelParameters.VerticalOffSet = -100; % (nm)
        ModelParameters.MonomerLength = 15; % 2.7; (nm)            
        ModelParameters.FilamentForce = 1; % pN (vector force along directin of filament)
        ModelParameters.FilamentMassThreshold = 9000; % (units of nMonomers) Threshold for adding new filaments to model
        ModelParameters.k_off_barbed = 0; % (s^-1)
        ModelParameters.k_off_pointed = 12; %6.5; %13; (s^-1)
        ModelParameters.k_on_barbed = 13; %12 s(-1)  % https://www.ncbi.nlm.nih.gov/books/NBK9908/    (Figure 11.3)  
        ModelParameters.k_cap = 0.6; % (s^-1)
        ModelParameters.k_branch = 0.5; %0.2  (s^-1)
        ModelParameters.branchWindowSize = 15; % (nm)
        ModelParameters.branchAngle = 70;
        ModelParameters.branchAngleSTD = 10; % Standard deviation of branching angle in degrees
    
           