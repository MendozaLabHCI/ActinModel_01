function [slope1,base,slope2,A,B,C] = MolecularClutchPeakParameters(PeakNumber)

    % PeakNumber = 1,2,or 3

        switch PeakNumber
            case 1 % Peak 1 (Tau_max = 3 2)
                slope1  =  6.8;
                base   =  2.0;
                slope2 =  1.0405;
                A      = -0.009727;
                B      =  0.3;
                C      =  0.5E-5;
            case 2 % Peak 2 (Tau_max = 7.5 s)
                slope1  =  6.8;
                base   =  0.5; % 2
                slope2 =  0.9228;  % 0.98; 
                A      = -0.007286;   % -0.01471;
                B      =  0.3;
                C      =  0.5E-5;
            case 3 % Peak 3 (Tau_max = 12 s)
                slope1  =  6.8;
                base   =  0.25; % 0.25
                slope2 =  0.8606;   % 0.9438;
                A      = -0.00612;  % -0.01728;
                B      =  0.3;
                C      =  0.5E-5;
            case 4 % Peak 4 (Tau_max = 16.5 s)
                slope1  =  6.8;
                base   =  0.165;    % change
                slope2 =  0.82;  % change
                A      = -0.005603; % change
                B      =  0.3;
                C      =  0.5E-5;    
            case 5 % Peak 5 (Tau_max = 21 s)
                slope1  =  6.8;
                base   =  0.1073;    % change
                slope2 =  0.7805;  % change
                A      = -0.0046; % change
                B      =  0.3;
                C      =  0.5E-5;  
        end
   
end
