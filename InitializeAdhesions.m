function Adhesions = InitializeAdhesions(Membrane,ModelParameters)
    
    SpringWidth = ModelParameters.SpringWidth;
    MembraneLength = Membrane.Nodes(end,1) - Membrane.Nodes(1,1);
    AdhesionRegionArea = MembraneLength*ModelParameters.AdhesionRegionDepth;
    nA = ceil((ModelParameters.AdhesionTotalStartNumber/(1000^2))*AdhesionRegionArea);
    nS = size(Membrane.Segments,1); 
    Adhesions.RegionNodes = zeros(4,2,nS); 
    
            for r = 1:nS % For each membrane segment create rectanglular region behind it, and create adhesions with a fixed density
               XY_LT = [ Membrane.Nodes(Membrane.Segments(r,1),1) - SpringWidth/2,... %  Left-Top xy-coordinates of current adhesion region
                         Membrane.Nodes(Membrane.Segments(r,1),2)];        
               XY_RT = [ Membrane.Nodes(Membrane.Segments(r,2),1) + SpringWidth/2,...  %  Right-Top xy-coordinates of current adhesion region
                         Membrane.Nodes(Membrane.Segments(r,2),2)]; 
               XY_RB = [ Membrane.Nodes(Membrane.Segments(r,2),1) + SpringWidth/2,...  %  Right-Bottom xy-coordinates of current adhesion region
                         Membrane.Nodes(Membrane.Segments(r,2),2) - ModelParameters.AdhesionRegionDepth];
               XY_LB = [ Membrane.Nodes(Membrane.Segments(r,1),1) - SpringWidth/2,... %  Left-Bottom xy-coordinates of current adhesion region
                         Membrane.Nodes(Membrane.Segments(r,1),2) - ModelParameters.AdhesionRegionDepth];   
                          
               % Plotting these points should create a vertical rectangle
%                plot([XY_LT(1); XY_RT(1); XY_RB(1); XY_LB(1); XY_LT(1)],...
%                     [XY_LT(2); XY_RT(2); XY_RB(2); XY_LB(2); XY_LT(2)],'.-')
%                hold on  

               % Record Current Region
               Adhesions.RegionNodes(1,:,r) = XY_LT;           
               Adhesions.RegionNodes(2,:,r) = XY_RT;            
               Adhesions.RegionNodes(3,:,r) = XY_RB; 
               Adhesions.RegionNodes(4,:,r) = XY_LB; 
               
%                delX = Adhesions.RegionNodes(2,1,r) - Adhesions.RegionNodes(1,1,r);
%                delY = Adhesions.RegionNodes(1,2,r) - Adhesions.RegionNodes(4,2,r);
            end
            
            Adhesions.XYPoints = NaN(nA,2);
            Adhesions.RegionLocation = NaN(nA,1);
            Adhesions.ActiveStatus   = false(nA,1);
            Adhesions.AttachedFilamentName = NaN(nA,1);
            
end

