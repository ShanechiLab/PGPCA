classdef TorusFitting < matlab.mixin.Copyable
%{
This class creates a "torus" object with methods to computes its "features" in one place.

(A) Model:  (R > r)

    [outAng, inAng] -> [(R + r*cos(inAng))*cos(outAng), (R + r*cos(inAng))*sin(outAng), r*sin(inAng)].

(B) Features:
    
    1.  Position:           Map 2D torus coordinate to 3D embedding space.
    2.  Tangent & normVec:  The local tangent space & normal vectors of 2D coordinate.
    3.  local area:         The local area factor computed by "coarea formula" and Jacobian matrix.
-----------------------------------------------------------------------------------------------------------
    
NOTE:
1.  This class is simple to speed up the PGPCA manifold-based simulation. Therefore, the model is really
      simplified. Further properties/methods need to be added for more general simulation.
2.  The name "TorusFitting" follows the 1D class "CurveFitting", even though the torus cannot be fitted
      by given samples right now (maybe a future method). Still, this naming is for consistency.
3.  For speeding up, the sample index = the last index when array dim > 2 due to MATLAB slicing method.
-----------------------------------------------------------------------------------------------------------

Properties:
--------------- Basic properties (to define a torus).
outR:           a +scalar. The outer radius of the torus.
inR:            a +scalar. The inner radius of the torus. Note that "outR > inR".
%}

%% Class main body.
    properties(SetAccess = protected, GetAccess = public)
        % Basic.
        outR;
        inR;
        dim = 3;
    end
    
    methods (Static)
        %% Overload "loadobj": since some variables shouldn't be initialized as an empty array.
        function obj = loadobj(obj)
            obj.RefreshProp(false);
        end
    end
    
    methods 
        %% Methods that defined in external m.file.
        %{
        This section lists templates of external functions. their I/O formats may not match the extermal 
          files 100%, but they must be consistent (e.g., the number of I/O).
        
        NOTE:
        1.  The 1st output of any functions must be the main output. This is needed for anonymous func.
        %}
        
        % ---------------------------------- Basic torus operator ----------------------------------
        % The torus embedding position.
        varargout = Position(obj, seleAng, varargin);
        % The local 1st order vectors (tangent & normal vector).
        varargout = LocalSpace(obj, spaceType, seleAng, varargin);
        % The local area factor based on Jacobian matrix.
        varargout = LocalArea(obj, seleAng, varargin);
        
        % ---------------------------------- Advance operator ----------------------------------
        % The various types of samples on the torus.
        varargout = GeneSample(obj, sampType, sampPara, varargin);
        % Three lattices "XX, YY, & ZZ" for plot the torus by "surf".
        varargout = GeneMeshgrid(obj, gridType, gridPara, varargin);
        
        % ---------------------------------- Analysis/Checking tool ----------------------------------
        % This function plots the manifold with selected features for visualization.
        varargout = PlotManifold(obj, varargin);  
        
        %% A initializer for properties (since there are "struct" properties).
        function obj = TorusFitting(outR, inR)
            % Assert: outR > inR, so the torus is well-defined.
            assert( outR > inR, 'The "outR" (%2.2e) must > "inR" (%2.2e).', outR, inR );
            % Assign properties.
            obj.outR = outR;
            obj.inR = inR;
            % Initialize the rest properties (if any).
            obj.RefreshProp(true);
        end
        
        %% This function checks all preperties. In short, it keeps different versions of this class consistent.
        function RefreshProp(~, hardRefreshSW)
            %{
            "hardRefreshSW" determines if all properties are initialized no matter they exist already or 
              not (default: false). This function is for back-up now.
            
            The 1st input should be "obj". For now it's "~" because it's unused.
            %}
            
            assert( islogical(hardRefreshSW), ...
                    'Input "hardRefresh" (class: %s) must be boolean', class(hardRefreshSW) );
        end
    end
end

