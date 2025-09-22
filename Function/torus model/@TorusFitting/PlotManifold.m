function [axisH, objHColl, torusInfo, tangInfo, normDInfo, asInfo] = PlotManifold(obj, varargin)
%{
2024/01/30: (A) Add option "localSpaceType" for plotting the torus coordinate (EuCOV/GeCOV).
-----------------------------------------------------------------------------------------------------------

This function plots the torus surface with following features. It's for visualization basically.

(A) Tangent space:      Along the outer (blue) and inner (red) radius. (w/ or w/o normalization)
(B) norm-D space:       A vector (green) orthogonal to the tangent space. (w/ or w/o normalization) 
(C) Area scaling:       It demonstrates by the dot size, because attaching a circle or torus is tricky.

NOTE: 
1.  In principle, there are two kinds of points:
    (a) interp point:       The points for outlining the torus surface.
    (b) highlight point:    The points for markering specific features (e.g., space, area scale).
2.  Basically, all functional settings (surfSet, vecSet, asDotSet, etc.) are attaching based. That is, the
      user-given settings are attached at the end of default ones. This avoids removing all settings when
      you only want to change one default property. It saves time.
-----------------------------------------------------------------------------------------------------------

Inputs:
------------------- Optional ('tag' + 'value' with (number))
------------------- (ABOUT FIGURE)
"axisH":            (1) (axisH) an axis handle. (def: auto-generated axis handle)
                        The axis handle for plotting.
------------------- (ABOUT TORUS SURFACE)
"interpP":          (2) (interpType) a string "number / sample". (def: "number")
                        (interpPara) a parameter. (def: [50, 20])
                        They decide the types of samples for plotting the manifold. The options are:
    "lattice":              interpPara = a row [1, 2] = [numIntpOut, numIntpIn]. The #sample along radius.
    "given":                interpPara = a cell row {1,3} = {XX, YY, ZZ}. The 2D arrays for func. "surf".
"torusSW":          (1) (torusSW) a boolean (def: true).
                        "true" means the torus will be plotted.
"surfSet":          (1) (surfSet) a cell row. (def: {"FaceColor", [1,1,1]*0.5, "FaceAlpha", 0.4})
                        The setting of "surf" function for plotting the torus.
------------------- (ABOUT TANGENT/norm-D SPACE)
"highLP":           (2) (highLType) a string "number / sample". (def: "number")
                        (highLPara) a parameter. (def: [50, 20])
                        They are points for tangent/norm-D vectors and area scaling. The options are:
    "lattice":              highLPara = a row [1, 2] = [numHighOut, numHighIn]. The #sample along radius.
    "given":                highLPara = a 2D array [numHighP, 2]. Each [j,:] is a highlight point.
"LocalSpaceType":   (1) (localSpaceType) a string ""





"tangVecInd":       (1) (tangVecInd) a row [1, numPTang] / empty []. (def: [1,2]) 
                        The selected tangent vector index to be plotted. [] means no vector.
"normDVecInd":       (1) (normDVecInd) a row [1, numPnormD] / empty []. (def: [1])
                        The selected norm-D vector index to be plotted. [] means no vector.
"tangColor":        (1) (tangColor) a string row / a 2D array [numPTang, 3] (def: ["blue", "red"])
                        The color of each tangent vector. Note that when #color = 1, all vectors have the
                          same color.
"normDColor":        (1) (normDColor) a string row / a 2D array [numPnormD, 3] (def: ["LimeGreen"])
                        The color of each norm-D vector. Note that when #color = 1, all vectors have the
                          same color.
"vecNormSW":        (1) (vecNormSW) a boolean (def: true)
                        "true" means all vectors are normalized to length = 1.
"vecLength":        (1) (vecLength) a positive scalar. (def: 0.5*obj.inR if vecNormSW = true, otherwise 1)
                        The length of the tangent and the norm-D vectors.
"vecSet":           (1) (vecSet) a cell row. (def: {'LineWidth', 2, 'AutoScale', 'off'})
                        It sets the "quiver3" properties for plotting tangent/norm-D vectors.
------------------- (ABOUT AREA SCALING)
"asFactor":         (1) (asFactor) a +scalar. (def: 10)
                        The scaling factor for area dot: "dotSize = areaScale * asFactor".
"asPSW":            (1) (asPSW) a boolean. (def: true)
                        "true" means the 3D dots for "area scaling" is plotted.
"asDotSizeFixSW":   (1) (asDotSizeFixSW) a boolean. (def: false if "vecNormSW = false", otherwise true).
                        "true" means "dotSize = asFactor". No scaling based on local Jacobian matrix.
"asDotSet":         (1) (asDotSet) a row cell. (def: {rgb("orange"), 'filled', 'marker', 'o'})
                        The setting for "scatter3" to show "area scaling".
------------------- (FINISH)

Outputs:
------------------- Essential
axisH:              an axis handle. The axis for the plot.
objHColl:           A struct. It collects all plotting object handles in this function. Empty [] means 
                      that object is not plotted. The fields are:
    torusH:             a handle. The torus surface.
    tangVecH:           a handle row [1, numPTang]. Each is a handle for one tangent direction.
    normDVecH:           a handle row [1, numPnormD]. Each is a handle for one norm-D direction.
    areaDotH:           a handle. The 3D dots for area scaling.
torusInfo:          A struct. It contains below fields about the torus. 
    interpZ:            a 2D array [numInterp, 2]. The interpolating radian coordinate [outAng, inAng].
    interpY:            a 2D array [numInterp, obj.dim]. The embedding coordinate of "interpZ".
    XX:                 a 2D array. The X coordinate for the torus surface.
    YY:                 a 2D array. The Y coordinate for the torus surface.
    ZZ:                 a 2D array. The Z coordinate for the torus surface.
tangInfo:           A struct. It contains fields about the tangent vectors as below.
    highLZ:             a 2D array [numHighL, 2]. The highlighting radian coordinate [outAng, inAng].
    highLY:             a 2D array [numHighL, obj.dim]. The embedding coordinate of "highLZ".
    tangVec:            a 3D array [numPTang, obj.dim, numHighL]. [i,:,j] is a tangent at highLZ(j,:).
    vecLength:          a +scalar. The scaling factor of the plotted tangent vector.
normDInfo:           A struct. It contains fields about the norm-D as below.
    normDVec:            a 3D array [numPnormD, obj.dim, numHighL]. [i,:,j] is an norm-D at highLZ(j,:).
    vecLength:          a +scalar. The scaling factor of the plotted norm-D vector.
asInfo:             A struct. It contains fields about the area scaling as below.
    asDotSizeFixSW:     a boolean. It's equal to optional input "asDotSizeFixSW".
    areaScale:          a column [numHighL, 1]. The area scaling value of each highlight point.
%}

%% Constants.
% Extract object's values for convenience.
% % % outR   = obj.outR;
inR    = obj.inR;
objDim = obj.dim;

% The maximal norm-D dimension.
normDimLimit = objDim - 2;

%% Defaults.
% Figure: axis.
axisH = [];

% Torus: surface interpolation.
interpType      = "lattice";
interpPara      = [50, 20];
torusSW         = true;
surfSet         = {"FaceColor", [1,1,1]*0.7, "FaceAlpha", 0.4, 'edgeColor', rgb('LightGray')};

% Tangent space: highlight points.
highLType       = "lattice";
highLPara       = [10, 5];
tangVecInd      = [1, 2];
tangColor       = ["blue", "red"];

% norm-D space: other settings.
normDVecInd      = [1:1:objDim-2];
normDColor       = "LimeGreen";

% Vector space: shared properties.
vecNormSW       = true;
vecLength       = [];
vecSet          = {'LineWidth', 2, 'AutoScale', 'off'};

% Area scaling: 3D dot.
asFactor        = 40;
asPSW           = true;
asDotSizeFixSW  = [];
asDotSet        = {rgb("orange"), 'filled', 'marker', 'o'};

%% Assign inputs.
% Optional inputs (if they exist).
numVarargin = length(varargin);
if(numVarargin > 0)
    tagInd = 1;
    while(tagInd <= numVarargin)
        % Find the option.
        tagName = varargin{tagInd};
        switch tagName
            case 'axisH';           tagShift = 1;       axisH          = varargin{tagInd + 1};
            case 'interpP';         tagShift = 2;       interpType     = varargin{tagInd + 1};
                                                        interpPara     = varargin{tagInd + 2};
            case 'torusSW';         tagShift = 1;       torusSW        = varargin{tagInd + 1};
            case 'surfSet';         tagShift = 1;       surfSet        = [surfSet, varargin{tagInd + 1}];
            case 'highLP';          tagShift = 2;       highLType      = varargin{tagInd + 1};
                                                        highLPara      = varargin{tagInd + 2};
            case 'tangVecInd';      tagShift = 1;       tangVecInd     = varargin{tagInd + 1};    
            case 'tangColor';       tagShift = 1;       tangColor      = varargin{tagInd + 1}; 
            case 'normDVecInd';      tagShift = 1;       normDVecInd     = varargin{tagInd + 1};    
            case 'normDColor';       tagShift = 1;       normDColor      = varargin{tagInd + 1};    
            case 'vecNormSW';       tagShift = 1;       vecNormSW      = varargin{tagInd + 1};
            case 'vecLength';       tagShift = 1;       vecLength      = varargin{tagInd + 1};    
            case 'vecSet';          tagShift = 1;       vecSet         = [vecSet, varargin{tagInd + 1}];    
            case 'asFactor';        tagShift = 1;       asFactor       = varargin{tagInd + 1};    
            case 'asPSW';           tagShift = 1;       asPSW          = varargin{tagInd + 1};   
            case 'asDotSizeFixSW';  tagShift = 1;       asDotSizeFixSW = varargin{tagInd + 1};
            case 'asDotSet';        tagShift = 1;       asDotSet       = [asDotSet, varargin{tagInd + 1}];    
            otherwise
                error('The given option is "%s", which is invalid.', tagName);
        end
        % Shift the index.
        tagInd = tagInd + tagShift + 1;
    end
end

%% Standardize parameters & sanity check.
% (if not given) Set default "vecLength".
if( isempty(vecLength) )
    if(vecNormSW); vecLength = 0.7*inR; else; vecLength = 0.3*inR; end
end

% (if not given) Set default "asDotSizeFixSW".
if( isempty(asDotSizeFixSW) )
    asDotSizeFixSW = vecNormSW;
end

% Assert: #tangent & norm-D vectors are valid.
numPTang = length(tangVecInd);
numPnormD = length(normDVecInd);
assert( numPTang <= 2, 'The #tangent vector (%d) must be <= 2.', numPTang );
assert( numPnormD <= normDimLimit, 'The #norm-D vector (%d) must be <= %d.', numPnormD, normDimLimit );

% Standardize all color options to match #color = #numPTang/numPnormD.
tangColor = ColorSetTransform(tangColor, numPTang);
normDColor = ColorSetTransform(normDColor, numPnormD);

%% Initialize variables for the outputs.
objHColl  = struct('torusH', [], 'tangVecH', [], 'normDVecH', [], 'areaDotH', []);
torusInfo = struct('interpZ', [], 'interpY', [], 'XX', [], 'YY', [], 'ZZ', []);
tangInfo  = struct('highLZ', [], 'highLY', [], 'tangVec', [], 'vecLength', vecLength);
normDInfo  = struct('normDVec', [], 'vecLength', vecLength);  
asInfo    = struct('areaScale', []);

%% Compute: interpolating points. (output: XX, YY, & ZZ)
% Select the method.
switch interpType
    case "lattice"
        % CASE: create a 2D uniform lattice for the torus.
        [XX, YY, ZZ, resColl] = obj.GeneMeshgrid("numGrid", interpPara);
        
        % Collect results.
        torusInfo.interpZ = resColl.gridZ;
        torusInfo.interpY = resColl.gridY;
        
    case "given"
        % CASE: surface arrays are given already. Just collect them.
        [XX, YY, ZZ] = deal(interpPara{:});
        
    otherwise
        error('The given "interpType" (%s) is invalid.', interpType);
end

% Collect results.
torusInfo.XX = XX;
torusInfo.YY = YY;
torusInfo.ZZ = ZZ;

%% Compute: highlight points. (output: highLZ & highLY)
% Select the method.
switch highLType
    case "lattice"
        % CASE: Generate uniform lattice points.
        [highLY, highLZ] = obj.GeneSample("lattice", highLPara);
        
    case "given"
        % CASE: The "highLZ" is given. Just collect it.
        highLZ = highLPara;
        highLY = obj.Position(highLZ);
        
    otherwise
        error('The given "highLType" (%s) is invalid.', highLType);
end

% Collect results.
tangInfo.highLZ = highLZ;
tangInfo.highLY = highLY;

%% Build the axis object.
if(isempty(axisH))
    figH = figure('position', [100,50,900,900], 'Color', 'white');
    axisH = subplot(1,1,1,'parent',figH);
    % Set the view.
    set(axisH, 'view', [33,36]);
end

% Set up the environment.
hold(axisH, 'on');
axis(axisH, 'equal');
box(axisH, 'off');

% The grid.
grid(axisH, 'on');

% Add light to view 3D object clearer.
camlight(axisH);
lighting(axisH, 'gouraud');

%% Plot: torus.
if( torusSW )
    objHColl.torusH = surf(axisH, XX, YY, ZZ, surfSet{:});
end

%% Plot: tangent space.
if( numPTang > 0 )
    % Compute tangent vectors.
    tangVec = obj.LocalSpace("tangent", highLZ, "normSW", vecNormSW);
    
    % Only plot the selected tangent vector.
    tangVec = tangVec(tangVecInd, :, :);
    
    % Plot each tangent vector one by one.
    tangVecPlot = permute(tangVec, [3,2,1]) * vecLength;
    objHColl.tangVecH = gobjects(1, numPTang);
    for(j = 1:1:numPTang)
        objHColl.tangVecH(j) = ...
            quiver3(axisH, highLY(:,1)       , highLY(:,2)       , highLY(:,3)       , ...
                           tangVecPlot(:,1,j), tangVecPlot(:,2,j), tangVecPlot(:,3,j), ...
                           'color', tangColor(j,:), vecSet{:});
    end
    
    % Collect results.
    tangInfo.tangVec = tangVec;
end

%% Plot: norm-D space.
if( numPnormD > 0 )
    % Compute tangent vectors.
    normDVec = obj.LocalSpace("norm-D", highLZ, "normSW", vecNormSW);
    
    % Only plot the selected orthonormal vector.
    normDVec = normDVec(normDVecInd, :, :);
    
    % Plot each tangent vector one by one.
    normDVecPlot = permute(normDVec, [3,2,1]) * vecLength;
    objHColl.normDVecH = gobjects(1, numPnormD);
    for(j = 1:1:numPnormD)
        objHColl.normDVecH(j) = ...
            quiver3(axisH, highLY(:,1)       , highLY(:,2)       , highLY(:,3)       , ...
                           normDVecPlot(:,1,j), normDVecPlot(:,2,j), normDVecPlot(:,3,j), ...
                           'color', normDColor(j,:), vecSet{:});
    end
    
    % Collect results.
    normDInfo.normDVec = normDVec;
end

%% Plot: area scaling.
if( asPSW )
    % Compute local area factor.
    areaScale = obj.LocalArea(highLZ, "compMd", "coarea");
    
    % Plot area by "scatter3".
    if( asDotSizeFixSW )
        asDotSize = asFactor;
    else
        asDotSize = areaScale * asFactor;
    end
    objHColl.areaDotH = scatter3(axisH, highLY(:,1), highLY(:,2), highLY(:,3), asDotSize, asDotSet{:});
    
    % Collect results.
    asInfo.asDotSizeFixSW = asDotSizeFixSW; 
    asInfo.areaScale = areaScale;

end

end

function colorMat = ColorSetTransform(colorMat, numVec)
%{
This AUXILIARY function transforms "colorMat" into its standard format. The class(colorMat) can be:

1.  char/string:    ["red", "blue"]  
2.  double:         [1,0,0; 0,1,0]

We transform the class into "double" to simplify further manipulation with dimension check.
%}

%% Process.
% Transformation based on the input class.
classColorMat = class(colorMat);
switch classColorMat
    case 'char' 
        % CASE: only a single color.
        colorMat = repmat(rgb(colorMat), numVec, 1);
        
    case 'string'
        % CASE: maybe multiple colors.
        assert( isrow(colorMat), '"colorMat" must be a row for transformation.' );
        numColorMat = length(colorMat);
        if( numColorMat == 1 )
            colorMat = repmat(rgb(colorMat), numVec, 1);
        else
            colorMat = cell2mat(arrayfun( @rgb, colorMat', 'UniformOutput', false ));
        end
            
    case 'double'
        % CASE: The format is clean already. Do nothing.

    otherwise
        error('The "colorMat" class (%s) is invalid.', classColorMat);
end

% Only use the #numVec colors.
colorMat = colorMat(1:numVec, :);

end

