function [axisH, objHColl, probInfo] = ...
    PlotProbDensity(gridCell, probYType, sampZ, paraStru, funStru, varargin)
%{
This function plots a 2D/3D probability density figure from the given PGPCA model. 

Model:  y = g(z) + K(z)*C*x + r (COV: R = (sigma^2)*I_n)

where:  z ~ manifold coordinate with a probability distribution.
        x ~ N(0, I_m), a multivariate normal distribution.

The types of the figure:

1.  2D:             The embedding space = R^2. The figure is a "surf" whose "Z" is density.
2.  3D:             The embedding space = R^3. The figure includes 3 slice planes "XY, YZ, & XZ".
-----------------------------------------------------------------------------------------------------------

The steps of this function:

1.  (if asked) prepare "funProbY".
2.  Compute p(y) over the given mesh grids.
3.  Plot surfaces following the case of "dim(Y)".
4.  Organize the figure (e.g., scatter dot, colorbar).
5.  Prepare the outputs.

The types of "probYType":

1.  "given":        [numY, 1] = funProbY(y); Key function "funProbY" is given directly.
2.  "parameter":    We need to build "p(y) = sum_z(p(y|z) * p(z))" from the given parameters below:

    (A) [numZ, 1]    = funProbZ(z);                     p(z) on the manifold as landmarks.
    (B) [numY, numZ] = funProbYlZ(y, z, paraStru);      p(y|z) with "paraStru", an additional parameter.
                  or = funProbYlZ(y, z);                p(y|z) with all settings in the function handle.

NOTE:
1.  To speed up the implementation, we don't include "space projection" for high-D embedding, so g(z) must
      be 2D/3D. Also, the slices are "elementary planes" without any shift/rotation.
2.  The 3D slices are actually from 3 "surf" functions, not MATLAB "slice" function. This saves the memory
      and computation loading a lot! Otherwise, a 32GB RAM can blow out sometimes!
3.  We include multiple "probYType" to handle various representations of the PGPCA model.
4.  Variables "cb*" are colorbar-related. We include them since colorbar is necessary for probability.
5.  Scatter dash lines are added in 3D to see the intersection between slices.
-----------------------------------------------------------------------------------------------------------

Inputs:
------------------- Essential
------------------- (BASIC ELEMENT)
gridCell:           a cell row {1,2/3} = {gridX, gridY, gridZ}. Every "grid*" is a vector of axis's grid.
probYType:          a string "given / parameter". The type of computing "p(y)". Details are above.
------------------- (ABOUT PGPCA PARAMETERS)
sampZ:              a 2D array [numZ, dimZ]. Each [j,:] is a landmark on the manifold g(z).
paraStru:           a struct. It includes learning parameters in following fields:
    matC:               a 2D array [dimL, dimC] / empty []. The dimensionality reduction matrix.
    mainVar:            a scalar >= 0. The average noise VAR of the "modeled space".
    sideVar:            a scalar >= 0 / empty []. The average noise VAR of the "unmodeled space".
funStru:            a struct. It includes all functional handles in following fields:
    funProbYlZ:         a func. handle. Probability p(y|z). Details are above.
    funProbZ:           a func. handle. Probability p(z). Details are above.
    funProbY:           a func. handle. Proabbility p(y). Details are above.
------------------- (FINISH)
------------------- Optional ('tag' + 'value' with (number))
------------------- (ABOUT AXIS)
"axisH":            (1) (axisH) an axis handle. (def: auto-generated axis handle)
                        The axis handle for plotting.
"axisSet":          (1) (axisHSet) a cell row. (def: {'xlim', gridX([1,end]), 'ylim', ..., 'zlim', ...})
                        The axis setting used in command "set(axisH, axisSet{:})".
"tickGap":          (1) (tickGap) a +scalar. (def: [])
                        The gap between ticks in X & Y axes. Z axis is included if it's 3D.
"grid":             (1) (gridSW) a boolean. (def: "true" if it's 3D. Otherwise "false".)
                        "true" means the grid is added to the figure.
"view":             (1) (viewAngle) a row [1, 2]. (def: [-33, 36] for 3D & [0, 90] for 2D)
                        The viewing angle of ths figure.
------------------- (ABOUT SCATTER DASH LINE)
"scatSW":           (1) (scatSW) a boolean. (def: "true" if it's 3D. Otherwise "false")
                        "true" means the dashed lines by "scatter3/scatter" are added.
"scatNumDot":       (1) (scatNumDot) a +scalar / row [1, 2/3]. (def: 25)
                        The number of dots in each dashed line along 2/3 axes. A scalar means all axes use
                          the same value, and a row vector is for each axis.
"scatSet":          (1) (scatSet) a cell row. (def: {'SizeData',        10, ...
                                                     'MarkerFaceColor', rgb('Darkgray'), ...
                                                     'MarkerEdgeColor', 'none'})
                        The setting for function "scatter3/scatter".
------------------- (ABOUT COLORBAR)
"cbMap":            (1) (cbMap) a string / a 2D array [numColor, 3]. (def: "turbo")
                        The color set of colorbar.
"cbSet":            (1) (cbSet) a cell row. (def: {})
                        The setting of colorbar properties.
"cbLimit":          (1) (cbLimit) a row [1,2] = [cbMin, cbMax]. (def: [])
                        The limit of colorbar. It's used by function "caxis".
"cbLabel":          (1) (cbLabel) a string. (def: "Probability")
                        The label of colorbar.
------------------- (FINISH)

Outputs:
------------------- Essential
axisH:              an axis handle. The axis for the plot.
objHColl:           a struct. It collects all plotting object handles in following fields:
    surfH:              a handle row [1, 1/3]. The surface handle [XY, YZ, XZ].
    scatH:              a handle row [1, 2/3]. The scatter handle for dotted lines.
    cbarH:              a colorbar handle. The colorbar for probability.
probInfo:           a struct. It collects all probability-related values in following fields:
    meshCell:           a cell array {1/3, 3}. Each {j,:} includes mesh grids for XY, YZ, or XZ planes.
    dotCell:            a cell row {1, 2/3}. The dotted lines along each axis.
    probY:              a cell row {1, 1/3}. The grid values for {XY, YZ, XZ} planes.
    funProbY:           a func. handle. The p(y) either given or created by this function.
    cbLimit:            a row [1, 2]. The current colorbar limit.
%}

%% Constant.
% Figure dimension (2D/3D).
dimAxis = length(gridCell);

% The limit [min, max] of each axis.
axisLimit = zeros(dimAxis, 2);
for(j = 1:1:dimAxis)
    axisLimit(j,:) = gridCell{j}([1,end]);
end
    
%% Default - fundamental.
% Figure axis.
axisH      = [];
tickGap    = [];

% Scatter.
scatNumDot = 25;
scatSet    = {'SizeData', 10, 'MarkerFaceColor', rgb('Darkgray'), 'MarkerEdgeColor', 'none'};

% Colorbar.
cbMap      = "turbo";
cbSet      = {};
cbLimit    = [];
cbLabel    = "Probability";

%% Default - axis dimension dependent.
switch dimAxis
    case 2  
        % CASE: the embedding space = R^2.
        axisSet     = {'xlim', axisLimit(1,:), ...
                       'ylim', axisLimit(2,:)};
        gridSW      = false;
        viewAngle   = [0, 90];
        scatSW      = false;
        
    case 3
        % CASE: the embedding space = R^3.
        axisSet     = {'xlim', axisLimit(1,:), ...
                       'ylim', axisLimit(2,:), ...
                       'zlim', axisLimit(3,:)};
        gridSW      = true;
        viewAngle   = [-33, 36];
        scatSW      = true;
        
    otherwise
        error('The axis dimension (%d) must be 2 or 3.', dimAxis);
end    
    
%% (if exist) optional inputs.
numVarargin = length(varargin);
if(numVarargin > 0)
    tagInd = 1;
    while(tagInd <= numVarargin)
        % Select the option.
        tagName = varargin{tagInd};
        switch tagName
            case "axisH";       tagShift = 1;       axisH      = varargin{tagInd + 1};
            case "axisSet";     tagShift = 1;       axisSet    = [axisSet, varargin{tagInd + 1}];
            case "tickGap";     tagShift = 1;       tickGap    = varargin{tagInd + 1};
            case "grid";        tagShift = 1;       gridSW     = varargin{tagInd + 1};    
            case "view";        tagShift = 1;       viewAngle  = varargin{tagInd + 1};
            case "scatSW";      tagShift = 1;       scatSW     = varargin{tagInd + 1};
            case "scatNumDot";  tagShift = 1;       scatNumDot = varargin{tagInd + 1};    
            case "scatSet";     tagShift = 1;       scatSet    = [scatSet, varargin{tagInd + 1}];
            case "cbMap";       tagShift = 1;       cbMap      = varargin{tagInd + 1};
            case "cbSet";       tagShift = 1;       cbSet      = [cbSet, varargin{tagInd + 1}];
            case "cbLimit";     tagShift = 1;       cbLimit    = varargin{tagInd + 1};
            case "cbLabel";     tagShift = 1;       cbLabel    = varargin{tagInd + 1};    
            otherwise
                error('The given option "%s" is invalid.', tagName);
        end
        % Shift the index.
        tagInd = tagInd + tagShift + 1;
    end
end

%% Initialize variables for the outputs.
objHColl = struct('surfH', [], 'scatH', [], 'cbarH', []);
probInfo = struct('meshCell', [], 'dotCell', [], 'probY', [], 'funProbY', [], 'cbLimit', []);

%% Compute: "funProbY".
% Select the type of "funProbY".
switch probYType
    case "given"
        % CASE: "funProbY" is given already. Use it directly.
        funProbY = funStru.funProbY;
        
    case "parameter"
        % CASE: Build "funProbY" by "funProbYlZ & funProbZ".
        funProbZ   = funStru.funProbZ;
        funProbYlZ = funStru.funProbYlZ;
        
        % Decide the form.
        probZ = funProbZ(sampZ);
        numInputYlZ = nargin(funProbYlZ);
        switch numInputYlZ
            case 2
                funProbY = @(y) funProbYlZ(y, sampZ) * probZ;
            case 3
                funProbY = @(y) funProbYlZ(y, sampZ, paraStru) * probZ;
            otherwise
                error('The #input of "funProbYlZ" (%d) must be 2 or 3.', numInputYlZ);
        end
        
    otherwise
        error('The given "probYType" (%s) is invalid.', probYType);
end

%% Compute "meshCell" for MATLAB "surf".
% Set variables depending on axis dimension.
switch dimAxis
    case 2
        % CASE: XY mesh.
        numSurf = 1;
        meshCell = {gridCell{1}, gridCell{2}, 0};
        
    case 3
        % CASE: XY, YZ, & XZ mesh.
        numSurf = 3;
        meshCell = {gridCell{1}, gridCell{2}, 0          ; ...
                    0          , gridCell{2}, gridCell{3}; ...
                    gridCell{1}, 0          , gridCell{3}};
        
    otherwise
        error('The axis dimension (%d) must be 2 or 3.', dimAxis);
end

%% Compute "probY" and update "meshCell".
probY    = cell(1, numSurf);

% Compute "probY" for each surface.
for(j = 1:1:numSurf)
    % Concatenate all points of a surface.
    [XX, YY, ZZ] = meshgrid(meshCell{j,:});
    meshPairAll = [XX(:), YY(:), ZZ(:)];

    % Compute the probability and transform it back.
    probY{j} = funProbY( meshPairAll(:, 1:dimAxis) );
    meshCell(j,:) = {squeeze(XX), squeeze(YY), squeeze(ZZ)};
    meshDim = size(meshCell{j,1});
    probY{j} = reshape( probY{j}, meshDim(1), meshDim(2) );
end

%% Build the axis object.
if(isempty(axisH))
    figH = figure('position', [100,50,900,900], 'Color', 'white');
    axisH = subplot(1,1,1,'parent',figH);
end

% Set up the environment.
hold(axisH, 'on');
axis(axisH, 'equal');
box(axisH, 'off');

%% Plot "probY" &  dotted lines.
surfH = gobjects(1, numSurf);
for(j = 1:1:numSurf)
    surfH(j) = surf(axisH, meshCell{j,1}, meshCell{j,2}, meshCell{j,3}, probY{j});
    set(surfH(j), 'EdgeColor', 'none');
end

% Collect the handle.
objHColl.surfH = surfH;

%% Compute the points of dotted lines.
% Decide the #dot in axes.
if( length(scatNumDot) == 1 )
    scatNumDot = ones(1, dimAxis) * scatNumDot;
end

% Compute "dotCell"
dotCell = cell(1, dimAxis);
for(j = 1:1:dimAxis)
    dotCell{j} = linspace(axisLimit(j,1), axisLimit(j,2), scatNumDot(j));
end

%% (If asked) plot the dotted line.
if( scatSW )
    scatH = gobjects(1, dimAxis);
    for(j = 1:1:dimAxis)
        % Prepare zeros for idle axes, and then plot dotted lines one by one.
        dotPoint = zeros(scatNumDot(j), dimAxis);
        dotPoint(:,j) = dotCell{j};
        
        % Select "scatter or scatter3".
        if( dimAxis == 2 )
            scatH(j) = scatter(axisH, dotPoint(:,1), dotPoint(:,2), scatSet{:});
        else
            scatH(j) = scatter3(axisH, dotPoint(:,1), dotPoint(:,2), dotPoint(:,3), scatSet{:});
        end
    end

    % Collect the handle.
    objHColl.scatH = scatH;
end

%% Colorbar.
colormap(axisH, cbMap);
cbarH = colorbar(axisH);

% (If asked) set the colorbar.
if( ~isempty(cbSet) )
    set(cbarH, cbSet{:});
end

% (If asked) set the limit.
if( ~isempty(cbLimit) )
    % CASE: set the colorbar limit.
    caxis(axisH, cbLimit);
else
    % CASE: record the current colorbar limit.
    cbLimit = get(cbarH, 'limit');
end

% Set the colorbar.
cbarH.Label.String = cbLabel;

% Collect the handle.
objHColl.cbarH   = cbarH;

%% Finish the figure.
% The grid.
grid(axisH, gridSW);

% General axis properties.
set(axisH, axisSet{:});

% (If asked) Set the view.
if( ~isempty(viewAngle) )
    set(axisH, 'view', viewAngle);
end

% (If asked) set axis tick based on "tickGap".
if( ~isempty(tickGap) )
    % Prepare ticks for each axis.
    tickRange = zeros(dimAxis, 2);
    tickRange(:,1) = ApproxWithUnit(axisLimit(:,1), "method", "floor", "unit", tickGap);
    tickRange(:,2) = ApproxWithUnit(axisLimit(:,2), "method", "ceil",  "unit", tickGap);
    
    % Set X & Y axis ticks.
    set(axisH, 'xtick', [tickRange(1,1) : tickGap : tickRange(1,2)], ...
               'ytick', [tickRange(2,1) : tickGap : tickRange(2,2)]);
    % If it's 3D, set Z axis, too.
    if( dimAxis == 3 )
        set(axisH, 'ztick', [tickRange(3,1) : tickGap : tickRange(3,2)]);
    end
end

%% Prepare outputs.
probInfo.meshCell = meshCell;
probInfo.dotCell  = dotCell;
probInfo.probY    = probY;
probInfo.funProbY = funProbY;
probInfo.cbLimit  = cbLimit;

end

