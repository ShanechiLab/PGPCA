function [XX, YY, ZZ, resColl] = GeneMeshgrid(obj, gridType, gridPara, varargin)
%{
This function generates three arrays "XX, YY, & ZZ" for plot the torus by "surf". "gridType" can be:

1.  "numGrid":      gridPara = a row [1, 2] = [numGridOut, numGridIn]. The #lattice along outR/inR.
2.  "given":        gridPara = a cell row {1, 2} = {gridOut, gridIn}. The angular samples along outR/inR.

NOTE:
1.  The type "given" includes the case that users want to sample more densely at some area.
2.  "reInitSW" is for convenience when the "given" grids aren't rounded already.
3.  Output "resColl.gridOut/gridIn" includes the repeated initial point for consistency.
-----------------------------------------------------------------------------------------------------------

Inputs:
--------------- Essential
gridType:       a string "numGrid / given". The mesh grid types. Details are above.
gridPara:       a parameter. It depends on "gridType". Details are above.
--------------- Optional ('tag' + 'value' with (number))
"reInitSW":     (1) (reInitSW) a boolean. (def: true if "gridType = numGrid". Otherwise false.)
                    "true" means the "given" gridOut/gridIn's first point need to be repeated at the end
                      so the torus surface is closed. (e.g., [0, pi/4, pi/2, 3*pi/4] -> add 0 at the end.)

Outputs:
--------------- Essential
XX:             a 2D array [numGridIn, numGridOut]. The X coordinate for 3D torus surface.
YY:             a 2D array [numGridIn, numGridOut]. The Y coordinate for 3D torus surface.
ZZ:             a 2D array [numGridIn, numGridOut]. The Z coordinate for 3D torus surface.
resColl:        a struct. It includes all other information in following fields:
    reInitSW:       a boolean. It's the same as the optional input "reInitSW".
    gridOut:        a column [numGridOut, 1]. The grid points along the outR.
    gridIn:         a column [numGridIn, 1]. The grid points along the inR.
    gridZ:          a 2D array [numGridAll, 2]. All grid points in 2D angular coordinate (Z).
    gridY:          a 2D array [numGridAll, obj.dim]. All grid points in the embedded space (Y).
%}

%% Assign inputs.
% Parser for optional inputs.
optionP = inputParser;
optionP.addParameter('reInitSW',  [],  @islogical);

% Parse the inputs.
optionP.parse(varargin{:});
reInitSW = optionP.Results.reInitSW;

%% Sanity check. 
% (If not given) Initialize "reInitSW" 
if( isempty(reInitSW) )
    if( gridType == "numGrid" )
        reInitSW = true;
    else
        reInitSW = false;
    end
end

% Check: if gridType == "numGrid", reInitSW must be true!
if( gridType == "numGrid" && reInitSW == false )
    warning('Since "gridType == %s", reInitSW is changed to "true".', gridType);
    reInitSW = true;
end

%% Compute: "gridOut & gridIn".
% Select the method.
switch gridType
    case "numGrid"
        % CASE: grid the outR/inR uniformly.
        [~, ~, resColl] = obj.GeneSample("lattice", gridPara);
        gridOut = resColl.lattZOut;
        gridIn  = resColl.lattZIn;
    
    case "given"
        % CASE: the gridOut/gridIn are given. Just collect them.
        gridOut = gridPara{1};
        gridIn  = gridPara{2};
        
    otherwise
        error('The given "gridType" (%s) is invalid.', gridType);
end
    
% Both grid vectors must be columns.
assert( iscolumn(gridOut), '"gridOut" must be a column.' );
assert( iscolumn(gridIn),  '"gridIn" must be a column.' );
    
% (if asked) repeat the initial value at the end.
if( reInitSW )
    gridOut(end+1) = gridOut(1);
    gridIn(end+1)  = gridIn(1);
end

%% Compute: XX, YY, & ZZ.
% Mesh the 2D angular coordinate.
[meshOut, meshIn] = meshgrid(gridOut, gridIn);
dimMesh = size(meshOut);

% Compute the 3D embedded position.
gridZ = [meshOut(:), meshIn(:)];
gridY = obj.Position(gridZ);

% Reshape the matrices.
XX = reshape(gridY(:,1), dimMesh(1), dimMesh(2));
YY = reshape(gridY(:,2), dimMesh(1), dimMesh(2));
ZZ = reshape(gridY(:,3), dimMesh(1), dimMesh(2));

%% Prepare outputs.
resColl = struct('reInitSW',    reInitSW, ...
                 'gridOut',     gridOut, ...
                 'gridIn',      gridIn, ...
                 'gridZ',       gridZ, ...
                 'gridY',       gridY);

end

