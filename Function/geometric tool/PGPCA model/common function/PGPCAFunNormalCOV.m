function obseW = PGPCAFunNormalCOV( matR, sampZ, probYlZType, funMatG )
%{
This function computes the Gaussian noise COV in the PGPCA model with any manifold. 

Input "probYlZType" can be:

(1) "EuCOV":    p(y|z) COV follows the extrinsic coordinate, the Euclidean one R^n.
(2) "GeCOV":    p(y|z) COV follows the instinsic coordinate, the geometric one.

Function "funMatG" I/O: 

[dimY, dimY, numZ] = funMatG(sampZ); Each [i,:,j] is the i-th coordinate at "z(j,:)".

NOTE:
1.  This function should keep simple to make it flexible. It's just a basic/auxiliary function.
-----------------------------------------------------------------------------------------------------------

Inputs:
------------------- Essential
matR:               a 2D square [dimY]. The basic COV for Gaussian noise.
sampZ:              a 2D array [numZ, dimZ]. Each row is a sample in the manifold intrinsic coordinate.
probYlZType:        a string "EuCOV / GeCOV". It decides the COV coordinate type. Options are above.
funMatG:            a func. handle. It computes the local coordinate at each z. Details are above.

Outputs:
------------------- Essential
obseW:              a 3D array [dimY, dimY, numZ]. Each [:,:,j] is the COV at "z(j,:)".
%}

%% Sanity check.
% Constant.
numZ = size(sampZ, 1);

% Assert: matR is symmetric.
assert( issymmetric(matR), '"matR" must be symmetric.' );

%% Compute COV.
switch probYlZType
    case "EuCOV"
        % CASE: COV follows the Euclidean coordinate R^n.
        obseW = repmat(matR, 1, 1, numZ);
        
    case "GeCOV"
        % CASE: COV follows the geometric coordinate on the manifold.
        matG = funMatG(sampZ);
        obseW = pagemtimes( pagemtimes(matG, 'transpose', matR, 'none'), matG );
        
    otherwise
        error('The given "probYlZType" (%s) is invalid.', probYlZType);
end

% Make sure that "obseW" is symmetric.
obseW = (obseW + permute(obseW, [2,1,3]))/2;

end

