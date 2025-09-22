function matK = PGPCATorusFunMatK( torusObj, sampZ )
%{
This is an auxiliary function for "PGPCA + Torus" model. It computes K(z) given z.

PGPCA model:  y = g(z) + K(z)*C*x + r (COV: R = (sigma^2)*I_n)

NOTE:
1.  This function should keep simple to make it flexible. It's just a basic/auxiliary function.
2.  We assume the local COV space is full-rank, so dimL = dimY.
-----------------------------------------------------------------------------------------------------------

Inputs:
------------------- Essential
torusObj:           a "TorusFitting" object. It's a user-defined torus class.
sampZ:              a 2D array [numZ, 2]. Each [j,:] = [outAng, inAng] is a 2D torus coordinate.

Outputs:
------------------- Essential
matK:               a 3D array [dimY, dimL, numZ]. Each [:,i,j] is a i-th coordinate vector at sampZ(j,:).
%}

%% Progress.
% Get the local coordinate axes.
matK = torusObj.LocalSpace("complete", sampZ, 'normSW', true);

% Change the index order to fit the output format.
matK = permute(matK, [2,1,3]);

end

