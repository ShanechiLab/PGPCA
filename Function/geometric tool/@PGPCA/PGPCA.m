classdef PGPCA < matlab.mixin.Copyable
%{
This class "probability geometric PCA (PGPCA)" computes PPCA-like dimensionality reduction on manifold.

To make it generalizable, we summarize the I/O format in this class, but all functions are PUBLIC. In this 
  way, this class is just a toolbox collecting all functions we need in one place.
-----------------------------------------------------------------------------------------------------------

In detail, this toolbox learns the following model from "data":

Model:  y = g(z) + K(z)*C*x + r (COV: R = (sigma^2)*I_n)

where:  z ~ manifold coordinate with a probability distribution.
        x ~ N(0, I_m), a multivariate normal distribution.

NOTE:   
1.  All y_t, z_t, & x_t are independent of other time index pair (y_i, z_i, & x_i)
-----------------------------------------------------------------------------------------------------------

Since it's data-based, there are multiple steps, expecially the cross-validation (CV)!

1.  E-step:     Find the representing value z_t^* for each y_t by intermediate value theorem (IVT).
2.  M-step:     Solve optimal parameter "C and sigma^2" with various conditions.
    
NOTE:   
1.  Unlike PCA, LPR, and PPCA, space partition is not a bonus property of PGPCA anymore (or at least we 
      cannot prove it!) However, there is no approximation in LL of PGPCA (E-step's IVT is precise), so
      this is a "characteristic" of "LL + this manifold model", not a flaw! Moreover, PGPCA still find 
      the optimal parameters "given the model dimension", so it still accomplishes the original goal. 
-----------------------------------------------------------------------------------------------------------

To compute cost LL in E-step, M-step, and the complete EM iteration most generally, the conditional prob.
  of "y" w.r.t "z" must follow the format below.
    
p(y|z) -> [probVal, resColl] = funProbYlZ( sampY, sampZ, paraStru, funStru );
    
NOTE:
1.  Letter "l" in name "funProbYlZ" is to replace "|" in p(y|z) symbol.   
2.  Inputs "paraStru & funStru" are structs whose fields are "system parameters/variables" & "relavant
      functional handles", respectively. This I/O includes most of probability models.
3.  Outputs "probVal" is a probability-related value/array (like LL). "resColl" collects the rest items. 
-----------------------------------------------------------------------------------------------------------
    
Since this class is just a toolbox, not a real PGPCA "object", there is no hard-coded properties and I/O
  rules. Still, in principle, all functions' "ESSENTIAL" I/O should follow rules below for consistency:

Inputs:
------------------- (RANDOM VARIABLE)
sampY:              a 2D array [numY, dimY]. Every row is an observation y_j.
sampZ:              a 2D array [numZ, dimZ]. Every row is an integral point in computing p(y).
------------------- (PGPCA VARIABLE/SETTING)
matPi:              a symmetric matrix [dimY]. The covariance of Y in "all dimensions".
matGa:              a symmetric matrix [dimL]. The covariance of Y in "K matrix subspace".
dimC:               an integer >= 0. The dimension of "matC", the latent state dimension.
------------------- (SYSTEM PARAMETER)
paraStru:           a struct. It contains all (learned) system parameters below.
    matC:               a 2D array [dimL, dimC] / empty []. The dimensionality reduction matrix.
    mainVar:            a scalar >= 0. The average noise VAR of the "modeled space".
    sideVar:            a scalar >= 0 / empty []. The average noise VAR of the "unmodeled space".
------------------- (SYSTEM FUNCTION)
funStru:            a struct. It collects all function handles in following fields:
    funValYZ:           a func. handle. The conditional value valYZ(y|z) (may not be probability!)
    funProbZ:           a func. handle. The probability p(z).
    funProbYlZ:         a func. handle. The probability p(y|z).
    funMean:            a func. handle. The mean vector "phi(z)" for each "z" sample.
    funMatK:            a func. handle. The "orthonormal" coordinate matrix "K(z)" for each "z" sample.              
-----------------------------------------------------------------------------------------------------------
  
Outputs:
------------------- (BASIC VARIABLE)
*NAME LIKE INPUT*:  It's updated by functions. The class/variable type is basically the same.
probVal:            a 2D array [numY, numZ] / a column [numY, 1]. A probability-relevant array/value.
------------------- (ADDITIONAL INFO.)
resColl:            a struct collecting all additional information in a function. It's mainly the last
                      output and useful in further manipulation/debugging.
res + *NAME*:       a struct similar to "resColl". It mainly collects "resColl" from other functions.
%}
%% Class definition.
    properties (Constant)   
        % The "relative" error ratio. (e.g., abs(a-b) < errRatio*abs(a) -> "a == b" (good enough))
        errRatio = 1e-8;
        % The "practical" probability limit to be considered = 0 due to numerical error.
        zeroProbLimit = 1e-50;
    end
    
    methods (Static)
        %{
        NOTE:
        1.  Most sub-functions are still public (no access protection). This is helpful in quick test. 
              But for the final result, we still suggest using the clean packaged U/I function.
        %}
        
        %% Main U/I functions for users.
        
        % The EM algorithm of PGPCA. It's the main user interface.
        varargout = EMAlgo(varargin);
        
        %% Component functions in PGPGA.
        
        % The E-step of PGPCA. Its I/O can handle general manifold.
        varargout = EStepAll(varargin);
        % The M-step of PGPCA. Its I/O can handle various model conditions.
        varargout = MStepAll(varargin);
        % The preparing function transfers modeling parameters & functions into variables ready-to-use.
        varargout = ParaFun2Model(varargin);
        
        %% Analyzing tools for PGPCA.
        
        % Analyze parameter's convergence.
        varargout = PlotParaConverge(varargin);
        % Plot the probability distribution in 2D/3D for visualization
        varargout = PlotProbDensity(varargin);   
        
        %% Accessory functions in PGPCA.
        
        % Multiple covariance of the PGPCA model.
        varargout = ModelCov(varargin);
        % The general transfer function of "modeling parameters/functions -> probability".
        varargout = ModelPara2Prob(varargin);
        % The paired conditional probability p(sampY|sampZ), where #sampY == #sampZ.
        varargout = PairCondiProb(varargin);
        % The evidence lower bound (ELBO) of log-likelihood under different conditions.
        varargout = ELBOForEM(varargin);
        
        % Compute the "low-D" system (may not be optimal) degenerated from the given "high-D" system.
        varargout = ModelParaHD2LD(varargin);
        
    end
    
    methods (Static, Access = protected)
        %% These are tools for those static & public functions.
        
        
        
    end
end