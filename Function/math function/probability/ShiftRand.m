function varargout = ShiftRand(varargin)
% 
% This function generates uniform random variables in the same way as the MATLAB "rand", but with shifting.
% 
% Assiging formulation:
% varargout = ShiftRand( [lowBound, upBound], {rest parameters} ): 
%       The first variable (a [1,2] array) determines the lower and upper bounds of the uniform random variable. 
%       The rest parameters are used in MATLAB "rand" directly. The output follows rand's output. 

% set parameters.
lowBound = varargin{1}(1);
upBound = varargin{1}(2);

% generate random variables.
[varargout{1:nargout}] = (upBound-lowBound)*rand(varargin{2:end})+lowBound;

end