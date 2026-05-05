function m_vec = mass_vector(N, delta, lambda)
% MASS_VECTOR  Floor mass profile for variable-mass buildings.
%   Computes a floor mass vector using the same parameterized reduction
%   formula as stiffness_vector, applied to mass instead of stiffness.
%   Models buildings where upper floors are lighter than lower floors.
%
%   Formula: m(j) = 1 - (1-delta) * (j/N)^lambda,  for j = 1..N
%
% INPUTS:
%   N : number of stories (scalar integer)
%   delta: mass ratio at the roof relative to the base (scalar, 0 < delta ≤ 1)
%            e.g. delta=0.2 means the roof floor has 20% of the base floor mass
%   lambda: exponent controlling the shape of mass reduction (scalar > 0)
%            lambda=1: linear,  lambda=2: parabolic
%
% OUTPUT:
%   m_vec: (N x 1) floor mass vector; m_vec(j) is the mass at floor j,
%            where j=1 is the ground floor and j=N is the roof
    j_vec = (1:N)';
    m_vec = 1 - (1 - delta) .* (j_vec / N).^lambda;
end