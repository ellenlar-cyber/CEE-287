function k_vec = stiffness_vector(N, delta, lambda)
% STIFFNESS_VECTOR  Story stiffness profile for variable-stiffness buildings.
%   Computes a story stiffness vector using Miranda's parameterized reduction
%   formula, where stiffness decreases from base to roof
%   Formula: k(j) = 1 - (1-delta) * (j/N)^lambda,  for j = 1..N
%
%   Special cases:
%     delta = 1  →  uniform stiffness k(j) = 1 for all j 
%     lambda = 1 →  linear reduction from base to roof
%     lambda = 2 →  parabolic reduction (stiffness stays high in lower
%                   stories and drops sharply near the roof)
%
% INPUTS:
%   N      - number of stories (scalar integer)
%   delta  - stiffness ratio at the roof relative to the base (scalar, 0 < delta ≤ 1)
%            e.g. delta=0.2 means the roof story is 20% as stiff as the base
%   lambda - exponent controlling the shape of stiffness reduction (scalar > 0)
%            lambda=1: linear,  lambda=2: parabolic
%
% OUTPUT:
%   k_vec  - (N x 1) story stiffness vector; k_vec(j) is the stiffness of
%            story j, where j=1 is the ground story and j=N is the top story
    j_vec = (1:N)';
    k_vec = 1 - (1 - delta) .* (j_vec / N).^lambda;
end