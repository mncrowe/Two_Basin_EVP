function [s,x] = Gen_EVP(A,B,n,s0,method)
% Solves the generalised EVP (s*A + B)*x = 0 for eigenvalue s and eigenvector x
%
% Inputs:
% - A,B: square matrices
% - n: number of eignevalues to find, default: 10
% - s0: initial guess for eigenvalue, default: 0
% - method: search method for eigs, default: 'lm'
%
% Note: this script finds the n eigenvalues closest to s0 by writing
% [(s0+s')*A + B]*x = 0  =>  [(1/s')*(s0*A+B) + A]*x = 0,
% and solving for eigenvalue 1/s'.

if nargin < 3; n = 10; end
if nargin < 4; s0 = 0; end
if nargin < 5; method = 'lm'; end

% determine eigenvalues and eigenvectors of modified system
[V,S] = eigs(sparse(A),sparse(-B-s0*A),n,method);

% convert to eigenvalues of original system and rescale eigenvectors
s = s0 + 1./sum(S).';
x = reshape(V,[length(A) n])./max(reshape(V,[length(A) n]));

% remove NaN and inf and sort in descending order by real part
x = x(:,and(~isnan(s),~isinf(s)));
s = s(and(~isnan(s),~isinf(s)));

[~,I] = sort(real(s),'descend');
s = s(I);
x = x(:,I);

end

