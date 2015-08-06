%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Particle Swarm Stepwise (PaSS) Algorithm                                     %
%                                                                              %
% genlin_hung.m                                                                %
% Transform Hung et al.'s data                                                 %
%                                                                              %
% Author: emfo<emfomy@gmail.com>                                               %
%                                                                              %
% Reference:                                                                   %
% Hung, H., Chen, P.-W., Wang, C.-C., Huang, S.-Y., & Tzeng, J.-Y. (2013).     %
%   Detection of Gene-Gene Interactions using Multistage Sparse and Low-Rank   %
%   Regression. http://arxiv.org/abs/1304.3769                                 %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%% Main function
function genlin_hung( srcname, srcroot, dstroot )
  if nargin < 1
    srcname = 'data1_1';
  end
  if nargin < 2
    srcroot = 'genlin_hung.mat';
  end
  if nargin < 3
    dstroot = '../../run/genlin.dat';
  end
  dstname = ['GenLin_Hung_', srcname, 0];

  % load data
  data = getfield(load(srcroot, srcname), srcname);
  X = data.x;
  Y = data.y;
  X_ind = data.x_ind;
  True_main = data.true_main;
  True_int = data.true_int;

  % Get size
  [n, q] = size(X);
  p = q*(q+1)/2;

  % Generate tensor effects of X
  k = q;
  X = [X, zeros(n, p-q)];
  for i = 1:q
    for j = i+1:q
      k = k+1;
      X(:, k) = X(:, i) .* X(:, j);
    end
  end

  % Convert indices
  k = 0;
  Idx = zeros(1, size(True_main, 2)+size(True_int, 2));
  for idx = True_main
    k = k+1;
    Idx(k) = find(X_ind == idx);
  end
  for idx = True_int
    k = k+1;
    i = find(X_ind == idx(1));
    j = find(X_ind == idx(2));
    Idx(k) = i*q - i*(i+1)/2 + j;
  end

  % Generate J
  J = false(1, p);
  J(Idx) = true;

  % Save data
  file = fopen(dstroot, 'wb');
  fwrite(file, length(dstname), 'integer*4', 'ieee-be.l64');
  fwrite(file, dstname, 'char*1', 'ieee-be.l64');
  fwrite(file, n, 'integer*4', 'ieee-be.l64');
  fwrite(file, p, 'integer*4', 'ieee-be.l64');
  fwrite(file, X, 'real*4', 'ieee-be.l64');
  fwrite(file, Y, 'real*4', 'ieee-be.l64');
  fwrite(file, J);
end
