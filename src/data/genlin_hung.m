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
    dstroot = 'genlin.dat';
  end
  dstname = ['GenLin_Hung_', srcname, 0];
  dstlen = length(dstname);

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
  k = 1;
  Idx = zeros(1, size(True_main, 2)+size(True_int, 2));
  for idx = True_main
    i = find(X_ind == idx(1));
    if ~isempty(i)
      Idx(k) = i;
      k = k+1;
    end
  end
  for idx = True_int
    i = find(X_ind == idx(1));
    j = find(X_ind == idx(2));
    if ~isempty(i) && ~isempty(j)
      Idx(k) = i*q - i*(i+1)/2 + j;
      k = k+1;
    end
  end
  Idx(k:end) = [];

  % Generate J
  J = false(1, p);
  J(Idx) = true;

  % Save data
  file = fopen(dstroot, 'wb');
  fwrite(file, dstlen, 'integer*4');
  fwrite(file, dstname, 'char*1');
  fwrite(file, n, 'integer*4');
  fwrite(file, p, 'integer*4');
  fwrite(file, X, 'real*4');
  fwrite(file, Y, 'real*4');
  fwrite(file, J);
  fclose(file);
end
