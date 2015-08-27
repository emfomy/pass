%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Particle Swarm Stepwise (PaSS) Algorithm                                     %
%                                                                              %
% genlin_power.m                                                               %
% Transform power data                                                         %
%                                                                              %
% Author: emfo<emfomy@gmail.com>                                               %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%% Main function
function genlin_power( srcroot, dstroot )
  if nargin < 1
    srcroot = 'genlin_power.mat';
  end
  if nargin < 2
    dstroot = 'genlin.dat';
  end
  srcname = 'p7';
  dstname = ['General_Linear_Power'];

  % Load data
  data = getfield(load(srcroot, srcname), srcname);
  X = data.x;
  Y = data.y;
  
  % Normalize data
  S = X(:, 1);
  X = normr(X);
  Y = Y .* X(:, 1) ./ S;

  % Get size
  [n, p] = size(X);

  % Generate J
  J = false(1, p);

  % Save data
  file = fopen(dstroot, 'w');

  fprintf(file, '# 1st  line:  data name\n');
  fprintf(file, '# 2st  line:  n p\n');
  fprintf(file, '# 3rd  line:  * J\n');
  fprintf(file, '# rest lines: Y X\n');
  fprintf(file, '# \n');
  fprintf(file, '# X: matrix, n by p, the regressors\n');
  fprintf(file, '# Y: vector, n by 1, the regressand\n');
  fprintf(file, '# J: vector, 1 by p, the chosen indices\n');
  fprintf(file, '# \n');

  fprintf(file, '%s\n', dstname);
  fprintf(file, '%d %d\n', n, p);

  fprintf(file, '%-16c', '*');
  fprintf(file, '%-16d', J);
  fprintf(file, '\n');
  for i = 1:n
    fprintf(file, '%-+16.6e', [Y(i), X(i, :)]);
    fprintf(file, '\n');
  end

  fclose(file);
end
