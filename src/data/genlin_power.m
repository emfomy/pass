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
  dstname = ['GenLin_Power', 0];
  dstlen = length(dstname);

  % Load data
  data = getfield(load(srcroot, srcname), srcname);
  X = data.x;
  Y = data.y;
  
  % Normalize data
  S = sqrt(sum(X.^2, 2));
  X = normr(X);
  Y = Y ./ S;

  % Get size
  [n, p] = size(X);

  % Generate J
  J = false(1, p);

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
