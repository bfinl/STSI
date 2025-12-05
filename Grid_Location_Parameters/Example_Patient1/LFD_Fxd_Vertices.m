% CURRY 8 M-file, Platform: PCWIN, Created on: 3/2/2018 5:20:36 PM
%
% curryloc contains mesh locations
% currytri contains triangle vertex indices
% currylfd(1:3,:) contains leadfield locations
% currylfd(4:6,:) contains leadfield orientations
% currylfd(7:end,:) contains leadfield
%
% load mat file
load ('LFD_Fxd_Vertices.mat');
% number of locations
nLoc = length ( curryloc );
% number of values and components per value
nTot = length ( currylfd );
nCom = ceil ( nTot / nLoc );
nVal = nTot / nCom;
% channel (>= 7) to be plotted
nRow = 7;
% basis vector (1..nCom) to be plotted
nBas = 1;
% prepare vector of values to be plotted (zero-padded,transposed)
V = zeros ( nLoc, 1 );
V(1:nVal,1) = currylfd(nRow,nBas:nCom:end)';
% plot using patch command
clf ( 'reset' );
nMin = min(curryloc,[],2)-10;
nMax = max(curryloc,[],2)+10;
axis ( [nMin(1),nMax(1),nMin(2),nMax(2),nMin(3),nMax(3)] );
axis equal;
hpatch = patch ( 'vertices',curryloc','faces',currytri','FaceVertexCData',V );
% get Matlab/Octave version
v = ver;
if strcmp(v(1).Name, 'MATLAB')
  set ( hpatch,'EdgeColor','none','FaceColor','interp','FaceLighting','phong','DiffuseStrength',0.8 );
  axis vis3d;
  camlight right;
  lighting phong;
else % Octave
  set ( hpatch,'EdgeColor','none','FaceColor','interp','DiffuseStrength',0.8 );
end
colorbar;
