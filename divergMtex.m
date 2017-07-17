

clear;

cs = crystalSymmetry('cubic');
ss = specimenSymmetry('mmm');
sS = slipSystem.fcc(cs);
phi1_id = 30;
phi_id = 45;
phi2_id = 30;
range = 10;
spacing = 2.5;
eul_ini = grid2euler(phi1_id, phi_id, phi2_id,range,spacing);


ori = orientation('Euler',eul_ini,cs);%,ss);

% 30 percent strain
q = 0; % Rolling velocity gradient
% epsilon = 0.3 * tensor.diag([1 -q -(1-q)],'name','strain');
epsilon = 0.02 * tensor(diag([1 -q -(1-q)]),'name','strain');
epsilon.CS = cs;

%
numIter = 1;
progress(0,numIter);

for sas=1:numIter

  % compute the Taylor factors and the orientation gradients
  [M,~,mori] = calcTaylor(inv(ori) * epsilon ./ numIter, sS.symmetrise,'silent');% [M,b,mori]

  % rotate the individual orientations
  ori = ori .* mori;
  progress(sas,numIter);
end
progress= 100;%
% plot the resulting pole figures

plotPDF(ori,Miller({1,1,1},{0,0,1},cs),'contourf','complete')
mtexColorbar