cs = crystalSymmetry('-43m'); 
ss = specimenSymmetry('mmm');
o = orientation(new_gmatrix,cs,ss);


%%
n = length(o);
                                        % n =20;
figure('position',[50 50 1200 500])
axesPos = subplot(4,4,1);
plotPDF(o(:,1),Miller(1,1,1,cs),'all',axesPos,'contourf','complete');
% plotPDF(o(:,1),Miller(1,1,1,cs),'all',axesPos,'MarkerSize',3);

jj=1;
hold on
for i=1:1:11
axesPos = subplot(4,4,jj);
plotPDF(o(:,i),Miller(1,1,1,cs),'all','parent',axesPos,'contourf','complete')
% plotPDF(o(:,i),Miller(1,1,1,cs),'all','parent',axesPos,'MarkerSize',3)
text(0,0,num2str(i-1));
jj=jj+1;
end
hold off

