function [] = main_lptayl_function(filename,incremental_strain,n_iteration,varargin)

%% 
% Syntax:
% main_lptayl_function('myfile.txt',0.05,20,'relaxConstraint','1ath','slip_type',...
% 'fp')



%%
g_file = [filename, '.txt']; 
% prompt = 'Please enter crystal strucutre f(fcc) and b(bcc)\n';
% cryst_str = input(prompt);
prompt = 'Please select type of slip\n f (full), p (only partial), fp (full+partial), ft(full+twin)\n';
slip_type = 'f';
relaxConstraint = 0;
% prompt = 'Please select constraints: 0 (Full)\t 1 (Lath)\t 2 (PanCake)\n';
% constraints = input(prompt);
% prompt = 'Please enter maximum strain e.g. 2,3 etc \n';
% max_strain = input(prompt);
% prompt = 'Please enter incremental strain e.g. 0.05,0.1 etc \n';
% incr_e = input(prompt);
% prompt = 'Please select criteria to solve taylor ambiguity \n 0 (minVar), 1(maxVar), 2(minPlasticSpin), 3(wintenberger) \n';
criteria = input(prompt);
tic



%% Initializing the cost matrix according to type of slip mode
c = ones(1,48);
% full slip: n = 1-12,25-36 % partial slip: n=13-24,37-48 % twin:n=13-24

switch slip_type
    case {'F', 'f'} % only perfect slip
        n = [13:24,37:48];
        c(n)=1e4; % All partial modes are made costlier        
    case {'P','p'} % only partial slip
        n = [1:12,25:36];
        c(n)=1e4; % All perfect modes are made costlier
    case {'FP','fp','Fp','fP'} %both perfect and partial slips
        
    case {'FT','ft','Ft','fT'}
        n = 37:48;
        c(n)=1e4; % All - partial modes are made costlier        
    otherwise
        warning('Unexpected slip type.\n {Default perfect slip is applied');        
 end
        


%% Initialization

% Criteria will be given as input in later development. But for present
% minvar is taken as the single unique solution selection criteria. Others
% are 'maxvar' and 'minplasticspin'
% criteria = 'minvar';
% criteria = 'maxvar';
% criteria = 'minplasticspin';

SlipSystem = SS_setFCC_fp_function();
A = getA(SlipSystem);
gmatrix0 = txtfile2ori(g_file);
gmatrix = gmatrix0;
lg =  length(gmatrix);
e_ext= [1,0,0;0,0,0;0,0,-1];
new_gmatrix(:,1) = gmatrix;

%% Calculation Block

for j=1:1:n_iteration

for i=1:1:lg 
    DCMtrx = matrix(gmatrix(i));  
    e_grain=DCMtrx'*e_ext*DCMtrx; 
    b = [e_grain(1,1);e_grain(2,2);2*e_grain(2,3);2*e_grain(1,3);2*e_grain(1,2)]; 
    AllSolutions = calcLPSolution(b,c,A); 
    [spin, UniqueSolution] = calculate_spin(AllSolutions,criteria,SlipSystem);
    spin = incremental_strain*spin;
    gmatrix(i) = applySpin(-spin,DCMtrx); 
    new_gmatrix(i,j+1) = gmatrix(i);
    ActiveSS(i,j).ss= UniqueSolution.B;
    ActiveSS(i,j).gamma = UniqueSolution.xb;
end
end

%% Output Block
prompt = 'Do you want to save results; yes or no \n';
reply = input(prompt,'s');
if strcmp(reply,'yes')
    folderpath = resmet_input(name,new_gmatrix);

    ActiveSSData = [folderpath '\' 'ActiveSSData.mat'];
    save(ActiveSSData,'ActiveSS');
    
    %%%Following lines to save o/p euler angles in file
    %     eu_gmat = Euler(new_gmatrix)/degree;
    %     euler_angles(name,new_gmatrix);
else
    return
end
%%

toc

