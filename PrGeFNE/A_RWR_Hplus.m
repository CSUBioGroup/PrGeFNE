function [P_G, P_D, Pt ] = A_RWR_Hplus(AdjGfG,AdjGfD,AdjDfD, P0_G,P0_D, restart, pro_jump, eta, NormalizationType, isdebug) 
% % % % % % % % % % % % % % % % % % % % 
if nargin<1 || isempty (AdjGfG)
    N_gene=300;
    N_disease = 10; 
    AdjGfG = rand(N_gene,N_gene)>0.2 ; 
    AdjGfD = rand(N_gene,N_disease)>0.2 ; 
    AdjDfD = rand(N_disease,N_disease)>0.2 ; 
    P0_G = zeros(N_gene,1); P0_G(1:3)=1; P0_G = P0_G./sum(P0_G);
    P0_D = zeros(N_disease,1); P0_D(1:2)=1; P0_D = P0_D./sum(P0_D);
    %
    path_tt = 'E:\__DataSet\Dataset_XiangJu\PPI_disease_gene_simulated_data';
    load([path_tt,filesep,'Simulation1-3PPI-Disease-Gene-dataset.mat']); 
    AdjGfG = DataGeneNet_1Net.matrix; 
    AdjDfD = DataDiseaseNet.matrixDisSim ;
    AdjGfD = (DataDiseaseNet.matrixDiseaseGene )' ;   
    %
    idx = find( AdjGfD ); n_pos = length( idx); 
    ind_fold = crossvalind('Kfold', n_pos, 2) ; 
    AdjGfD1 = AdjGfD; AdjGfD1(idx(ind_fold~=1) ) =0; 
    AdjGfD2 = AdjGfD; AdjGfD2(idx(ind_fold==1) ) =0; 
    AdjGfD  = AdjGfD2; % used to train the model.  
    DisIDset = 1: size( AdjDfD,1 ); 
    DisIDset = [1:10]; 
    P0_G = AdjGfD(:,  DisIDset   );   
    P0_D = speye( size( AdjDfD) );  P0_D=P0_D(:,DisIDset); 
   
    restart = 0.7;   pro_jump = 0.5;  eta =0.5;   
    warning('TestTestTestTestTestTestTestTestTestTestTestTestTestTestTest'); 
    isdebug =true; 
    istest =true; 
     NormalizationType = 'ProbabilityNormalizationColumn'; %%for 'Random Walk' RWR, RWRH  RWRM  RWRMH and more   
%     NormalizationType = 'LaplacianNormalization'; %%  for prince and more....    
end
% % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % 
% Input % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % %
% AdjGfG: associatins from (f) genes (G) to Genes (G)  
% AdjGfD: associatins from Diseases (D) to Genes (G) GfD
% AdjDfD  associatins from Diseases (D) to Disease (G) 
% P0_G: column vector (set) initial probabilities in Gene network
% P0_D: column vector (set) initial probabilities in Disease network
% P0_G and P0_D must have the same # of columns. 
% gamma/restart: restarting Probability  
% pro_jump: jumping Probability between different networks
% eta: ratio of Probability in second network
% NormalizationType = 'ProbabilityNormalizationColumn'; %%for 'Random Walk' RWR, RWRH  RWRM  RWRMH and more   
% NormalizationType = 'LaplacianNormalization'; %%  for label propagation, prince and more....    
% Ouput % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % %
% P_G: stable probabilities in Gene network 
% P_D: stable probabilities in Disease network  
% Pt: stable probabilities in Gene+disease network 
% % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % %
% A_RWR_Hplus: include dynamical processes: RWR, label propagation
% RWRH random walk with restart algorithm in heterogeneous network. 
% According to reference: Yongjin Li, Jagdish C. Patra, Genome-wide inferring gene¨Cphenotype relationship 
% by walking on the heterogeneous network, Bioinformatics, 2010, 26: 1219-1224.
% Different from lapRWRH used 'LaplacianNormalization'; for label propagation
% Zhi-Qin Zhao, et al, Laplacian normalization and random walk on heterogeneous networks for disease-gene prioritization, 
% Computational Biology and Chemistry, 2015, 57: 21-28.
% By Ju Xiang, 
% Email: xiang.ju@foxmail.com, xiangju@csu.edu.cn  
% 2019-3-11   
    % global   Global_Var_RWRH 
    if ~exist('pro_jump','var') || isempty (pro_jump)
        pro_jump = 0.5; 
    end   
    if ~exist('eta','var') || isempty (eta)
        eta = 0.5; 
    end   
    if ~exist('NormalizationType','var') || isempty (NormalizationType)
        NormalizationType = 'ProbabilityNormalizationColumn' ; 
    elseif ~ismember(NormalizationType,{ 'LaplacianNormalization', 'column','col',  'ProbabilityNormalizationColumn','ProbabilityNormalizationCol' , 'row','ProbabilityNormalizationRow','NormalizationRow'   } )
        error(['NormalizationType is wrong: ',char(string(NormalizationType)) ]);
    end   
   if ~exist('isdebug','var') || isempty (isdebug)
        isdebug = false;  
    end   
    %  
    if isempty( AdjDfD )
       AdjDfD = speye(N_disease);  
       warning('AdjDfD is empty.');
    end
    % 
    [N_gene, N_disease] = size( AdjGfD );
    if size(P0_G,1)~=N_gene; error( 'P0_G must be column vector(set), with length of the number of genes.'  );end
    if size(P0_D,1)~=N_disease; error( 'P0_D must be column vector(set), with length of the number of diseases.'  );end
    %
    [ M , IsNormalized ] = getNormalizedMatrix_Heter(AdjGfG,AdjGfD,AdjDfD, pro_jump,  NormalizationType, []) ; 
% %     %
% %     idxDis_WithDiseaseGene =  sum( AdjGfD, 1)~=0;   % mark diseases with disease-genes
% %     idxGene_WithDisease    = (sum( AdjGfD, 2)~=0)';   % mark genes that are associated with diseases
% %     %
% %     % WAdj = getNormalizedMatrix(Adj, NormalizationType, SetIsolatedNodeSelfLoop )
% %     M_GfG = getNormalizedMatrix(AdjGfG   , NormalizationType, true ); 
% %     M_DfD = getNormalizedMatrix(AdjDfD   , NormalizationType, true ); 
% %     M_GfD = getNormalizedMatrix(AdjGfD   , NormalizationType, false );  % probabilities from disease space to gene space 
% %     M_DfG = getNormalizedMatrix(AdjGfD'  , NormalizationType, false );  % probabilities from gene space to disease space
% %     %
% %     M_GfG(:,idxGene_WithDisease)       = (1-pro_jump).*M_GfG(:,idxGene_WithDisease); 
% %     M_DfD(:,idxDis_WithDiseaseGene )   = (1-pro_jump).*M_DfD(:,idxDis_WithDiseaseGene ) ; 
% %     M_GfD                           = pro_jump.*M_GfD; % Disease-columns without disease-genes is all zeros. So no use idxDis_WithDiseaseGene
% %     M_DfG                           = pro_jump.*M_DfG; % Gene-columns without diseases is all zeros. So no use idxGene_WithDisease
% %     %
% %     M = [ M_GfG, M_GfD; M_DfG, M_DfD    ] ;
% %     IsNormalized = true; 
    % 
    % %     P0 = [ (1-eta)*P0_G; eta*P0_D]; 
% %     n_DisID = length( DisIDset ); 
% %     P0_G = getNormalizedMatrix( reshape(P0_G,[],n_DisID)   , 'col', false  );  
% %     P0_D = getNormalizedMatrix( reshape(P0_D,[],n_DisID)   , 'col', false  );  
%     P0_G = getNormalizedMatrix(   P0_G  , 'col', false  );  
%     P0_D = getNormalizedMatrix(   P0_D  , 'col', false  );  
    if any( strcmpi(NormalizationType,{'col','ProbabilityNormalizationColumn','NormalizationColumn', 'Column'}) )
        P0_G = P0_G./(sum(P0_G,1)+eps);  
        P0_D = P0_D./(sum(P0_D,1)+eps);  
    end
    P0 = [ (1-eta)*P0_G; eta*P0_D];   
    % [Pt, WAdj ]= A_RWRplus(Adj, r_restart, P0, N_max_iter, Eps_min_change, IsNormalized,  NormalizationType)
    Pt = A_RWRplus(M, restart, P0 , [],[], IsNormalized);   
    P_G = Pt(1:N_gene,:);
    P_D = Pt(N_gene+1:end,:);
% %     Pset =reshape(Pt,n_node,L_net)  
% % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % 
% % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % 
% % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % 
% % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % 
% % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % 
    if exist('istest','var') && istest  
        AUCset = [];  
        AdjGfD_test  = AdjGfD1(:, DisIDset  );  
        AdjGfD_train = AdjGfD2(:, DisIDset  );  
        %AdjGfD_initial = AdjGfD2(:, DisIDset  );  
        %AdjGfD_all = AdjGfD(:, DisIDset  );  
        for ii = 1:size( AdjGfD_test, 2)
            idx_test = AdjGfD_train(:,ii)~=1;   
            label = AdjGfD_test(idx_test,ii); 
            if any( label) 
                [X,Y,T,AUC] = perfcurve(label,P_G(idx_test,ii),1);
            else
    %             nnz( AdjGfD_initial(:,ii) )
                AUC = -1;
            end
%             [ nnz( AdjGfD_initial(:,ii) )   nnz( AdjGfD_all(:,ii) ) ] 
            AUCset(ii,1) = AUC  ;
        end
        disp('*******************************************************')
        disp( AUCset'  )  
        AUCsetmean =mean(AUCset(AUCset~=-1))       
    end

    if isdebug
%         sum(M)  
        disp('*******************************************************')
        allP=sum(Pt) 
        ss = sum(M(:)) 
        size(M)  
        P_G=[]; P_D=[];Pt=[] ; 
    end
end
