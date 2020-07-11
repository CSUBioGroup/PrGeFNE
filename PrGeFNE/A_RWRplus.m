function  [Pt, WAdj ]= A_RWRplus(Adj, r_restart, P0, N_max_iter, Eps_min_change, IsNormalized,  NormalizationType)
% A_RWRplus is a generalization of RWR algorithm.
% Including various propagation algorihtms with initial regularization: classical RWR, Label propagation,and so on
% Including a Solver_IterationPropagation, which can be used directly when IsNormalized is TRUE.
% % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % 
% need getNormalizedMatrix  
% Input 
% Adj
% r_restart
% P0   It should be normalized, though it is neccesary. 
% N_max_iter
% Eps_min_change
% IsNormalized   Whether Adj has been normalized: True or False  
% NormalizationType: including two types of methods by different Normalization Types
% (1) random walk with restart {'ProbabilityNormalizationColumn','ProbabilityNormalizationCol','col','column'}
% (2) 'LaplacianNormalization'   similar to PRINCE(PLOS Computational Biology, 2010, 6: e1000641.)               
% it is equivalent to PRINCE if assigned extended P0 with % disease simialrity and logistic function
% (3)label propagation with memory restart { 'row' ,'ProbabilityNormalizationRow'} 
% % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % 
% Ouput
% Pt  
% WAdj   normalized Adj 
% % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % 
% % % By Ju Xiang, 
% % % Email: xiang.ju@foxmail.com, xiangju@csu.edu.cn  
% % % 2019-3-11   
% % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % 
% % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % 
% % Adj = rand(5); Adj([1 2 ],:)=0.1;   Adj(:, [1 2 ])=0.1;  P0 = [0 0 1 0 0 ]';r_restart=0.7; NormalizationType ='col';
% % Adj = rand(5);   P0 = [0 0 1 0 0 ]';r_restart=0.7; NormalizationType ='col';
% % % Adj = rand(5);   P0 = [0 0 1 0 0; 0 0 0.5 0.5 0 ]';r_restart=0.7; NormalizationType ='col';
    if ~exist('Adj','var')|| isempty (Adj) 
        warning('TestTestTestTestTestTestTestTestTestTestTestTestTestTestTest'); 
        path_tt = 'E:\__DataSet\Dataset_XiangJu\PPI_disease_gene_simulated_data';
        load([path_tt,filesep,'Simulation1-3PPI-Disease-Gene-dataset.mat']); 
        Adj = DataGeneNet_1Net.matrix; 
        AdjDfD = DataDiseaseNet.matrixDisSim ;
        AdjGfD = (DataDiseaseNet.matrixDiseaseGene )' ;   
        %
        idx = find( AdjGfD ); n_pos = length( idx); 
        ind_fold = crossvalind('Kfold', n_pos, 2) ; 
        AdjGfD1 = AdjGfD; AdjGfD1(idx(ind_fold~=1) ) =0; 
        AdjGfD2 = AdjGfD; AdjGfD2(idx(ind_fold==1) ) =0; 
        AdjGfD  = AdjGfD2; % used to train the model.
        DisIDset = 1: size( AdjDfD,1 ); 
        P0 = AdjGfD(:,  DisIDset   );
        NormalizationType = 'col';
    %     NormalizationType = 'row';
    %     NormalizationType = 'LaplacianNormalization';
        r_restart = 0.7;
        istest = true; 
    end
%% % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % 
    if ~exist('N_max_iter','var') || isempty(N_max_iter) || (isnumeric( N_max_iter) && N_max_iter<=1 ) 
        N_max_iter =100; 
    elseif ~isnumeric( N_max_iter)  
        error('N_max_iter should be isnumeric!!!!!') ;
    end
    %
    if ~exist('Eps_min_change','var') || isempty(Eps_min_change) 
        Eps_min_change =10^-6; 
    elseif isnumeric( Eps_min_change) && Eps_min_change>=1 
        warning('The Eps_min_change is nomenaning. Reset Eps_min_change to be 10^-6.'); 
        Eps_min_change =10^-6;  
    elseif ~isnumeric( Eps_min_change)  
        error('Eps_min_change should be isnumeric!!!!!') ;
    end
    
    if ~exist('IsNormalized','var') || isempty(IsNormalized) 
        IsNormalized = false;  % Adj has been normalized for fast run.   
    end
    
    if ~exist('NormalizationType','var') || isempty(NormalizationType) 
        NormalizationType = 'ProbabilityNormalizationColumn'; %%for 'Random Walk' RWR, RWRH  RWRM  RWRMH and more   
    end        
    % % % 取消稀疏矩阵的转换，在调用的外部根据情况指定矩阵形式,除非极端情况 % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % %    
    P0  = full(double(P0)); 
    % 
    % AdjIsSparse = isparse(Adj); 
    AdjDensity =nnz(Adj)/numel(Adj); 
    if  (size(P0,2)==1 && AdjDensity<0.3 ) || (size(P0,2)>1 && AdjDensity<0.05 )
        Adj = sparse(Adj);          
    elseif (size(P0,2)==1 && AdjDensity>0.3 ) || (size(P0,2)>1 && AdjDensity>0.05 )
        Adj = full(Adj); 
    else 
        % no operation 
    end
    %
    if IsNormalized 
        WAdj = Adj; 
    else
        % WAdj = getNormalizedMatrix(Adj, 'col', true );
        switch NormalizationType
            case {'ProbabilityNormalizationColumn','ProbabilityNormalizationCol','col','column'}
                % random walk with restart
                WAdj = getNormalizedMatrix(Adj, 'col', true );
                %%%P0 = P0./(sum(P0,1)+eps);    % total probability is 1. 
                
            case {'LaplacianNormalization'  ,'lap'}
                % propagation similar to PRINCE(PLOS Computational Biology, 2010, 6: e1000641.)
                % it is equivalent to PRINCE if assigned extended P0 with
                % disease simialrity and logistic function
                % A_PRINCEplus is better. 
                WAdj = getNormalizedMatrix(Adj, 'LaplacianNormalization', true );   
                
            case { 'row' ,'ProbabilityNormalizationRow'} 
                % label propagation with memory restart  
                WAdj = getNormalizedMatrix(Adj, 'row', true );  
                
            otherwise
                error(['NormalizationType is wrong: ',char( string(NormalizationType) )]); 
        end        
    end
  
    % sum(WAdj,1)
% %     sss=1111111111
% %     size(WAdj)
% %     size(P0)
    % P0 = reshape( P0,[], 1); 
% % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % %     
% % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % %     
% % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % %     
% % Solver_IterationPropagation
% % It can be used directly when IsNormalized is TRUE.  
    Pt = P0;
    for T = 1: N_max_iter
        Pt1 = (1-r_restart)*WAdj*Pt + r_restart*P0;
        if all( sum( abs( Pt1-Pt )) < Eps_min_change )
            break;
        end
        Pt = Pt1;
    end
    Pt = full(Pt);
% % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % %     
% % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % %     
% % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % %     
% % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % %     
%% % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % 
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
                [X,Y,T,AUC] = perfcurve(label,Pt(idx_test,ii),1);
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
        Pt =[];
    end
    
    
    
end
