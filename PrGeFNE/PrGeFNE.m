% % function RW_RDGN(data_path, cv_subdir , topk,        file_net, file_node ,file_emb   )
% % function RW_RDGN( topk, file_dis_gene_net, file_node ,file_emb, file_pred_results, knn   )
function [ResScalartable, ResCurvetable,ResTopKtable, disease_gene_tableset ] = PrGeFNE( topk_save, fileset , knn, StartID, path_RC , NormalizationType, topk_perf_set , DisNetType, GeneNetType, UseMethodType  )
% % UsePrinceStrategy = 0;   % ÊÇ·ñ Ê¹ÓÃPRINCE ²ßÂÔ ´«²¥ ¼²²¡ÏàËÆÐÔÐÅÏ¢ÓàPPIÍøÂç 
if ~exist('topk_save','var')
    topk_save = 100;  
    % fileset.gene_gene = 'dataset\blab_ppi2016.txt'; 
    fileset.dis_gene = 'data_CV_DisGene_Authors\cv10_of0.txt'; 
    fileset.nodes    = 'data_OrigHerNet\DGG\nodes.txt'; 
    fileset.edges    = 'data_OrigHerNet\DGG\edges.txt'; 
    fileset.embtxt   = 'data_OrigHerNet\DGG\emb.txt';  
    fileset.results  = 'data_CV_Results_EmbOrigNet_author\DGG\cv10_of0\tempRW_prediction_results.txt';  
    tic  
end
if ~exist('StartID','var')
    StartID = 0;  
end
%%
if ~exist('knn','var')|| isempty( knn )
% %     alpha = 20   % ï¿½ï¿½ï¿½ï¿½ï¿?ï¿½ï¿½ï¿½ï¿½ 
% %     beta = 100    
    knn.disease = 20; 
    knn.gene    = 100; 
end
if ~exist('path_RC','var')|| isempty( path_RC )
    path_RC = 'temp_RC';  
    if  ~exist(path_RC,'dir'); mkdir( path_RC ); end
end
if ~exist('topk_save','var')
    topk_save = 100;          % used for saving dis-gene list 
end
% if ~exist('topk_perf_set','var')
%     topk_perf_set = [];        % used for perfcurves
% end
if ~exist('DisNetType','var') || isempty(DisNetType)
    DisNetType = 'None'; 
end
if ~exist('GeneNetType','var') || isempty(GeneNetType)
    GeneNetType = 'None'; 
end

%%
%
% file_gene_gene    = APStr( fileset.gene_gene ); 
file_dis_gene_net = APStr(fileset.dis_gene ); 
file_node         = APStr(fileset.nodes ); 
file_edge         = APStr(fileset.edges ); 
file_emb          = APStr(fileset.embtxt ); 
file_pred_results = APStr(fileset.results );  

% % % % % % % % % % % % % % % % %  
% out_file_dis_dis_net   = [path_RC,filesep,'RC_dis_dis_net.txt'] ;
% out_file_gene_gene_net = [path_RC,filesep,'RC_gene_gene_net.txt'] ;
% % out_file_dis_gene_train= [path_RC,filesep,'RC_dis_gene_train.txt'] ;
% % % % % % % % % % % % % % % % % % % % % % % % % % 
% ï¿½ï¿½È¡ feature vectors ï¿½ï¿½È¡ emb ï¿½Ä¼ï¿½
% % firstline = dlmread(file_emb,' ',[0,0,0,1]) ; % ï¿½ï¿½ï¿½ï¿½ï¿½ï¿½ï¿½ï¿½ï¿½ï¿½ï¿½ï¿½ï¿½ï¿½ï¿½ï¿½ï¿½ï¿½ï¿½ï¿½ 
% % n_rows_emb = firstline(1);
% % n_features = firstline(2);
% node_features = dlmread(file_emb,'',1,0) ;  
% % nodeIDs_emb = strtrim(  cellstr( num2str( node_features(:,1) )  )  ); 
% % nodeIDtable = array2table((1:n_rows_emb)', 'RowNames', nodeIDs_emb );  
% %   
TableNode=readtable( file_node ,'ReadVariableNames',false ,  'Format', '%s%d%s');
n_node   = size(TableNode,1); 
ind_gene = find( strcmp(TableNode{:,3},'gene')            ) ; 
ind_dis  = find( strcmp(TableNode{:,'Var3'},'disease')    ) ;   
% % Ö¸ï¿½ï¿½ï¿½ï¿½ï¿½ï¿½ï¿½Ð±ï¿½ 
% gene_table  = TableNode( ind_gene ,:); n_gene = size( gene_table  , 1 ); 
% % gene_id_set = gene_table{:,2};              %   length(gene_id_set)
% disease_table  = TableNode( ind_dis,:); n_dis = size(disease_table,1 );
% % disease_id_set = disease_table{:,2};        %   length(disease_id_set)
% 
edges = load( file_edge )-StartID+1; 
edges = [edges; edges(:,[2 1 ] ) ] ; 
Mat = ( sparse( edges(:,1), edges(:,2), 1,  n_node, n_node )>0  );   
%
TableNode.Properties.RowNames = TableNode{:,1};          % for extracting numeric ID of disease 
disease_table                 = TableNode(ind_dis , :);  disease_table.IDnew = (1:size(disease_table,1))';  
gene_table                    = TableNode(ind_gene , :); gene_table.IDnew    = (1:size(gene_table,1))';
% 
node_features = []; 
keepweight = true; 
switch lower(DisNetType)
    case lower({'None','NoERC'})
        Mdd = Mat(ind_dis,ind_dis);
        Mdd(logical(speye( size(Mdd) )) )= 1; 
% %         size(Mdd  )
% %         ss  = sum( Mdd(:)  ) 
% %         if sum( Mdd(:)  )==0; Mdd = speye(Mdd); end          
        
    case lower('DisERC')
        if isempty( node_features ) 
            node_features = dlmread(file_emb,'',1,0) ; 
        end
        % node_type1 = 'disease';  %node_type2 = 'disease'; 
        [Mdd ]=getLVRReconstructedNet(node_features,TableNode, disease_table, [], knn.disease, []  , keepweight);         
        
    case lower('DisFN')  
        if isempty( node_features ) 
            node_features = dlmread(file_emb,'',1,0) ; 
        end 
        Mdd = Mat(ind_dis,ind_dis);   % if sum( Mdd(:)  )==0; Mdd = speye(Mdd); end  
        Mdd(logical(speye( size(Mdd) )) )= 1; 
        Wall{1} = double( Mdd );
        Wall{2} = double( getLVRReconstructedNet(node_features,TableNode, disease_table, [], knn.disease, []  , keepweight) );         
%         [Mdd] = NetworkFusion_SNF(Wall ); 
% %         [Mdd] = NetworkFusion_SNF(Wall,K,t,ALPHA); 
%         Mdd = Wall{1} + Wall{2};
%         issymmetric( Wall{1} )
%         issymmetric( Wall{2} )
        Wall{1} = Wall{1}*Wall{1}'; 
        Wall{2} = Wall{2}*Wall{2}'; 
        [Mdd] = NetworkFusion_karchermean(Wall{:} );
        
    otherwise; error('No defintion')
        
        
end
% % % % % % % % % 
switch lower(GeneNetType)
    case lower({'None','NoERC'})
        Mgg = Mat(ind_gene,ind_gene); 
        Mgg(logical(speye( size(Mgg) )) )= 1; 
%         size(Mgg  )
%         ss  = sum( diag(Mgg)   )         
%         if sum( Mgg(:)  )==0; Mgg = speye(Mgg); end 
        
    case lower('GeERC')
        if isempty( node_features ) 
            node_features = dlmread(file_emb,'',1,0) ; 
        end 
% %         knn.gene = sum(Mat(ind_gene,ind_gene),2);
        [Mgg ]=getLVRReconstructedNet(node_features,TableNode, gene_table, [], knn.gene, []  , keepweight);         
                 
    case lower('GeFN')  
        if isempty( node_features ) 
            node_features = dlmread(file_emb,'',1,0) ; 
        end 
        Mgg = Mat(ind_gene,ind_gene); 
        Mgg(logical(speye( size(Mgg) )) )= 1; 
%         if sum( Mgg(:)  )==0; Mgg = speye(Mgg); end 
        Wall{1} = double(Mgg);
% %         knn.gene = sum(Mgg,2);
        Wall{2} = double(getLVRReconstructedNet(node_features,TableNode, gene_table, [], knn.gene, []  , keepweight));           
%         [Mgg] = NetworkFusion_SNF(Wall ); 
        Wall{1} = Wall{1}*Wall{1}'; 
        Wall{2} = Wall{2}*Wall{2}'; 
        [Mgg] = NetworkFusion_karchermean(Wall{:} ); 
%         Mgg = Wall{1} + Wall{2}; 
    otherwise; error('No defintion')
       
         
end

% % % % % % % % % % % % % 
% keepweight = true; 
% % node_type1 = 'disease';  %node_type2 = 'disease'; 
% getLVRReconstructedNet(node_features,node_table, disease_table, [], knn.disease, out_file_dis_dis_net  , keepweight); 
% % % % % % % % % % % % % % % % 
% % node_type1 = 'gene';     %node_type2 = 'gene'; 
% getLVRReconstructedNet(node_features,node_table, gene_table, [], knn.gene, out_file_gene_gene_net  , keepweight ); 
%  ï¿½ï¿½È¡ train ï¿½ï¿½Ïµ 
% edgelist_dis_gene =readtable(file_dis_gene_net,'ReadVariableNames',false   , 'Format', '%s%s%s');
% ind = strcmpi( 'train' , edgelist_dis_gene{:,3});
% edgelist_dis_gene_train = edgelist_dis_gene(ind, :  ); 
% writetable(edgelist_dis_gene_train, out_file_dis_gene_train, 'Delimiter', '\t' ,'WriteVariableNames',false);
% % % % % % % % % % % % % % % % % 
% ï¿½ï¿½ï¿½ï¿½ï¿½ÛºÏµï¿½ disease-gene ï¿½ì¹¹ï¿½ï¿½ï¿½ï¿½ 
% HerNetFileInfoSet.dis_dis  ={ out_file_dis_dis_net ,  'disease', 'disease' }   ;  
% HerNetFileInfoSet.gene_gene={ out_file_gene_gene_net, 'gene',    'gene' }   ;  
% HerNetFileInfoSet.dis_gene ={ file_dis_gene_net, 'disease', 'gene', 'train' }   ;  % only use  dis-gene train
% %
% StartID = 0 ;
% SaveDir = '' ;
% IsWeighted = false; 
% [netfile, Mat, TableNode  ] = getHeterNet_tabdelimiter(HerNetFileInfoSet ,StartID, SaveDir, IsWeighted,  false) ; 
% ind_gene = find( strcmp(TableNode{:,'type'},'gene')       ) ; 
% ind_dis  = find( strcmp(TableNode{:,'type'},'disease')    ) ;   
%
% % ind_gene = TableNode{ strcmp(TableNode{:,'Var3'},'gene')   , 1}; 
% % ind_dis  = TableNode{ strcmp(TableNode{:,'Var3'},'disease')   , 1};  
% Mdd = Mat(ind_dis,ind_dis);
% if sum( Mdd(:)  )==0; Mdd = speye(Mdd); end 
Mgd = ( Mat(ind_gene,ind_dis) );   [n_gene, n_dis] = size( Mgd ) ;
% Mdd = sparse(  squareform(  1-pdist( Mgd' ,'cosine' ) )>0.8   ); 
% % Mdd = corr( Mgd )>0.8 ; 
n_dis_gene_all = nnz( sum(Mgd,2)>0 )
% Mgg = Mat(ind_gene,ind_gene);   
% TableNode.Properties.RowNames = TableNode{:,1};   % for extracting numeric ID of disease 
% disease_table = TableNode(ind_dis , :);  disease_table.IDnew = (1:size(disease_table,1))';  
% gene_table    = TableNode(ind_gene , :); gene_table.IDnew    = (1:size(gene_table,1))';
% %  class( disease_table{:,'IDnew'} )
%   È¡ test diseases    
T_dis_gene =readtable(file_dis_gene_net,'ReadVariableNames',false   , 'Format', '%s%s%s');
% diseasesIDstr = T_dis_gene{ strcmpi( 'test' , T_dis_gene{:,3}),  1  };
% diseasesIDstr = unique( diseasesIDstr );  
% % % % % % % % % % % % % 
%for AUC  
testlist                         = T_dis_gene(strcmp(T_dis_gene{:,3},'test'),:);
[disease_test_set,ia_dis,ic_dis] = unique(testlist{:,1});  testdiseasesIDstr = disease_test_set ; 
[gene_test_set,ia_gene,ic_gene]  = unique(testlist{:,2}); 
IDdisease_test = disease_table{disease_test_set,'IDnew'};    ind_dis_test  = IDdisease_test(ic_dis); 
IDgene_test    = gene_table{gene_test_set,'IDnew'};          ind_gene_test = IDgene_test(ic_gene); 
M_gene_dis_test = sparse(ind_gene_test,ind_dis_test, 1, n_gene, n_dis ) ; 
M_gene_dis_test = M_gene_dis_test(:, IDdisease_test );   
% ind_pos_pairs_geneVSdis_test = sub2ind([n_gene, n_dis ] , ind_gene, ind_dis) ; 
% % % % % 
% % % % % % % % % % % % % % % % % % % % % % % % 
% %     pai = 0.7
% %     C = 1e-10
% %     tao = 0.8
restart  = 0.7 ;
eps_c    = 1e-10; 
pro_jump = 0.8 ; 
eta      = 0.5; 
% % NormalizationType = 'col'; 
% % ids = disease_table{diseasesIDstr,2}; 
P0_G = Mgd(:, IDdisease_test  );  % for indicator of train  
if strcmpi(UseMethodType,'PRINCE')
    % % [Fset ] = A_PRINCEplus(AdjGfG,AdjGfD,AdjDfD, DisIDset, restart,LogisticFuncParaC,options, submethod, IniProExtendType  )     
    LogisticFuncParaC = -15 ; 
    options.IsTableOutput = 1; 
    % 'RWRextend'   'LPM'
    if strcmpi( NormalizationType, 'lap'  )
        submethod = 'PRINCE';    
    elseif strcmpi( NormalizationType, 'col'  )
        submethod = 'RWRextend';    
    elseif strcmpi( NormalizationType, 'row'  )
        submethod = 'LPM';   
    else; error('no definition. ')
    end        
    IniProExtendType = 'PrinceIniProExtend';    
    [P_G ] = A_PRINCEplus(Mgg,Mgd,Mdd, IDdisease_test, restart,LogisticFuncParaC,options, submethod, IniProExtendType  ) ; 
%     size(M_gene_dis_test)
%     size(P_G)
    
elseif strcmpi(UseMethodType,'CIPHER_DN')
    Type_Topological_Closeness = 'DN';
    [P_G ] = A_CIPHER(Mgg,Mgd,Mdd, IDdisease_test, 1,Type_Topological_Closeness, 0) ; 
elseif strcmpi(UseMethodType,'CIPHER_SP')
    Type_Topological_Closeness = 'SP';
    [P_G ] = A_CIPHER(Mgg,Mgd,Mdd, IDdisease_test, 1,Type_Topological_Closeness, 0) ; 
    
elseif strcmpi(UseMethodType,'BiRW')  
    P0_G_train =  Mgd ;    
% %     P0_G_train =  P0_G_train + any(Mgd,2)*0.001;    
% %     alpha = restart; 
    alpha = 0.8; 
% % %     [P_G] = A_BiRW_mn(Mgg, Mdd, Mgd, 4, 4, alpha, [], []); 
    [ P_G] = A_BiRWplus(Mgg, Mdd, P0_G_train, 4, 4, alpha, NormalizationType) ; 
    P_G = P_G(:, IDdisease_test ); 
    
elseif strcmpi(UseMethodType,'KATZ')  
    P0_G_train =  Mgd ;    
% %     P0_G_train =  P0_G_train + any(Mgd,2)*0.001;     
    P_G = A_KATZ_H_ALL_DGP(Mgg,{P0_G_train}, Mdd, [] );  
    P_G = P_G(:, IDdisease_test ); 
    
elseif strcmpi(UseMethodType,'RWRH')
    P0_G_train =  Mgd(:, IDdisease_test  ) ;    
    % %     P0_G_train =  P0_G_train + any(Mgd,2)*0.001;
    % %     P0_G_train =  P0_G_train + sum(Mgd,2)*0.001;
    % % fd = sum(Mgd,2)./size(Mgd,2) 
    % % fg = sum(Mgd,1)./size(Mgd,1);
    % % P0_G_train =  P0_G_train + fd;
    P0_D = eye( size(Mdd) ); P0_D = P0_D(:, IDdisease_test  );   
    % % [P_G, P_D ] = A_RWR_Hplus(Mgg,Mgd,Mdd, P0_G,P0_D, restart, pro_jump, eta, NormalizationType, []) ;
    if strcmpi(NormalizationType,'lap'); NormalizationType = 'LaplacianNormalization' ;  end  
    [P_G ] = A_RWR_Hplus(Mgg,Mgd,Mdd, P0_G_train,P0_D, restart, pro_jump, eta, NormalizationType, []) ;

elseif strcmpi(UseMethodType,'RWRHini')
    P0_G_train =  Mgd ;    
    fd = sum(P0_G_train,2)./size(P0_G_train,2) ;
    fg = sum(P0_G_train,1)./size(P0_G_train,1);
    P0_G_train =  P0_G_train + (fg+fd)/2;
    P0_G_train =  P0_G_train(:, IDdisease_test  ) ;    
% %     P0_G_train =  P0_G_train + any(Mgd,2)*0.001;
% %     P0_G_train =  P0_G_train + sum(Mgd,2)*0.001;
% % P0_G_train =  P0_G_train + any(Mgd,2)*0.001;

    P0_D = eye( size(Mdd) ); P0_D = P0_D(:, IDdisease_test  );   
    % % [P_G, P_D ] = A_RWR_Hplus(Mgg,Mgd,Mdd, P0_G,P0_D, restart, pro_jump, eta, NormalizationType, []) ;
    if strcmpi(NormalizationType,'lap'); NormalizationType = 'LaplacianNormalization' ;  end  
    [P_G ] = A_RWR_Hplus(Mgg,Mgd,Mdd, P0_G_train,P0_D, restart, pro_jump, eta, NormalizationType, []) ;
    
    
elseif strcmpi(UseMethodType,'RWR')  
    P0_G_train =  Mgd(:, IDdisease_test  ) ;
% %     P0_G_train =  P0_G_train + any(Mgd,2)*0.001;
    [P_G ]= A_RWRplus(Mgg, restart, P0_G_train, [], [], [],  NormalizationType); 
    
elseif strcmpi(UseMethodType,'DK')  
    P0_G_train =  Mgd(:, IDdisease_test  ) ;
% %     P0_G_train =  P0_G_train + any(Mgd,2)*0.001; 
    AdjSet.PPI = Mgg; 
    % %     [~, P_G, TableScores ] = A_DiffusionKernel_M(AdjSet, P0_G_train, p_beta, pro_jump, UsedGlobalVar, options, NormalizationType); 
    [~, Pm ] = A_DiffusionKernel_M(AdjSet, P0_G_train, [], 0.5, [], [], []); 
    P_G = Pm.mean ; 
    
elseif strcmpi(UseMethodType,'GUILD')  
    P0_G_train =  Mgd(:, IDdisease_test  ) ;
    methodname = 'NetCombo'; 
% %     P0_G_train =  P0_G_train + any(Mgd,2)*0.001;
%     [P_G ]= A_RWRplus(Mgg, restart, P0_G_train, [], [], [],  NormalizationType); 
    [P_G ] = A_DGP_GUILD_singlenet(Mgg, P0_G_train, methodname, [], []) ; 

else
    error('no defintion')
end
% % % P_G = P_G
P_G(  P0_G~=0  ) = -inf  ; 
scores_gene_dis  = P_G ; 
% save top k to file
disease_gene_tableset = []; 
if ~isempty(file_pred_results) 
    [scores,II] = sort( P_G ,1, 'descend' ); 
    topk_save = min(n_gene,topk_save);
    II      = II(1:topk_save, :) ; 
    scores  = scores(1:topk_save, :) ; 
    n_test_dis = length( testdiseasesIDstr );   disease_gene_tableset =repmat({[]},n_test_dis , 1 ); 
    for ii=1:n_test_dis 
        gene_topk = gene_table{ II(:,ii) , 1};  
        disease_gene_tableset{ii} = table( repmat(testdiseasesIDstr(ii),topk_save,1) ,gene_topk, (scores(:,ii)  ) ) ; 
    end
    disease_gene_tableset = cat(1, disease_gene_tableset{:} ); 
    disease_gene_tableset.Properties.VariableNames = {'disease','gene','score' } ; 
	disease_gene_tableset.Properties.VariableNames 
    % disease_gene_tableset.Properties.VariableNames = {'disease','gene','score','rank'}; 
    writetable(disease_gene_tableset, file_pred_results , 'Delimiter', '\t'  ,'WriteVariableNames',true);
end
% toc 
% save  ttttttttttttttt.mat
%% calculate AUC ....  % % % % % % % % % % % % % % % % % % % % % % %  
if ~isempty(topk_perf_set) && ~isempty(M_gene_dis_test)
    [ResScalartable, ResCurvetable ,ResTopKtable, n_gene,n_dis_test] = getResPerfFromMatrixGeneDis(M_gene_dis_test,scores_gene_dis, topk_perf_set, topk_save) ;  n_dis
else
    ResScalartable =table;  
    ResCurvetable =table; 
end

end
