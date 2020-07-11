% % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % 
function Adjknn = getAdjKnnColumns( SimMatrix,  k_neighbors_vec , symmetrized, keepweight )
% Input: 
% SimMatrix  similarity matrix
% k_neighbors  vector with number of neighbors of each node, default: 5 
% symmetrized  1 or 0
% keepweight   keep similarity as weight of edges
% Output: Adjknn   matrix
% 从相似性矩阵提取  k近邻矩阵
% Ju Xiang 
% 2019-5
    if isempty(SimMatrix) 
        % sort_dim = 2;
        SimMatrix =rand(10);   % for testing only 
        warning('test');
    end 
    SZ =size( SimMatrix); 
    % 
    if isscalar(k_neighbors_vec)
        k_neighbors_vec = k_neighbors_vec(1)*ones( SZ(1), 1 ); 
    end
    if isempty(symmetrized) 
        symmetrized = true;
    end
    if isempty(keepweight) 
        keepweight = false;
    end
     
    if any( k_neighbors_vec>SZ(1)-1 )
        k_neighbors_vec(  k_neighbors_vec>SZ(1)   ) = SZ(1)-1;
        warning( ['There is k_neighbors:','>', num2str(SZ(1)-1),' the maximal number of neighbors'] );
    end
    
    SimMatrix = SimMatrix + rand( SZ ).*( min(abs(SimMatrix(:)))/100000000 );  % 添加随机扰动  
    SimMatrix(   sub2ind( SZ, 1:SZ(1),1:SZ(1) )      ) = -inf;   % 使 自连接排在最末尾 
    % SimMatrix(   ( eye( SZ ) )==1       ) = -inf; 
    [~,II] = sort( SimMatrix ,2, 'descend' );  
    Adjknn = zeros( SZ );
    for ii=1: SZ(1)
        knn = II(ii,1: k_neighbors_vec(ii) ); 
        if keepweight
            Adjknn(ii,  knn ) = SimMatrix(ii,  knn );    
        else
            Adjknn(ii,  knn ) = 1; 
        end
    end 
    
    if symmetrized
        [i,j,v] = find( Adjknn ); 
        ind = sub2ind( SZ, i ,j );
        Adjknn = Adjknn' ; 
        Adjknn(ind) = v; 
        % %     Adjknn = Adjknn';
        % %     Adjknn(ind) = 
        % Adjknn = full
    end 
end


