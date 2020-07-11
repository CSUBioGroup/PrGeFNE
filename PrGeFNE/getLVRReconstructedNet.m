function [SimMatrix, node_table1, node_table2 ]=getLVRReconstructedNet(node_features,node_table, node_type1, node_type2, knn, out_file ,keepweight )
    if isa(node_type1,'char') 
        node_table1 = node_table(strcmpi(node_type1, node_table{:,3}),  : ); 
    elseif isa(node_type1,'table') 
        %node_type1 is a sub-node_table rcording all information of subset of nodes
        node_table1 = node_type1;  
    else 
        error( 'No defintion' )
    end
    if  isa(node_type2,'char') 
        node_table2 = node_table(strcmpi(node_type2, node_table{:,3}),  : );
    elseif isa(node_type2,'table')
        %node_type1 is a sub-node_table rcording all information of subset of nodes 
        node_table2 = node_type2; 
    elseif isempty( node_type2  )
        node_table2 = []; 
    else 
        error( 'No defintion' )
    end
    % % % % % % % % % % % % % % % % % 
    SimMatrix=getLVRSim(node_features, node_table1, node_table2  ); 
    if exist('knn','var') && ~isempty(knn)
        %  sim = Adjknn = getAdjKnnColumns( SimMatrix,  k_neighbors_vec , symmetrized, keepweight ); 
        k_neighbors_vec = knn ; 
        symmetrized     = true ; 
        % keepweight      = true ; 
        SimMatrix = getAdjKnnColumns( SimMatrix,  k_neighbors_vec , symmetrized, keepweight ); 
    end    
    if exist('out_file','var') && ~isempty(out_file)
        [i,j,v] =find( SimMatrix );
% %         listtable = table(i,j,v);
        if isempty( node_type2  ); node_table2 = node_table1;  end    
        listtable = table(node_table1{i,1},  node_table2{j,1},  v);
        writetable(listtable, out_file , 'Delimiter', '\t'  ,'WriteVariableNames',false);
    end    
end
