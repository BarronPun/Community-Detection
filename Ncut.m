function labels = Ncut (Mat)  % return the labels of each nodes
    M = Mat;
    num = size(M,1); %rows and cols of M, number of nodes
    D = zeros(num,num);%Diagonal matrix corresponding to M
    % initialize D, Dii is the sum of the ith row's elements
    for i = 1:num
       temp = 0;
       for j = 1:num % sum of ith row's elements
          temp = temp + M(i,j); 
       end
       D(i,i) = temp;
    end
    
    %targetM is the Matrix whose eigenvalue and eigenvector are going to be
    %calculated
    targetM = D^(-1/2)*(D-M)*D^(-1/2);
    
    %calculate the eigenvalues and the corresponding eigenvectors
    %Stored in eigvalueM and eigvectorM Respectively
    %'doc eig' for more information in command window
    [eigvectorM, eigvalueM] = eig(targetM);
    
    %Extract the eigenvalues from eigvalueM, store it in eigvalues
    eigvalues = zeros(1,num);
    for i = 1:num
        eigvalues(i) = eigvalueM(i,i);
    end
    
    %sort the eigvalues in increasing order
    %also sort the indexes of elements of eigvalues
    indexes = zeros(1,num);
    for i = 1:num
        indexes(i) = i;
    end
    for i = 1:num
       for j = i+1:num
           if eigvalues(i) > eigvalues(j) %swap eigvalues and indexes
               temp = eigvalues(i);
               eigvalues(i) = eigvalues(j);
               eigvalues(j) = temp;
               
               temp = indexes(i);
               indexes(i) = indexes(j);
               indexes(j) = temp;
           end
       end
    end
    
    %Get the top n small eigenvalues' corresponding eigenvectors
    %Construct the dataset 'embedM' with the eigenvectors
    n = 15;
    embedM = zeros(num,n);
    for i = 1:n
        embedM(:,i) = eigvectorM(:,indexes(i));
    end
    
    %Use k-means to classify the the data in embedM
    %each row of embedM can represent a data point
    %choose pre_k = 10
    pre_k = 10;
    [~,result] = K_means(embedM,pre_k);
    prelabels = result(:,n+1);  % prelabels is the labels after k-means
    
    %------------%
    %Merge two small community in greedy manner at a time
    %until k cummunity left, here k = 5
    k = 5;
    while pre_k > k
        maxlabel = max(prelabels); %max label is also the number of types
        ncut = 0;
        choicei = 1;
        choicej = 2;
        for i = 1:maxlabel
           for j = i+1:maxlabel
               if i==1&&j==2 %initial, merge 1 and 2, and get the ncut
                   templabels = prelabels; %temporary labels
                   for x = 1:num %merge
                      if templabels(x) == 2
                          templabels(x) = 1;
                      elseif templabels(x) > j
                          templabels(x) = templabels(x) - 1;
                      end 
                   end
                   ncut = currentNcut(M,templabels);
                   choicei = i;
                   choicej = j;
               else  % merge other pairs of nodes
                   templabels = prelabels;
                   for x = 1:num %merge
                      if templabels(x) == j
                         templabels(x) = i;
                      elseif templabels(x) > j
                         templabels(x) = templabels(x) - 1;
                      end
                   end
                   temp = currentNcut(M, templabels);
                   if temp < ncut % smaller than the current ncut, update
                      ncut = temp;
                      choicei = i;
                      choicej = j;
                   end
               end
           end
        end
        % We have got the best choice of i and j, choicei and choicej
        % update prelabels according to choicei and choicej
        for i = 1:num
           if prelabels(i) == choicej
               prelabels(i) = choicei;
           elseif prelabels(i) > choicej
               prelabels(i) = prelabels(i) - 1;
           end
        end
        %After each merge, pre_k should be subtracted by 1
        pre_k = pre_k - 1;
    end
    
    labels = prelabels;  % return labels
    
end

function ncut = currentNcut(V, labels)
    % The smallest label here is 1
    % So we should find the maximum label in lables
    % The maximum label is also the number of the label types
    maxlabel = max(labels);
    len = length(labels); %length of labels
    num = size(V,1);     %cols or rows of V, number of nodes
    %ncut = cut(A1, V-A1)/assoc(A1, V)+...+cut(Ak, A-Ak)/assoc(Ak, V)
    ncut = 0;
    for i = 1:maxlabel 
        cut = 0;
        assoc = 0;
        index1 = 1;
        index2 = 1;
        for j = 1:len
            if labels(j) == i %if equal the current label
                sub_index(index1) = j; % record this node
                index1 = index1 + 1; 
            else %if not, it belongs to the complementary sets
                com_of_sub(index2) = j; %com_of_sub is the complementary set
                index2 = index2 + 1;
            end
        end
        sublen = length(sub_index);
        comlen = length(com_of_sub);
        
        %calculate the value of cut(Ai,V-Ai)
        for j = 1:sublen
           for k = 1:comlen
               cut = cut + V(sub_index(j),com_of_sub(k));
           end
        end
        
        %calculate the value assoc(Ai, V)
        for j = 1:sublen
           for k = 1:num
              assoc = assoc + V(sub_index(j),k); 
           end
        end
        
        %calculate the ith cut/assoc and add it to ncut
        ncut = ncut + cut/assoc;
    end

end