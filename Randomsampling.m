% Inputs:
% c: the number of samples to take
% pk: a vector of length n representing the probability mass function of the rows of A and columns of B
% A: matrix A
% B: matrix B
% Outputs:
% E: the relative error of the approximation
% M: an approximation of the product AB, computed as a sum of one-rank products


function M = Randomsampling(T,G,c)
    n=size(T,1);
    
    rowsums = sum(T,2); % compute the row sums of A
    colsums = sum(G,1); % compute the column sums of B
   
    pk = ((rowsums/n))*(colsums/n); % define pk as the outer product of the row and column sums, normalized by n
    
    pk = pk/sum(pk); % normalize pk so that it sums to 1
% it assumes that all the rows of T and columns of G are equally important
    M=zeros(size(T,1),size(G,2));
    for t = 1:c
        k=find(rand(1)<cumsum(pk),1,"first");

        one_rank_prod = T(k,:)'*G(:,k)'/(c*pk(k));
        M = M + one_rank_prod;
    end


