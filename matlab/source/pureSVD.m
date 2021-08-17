function [ R ] = pureSVD( Y, f )
%pureSVD truncated singular value decomposition
[m,n] = size(Y);
 %[Ynorm, Ymean] = normalizeRatings(Y);
if f <= m
    [U,S,Q] = svd(Y);
    U = U(:,1:f);
    Sigma = S(1:f,1:f);
    Q = Q(:,1:f);
    R = U*Sigma*Q';
else
    error(' The truncated factor f, cannot be bigger than the number of users'); 
end
 %R = R + repmat(Ymean, 1, n);

end

