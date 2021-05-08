function [ W, H, J ] = cNMF( X, k, alpha, masks, maxiter )
%% Non-negative matrix decomposition algorithm
%  Inputs;
%     X: drug-virus matrix of n x m.
%     k: number of latent features.
%     masks: cellarray containing matrices masking the entries of X    
%     maxiter: maximun number of iterations 
%  Outputs:
%     W: nxk matrix (drug features or drug signatures)
%     H: kxm matrix (virus features or virus signatures)
% Copyright: Diego Galeano, 2020.

            % Termination tolerance on change in the elements of W and H.
            % Default is 1e-3.
            tolx = 1e-3;   
            
            % variance
            variance = 0.01;
            
            [nrows, rcols] = size(X);
       
             % Initialization
             Wini = rand(nrows, k)*sqrt(variance);      
             Hini = rand(k, rcols)*sqrt(variance);        
             
             % normalization
             Hini = Hini./repmat(sqrt(sum(Hini.^2,2)),1,rcols); 
             
             % epsilon based on machine precision. 
             sqrteps = sqrt(eps);
            
             % cost function
             J = zeros(maxiter, 1);
                    
             for iter = 1:maxiter
                 numer = 0; denom = 0;
                 for j = 1:length(masks)
                     numer = numer + alpha(j).*(masks{j}.*X);
                     denom = denom + alpha(j).*(masks{j}.*(Wini*Hini));
                 end
                 numer = numer*Hini';            
                 denom = denom*Hini' +  eps(numer); 
                 
                                 
                 W = max(0,Wini .* (numer./denom)); 
                 
                 numer = 0; denom = 0;
                 for j = 1:length(masks)
                     numer = numer + alpha(j).*(masks{j}.*X);
                     denom = denom + alpha(j).*(masks{j}.*(W*Hini));
                 end
                 numer = W'*numer;            
                 denom = W'*denom +  eps(numer); 
                 
                 
                 H = max(0,Hini.* (numer./ denom));
                
                 
                % Compute cost function
                for j = 1:length(masks)
                      J(iter) = J(iter) + 0.5*alpha(j)*norm(masks{j}.*(X - W*H), 'fro')^2;
         
                end
                
                
                 % Get norm of difference and max change in factors                
                 dw = max(max(abs(W-Wini) / (sqrteps+max(max(abs(Wini))))));
                 dh = max(max(abs(H-Hini) / (sqrteps+max(max(abs(Hini))))));
                 delta = max(dw, dh);
                 
                 % Check for convergence
                
                 if iter > 1                      
                     if delta <= tolx               
                          
                       % fprintf('\n iter %d delta %e\n', iter, delta);
                        J(iter+1:end) = [];
                        break;               
                    
                         
                     end
                 end
                
                 % Remember previous iteration results              
                 Wini = W;       
                 Hini = H./repmat(sqrt(sum(H.^2,2)),1,rcols); % normalize
             end
         
             
            

end
