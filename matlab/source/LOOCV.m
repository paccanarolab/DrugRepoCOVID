function [] = LOOCV(X, rows, cols, method, input_representation, f, alpha_zeros, alpha_invitro, typeevaluation)
%UNTITLED Summary of this function goes here
%   Detailed explanation goes here
[~, nvirus] = size(X);


fprintf('\n number of associations for the LOOCV %d\n', length(rows));

% model parameters
n = 1:length(rows);
alpha = 0; masks = 0;
if strcmp(method, 'cNMF')    
     allparams = combvec(n, f, alpha_zeros, alpha_invitro);  
     alpha = [0.1, 0.1, 0.1, 0.1, 0.16, 0.27, 0.73, 1, 1];
        % create the masks
        masks = cell(9,1);
        for i = 1:9    
            masks{i} = X == i-1;
        end

      size(allparams, 2)
else
    allparams = combvec(n, f);
end


if strcmp(input_representation, 'binary')
    Xcomplete = X;
    X = double(X > 0);
end

rankPerDrug = zeros(size(allparams, 2),1);
rankPerVirus = zeros(size(allparams, 2), 1);
%scoresPerVirus = cell(size(allparams, 2), 1);
%labelTestSet = cell(size(allparams, 2), 1);
typeEntry = zeros(size(allparams, 2), 1);

maxiter = 2000;
m = size(allparams, 2);
fprintf('\t Completion: ');
showTimeToCompletion; startTime=tic;
p = parfor_progress( m );

parfor i = 1:m
        Xtrain = X;
        idxEntry = allparams(1,i);
        Xtrain(rows(idxEntry), cols(idxEntry)) = 0;
        
        typeEntry(i) = Xcomplete(rows(idxEntry), cols(idxEntry));
        
        % model
        nofeatures = allparams(2,i);
           
        if strcmp(method, 'tSVD')

            [ Xhat ] = pureSVD( Xtrain, nofeatures );
        
        elseif strcmp(method, 'NMF')

            [W,H] = nnmf(Xtrain, nofeatures);
             Xhat = W*H;
        
        elseif strcmp(method, 'cNMF')

            myalpha = alpha;

            myalpha(1) = allparams(3,i);
            myalpha(2:4) = allparams(4,i); 

            [ W, H, ~ ] = cNMF( Xtrain, nofeatures, myalpha, masks, maxiter );
             Xhat = W*H;
        
        elseif strcmp(method, 'random')
            
             Xhat = rand(size(Xtrain));
            
        end
        
        % save ranking per drug
        idx = Xtrain(rows(idxEntry), :) > 0; % get those in training.
        scores = Xhat(rows(idxEntry),:); % predicted scores
        scores(idx) = -1; % remove training
        [~,I] = sort(scores, 'descend');
        
        rankPerDrug(i) = find(I == cols(idxEntry));  
        
                
        % save ranking per virus
        idx = Xtrain(:, cols(idxEntry)) > 0; 
        scores = Xhat(:, cols(idxEntry));
        scores(idx) = -1; % remove training
        [~,I] = sort(scores, 'descend');
        
        rankPerVirus(i) = find(I == rows(idxEntry));  
        %scoresPerVirus{i} = scores;  
        % parfor progress 
         p = parfor_progress;
         showTimeToCompletion( p/100, [], [], startTime );
        
end

    
save(strcat('output\', method, typeevaluation, '.mat'), 'rankPerDrug', 'rankPerVirus', 'typeEntry', 'allparams');
end

