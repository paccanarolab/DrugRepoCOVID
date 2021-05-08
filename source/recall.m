    function [maxrecallPerDrug, maxrecallPerVirus, paramPerDrug, paramPerVirus] = recall(  method, topN, type)
    %UNTITLED2 Summary of this function goes here
    %   Detailed explanatio
       
    load(strcat('output\', method,type, '.mat'), 'rankPerDrug', 'rankPerVirus', 'typeEntry', 'allparams');
    
    
    unique_params = unique(allparams(2:end,:)', 'rows');


    recallDrugs = zeros(length(unique_params),length(topN));
    recallViruses = zeros(length(unique_params),length(topN));

   
    for i = 1:size(unique_params, 1)
        % index for a configuration of parameters
        idx = ones(1, size(allparams, 2)) > 0;        
        
        if size(unique_params, 2) > 1
            for j = 2: size(allparams, 1)
                idx = idx & allparams(j, :) == unique_params(i, j-1);
            end
        else
            idx = allparams(2, :) == unique_params(i);
        end
        for j = 1:length(topN)
            recallDrugs(i,j) = sum(rankPerDrug(idx) <= topN(j)) / sum(idx);
            recallViruses(i,j) = sum(rankPerVirus(idx) <= topN(j)) / sum(idx);
        end
    end

    % max value
    [~, Idx] = max(mean(recallDrugs, 2));
    [~, Idy] = max(mean(recallViruses, 2));
    
    maxrecallPerDrug = recallDrugs(Idx, :);
    maxrecallPerVirus = recallViruses(Idy, :);
    
    paramPerDrug = unique_params(Idx,:);
    paramPerVirus = unique_params(Idy,:);
    
 end

