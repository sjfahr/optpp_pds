cd /FUS4/data2/sjfahrenholtz/MATLAB/Tests/direct_search
load GPU_global_search.mat

cd /FUS4/data2/sjfahrenholtz/gitMATLAB/opt_new_database/PlanningValidation
exclude_index_lo_best = find( (datasummary(:,4)<700)==1);
exclude_index_hi_best = find( (datasummary(:,4)>700)==1);

% Best dice above 700 when the global best dice is below 700
total_hi = cell ( (length(total) - length(exclude_index_lo_best)),5);
index_hi = zeros( (length(total) - length(exclude_index_lo_best)),1);
val_hi = index_hi;
mu_blip = index_hi;
diff_hi = zeros ( (length(total) - length(exclude_index_lo_best)),2);
clear ii jj kk
for ii = 1:length(exclude_index_lo_best)
    jj = exclude_index_lo_best(ii);
    
    for kk = 1:5
        if kk == 2  || kk == 3
            total_hi{ii,kk} = total{jj,kk}(702:end,:);
        else

            total_hi{ii,kk} = total{jj,kk};
        
        end
    end
    
    [val_hi(ii) index_hi(ii)] = max(total_hi{ii,3}(:,7)); % Best dice above 700 when the global best dice is below 700
    mu_blip(ii) = total_hi{ii,2}(index_hi(ii),1);
    diff_hi(ii,:) = [ (val_hi(ii) - total{jj,5}(1))  (mu_blip(ii) - total{jj,5}(2))];
end

% Best dice below 700 when the global best dice is above 700
total_lo = cell ( (length(total) - length(exclude_index_hi_best)),5);
index_lo = zeros( (length(total) - length(exclude_index_hi_best)),1);
val_lo = index_lo;
mu_blip = index_lo;
clear ii jj kk
for ii = 1:length(exclude_index_hi_best)
    jj = exclude_index_hi_best(ii);
    
    for kk = 1:5
        if kk == 2 || kk == 3
            total_lo{ii,kk} = total{jj,kk}(1:701,:);
        else
            total_lo{ii,kk} = total{jj,kk};
        
        end
    end
    
    [val_lo(ii) index_lo(ii)] = max(total_lo{ii,3}(:,7)); % Best dice below 700 when the global best dice is above 700
    mu_blip(ii) = total_lo{ii,2}(index_lo(ii),1);
    diff_lo(ii,:) = [ (val_lo(ii) - total{jj,5}(1))  (mu_blip(ii) - total{jj,5}(2))];
end

% Difference
