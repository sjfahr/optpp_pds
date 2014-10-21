% This function quickly does descriptive statistics for a 1-D distribution

function [stats] = Descriptive_statistics_LOOCV ( dist );

stats.mean = mean(dist);
stats.median = median(dist);
stats.std = std(dist);
stats.max = max(dist);
stats.min = min(dist);
stats.skew = skewness(dist);
stats.kurt = kurtosis(dist);
stats.n = length(dist);
stats.pass = sum (dist > 0.7);

end