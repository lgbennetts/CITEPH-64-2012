function s = stats(u)
%
%  s = stats(u)
%
%  Returns s = [mean, standard deviation, skewness, kurtosis] of u
%  Note: if u(samples,cases) is a matrix then s takes the form s(cases,4)


[ns,nc] = size( u ) ;                            % [#_samples, #_cases]

Is = ones(ns,1) ;

mean_u = mean(u) ;                               % mean

std_u  = std(u) ;                                % std deviation

skew_u = mean( (u-Is*mean_u).^3 ) ./ std_u.^3;   % skewness

kurt_u = mean( (u-Is*mean_u).^4 ) ./ std_u.^4;   % kurtosis

s = [mean_u' std_u' skew_u' kurt_u'] ;
