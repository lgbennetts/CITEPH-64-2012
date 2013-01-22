function freq = FindFreq_FS(parameter_vector, kk, HH)

% freq = FindFreq_FS(parameter_vector, kk, HH)
%
% - function that produces the frequency from a given waveno in the
% free-surface region - %

freq = kk.*tanh(kk.*HH);

freq = sqrt(parameter_vector(5)*freq);

return