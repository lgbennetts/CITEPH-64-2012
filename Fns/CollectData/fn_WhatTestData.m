function fn_WhatTestData

cprintf(0.6*[1,1,1],'Hs [mm], Tp [s]\n')
cprintf(0.6*[1,1,1],'<.79 conc>:\n')

c79_prams = conc79_testspecs();

N = length(c79_prams);

for loop=1:N
 dum(loop,1) = 10*c79_prams(loop).wave_height;
 dum(loop,2) = 10\c79_prams(loop).period;
end

disp(dum)
clear dum

cprintf(0.6*[1,1,1],'<.39 conc>:\n')

c39_prams = conc39_testspecs();

N = length(c39_prams);

for loop=1:N
 dum(loop,1) = 10*c39_prams(loop).wave_height;
 dum(loop,2) = 10\c39_prams(loop).period;
end

disp(dum)

return