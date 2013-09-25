%% C112083_Wave_calibration+Irregular_Spectrum.m
%% Author: Timothy Williams
%% Date:   20130724, 09:39:18 CEST
%% from C112083_Wave_calibration+Irregular_Spectrum.xlsx
function c39_prams = conc39_testpecs();





names = {'1'; 
	 '2'; 
	 '2a'; 
	 '2b'; 
	 '3'; 
	 '4'; 
	 '5'; 
	 '6'; 
	 '7'; 
	 '8'; 
	 '9'; 
	 '10'}; 








types = {'Regular';	    
	 'Regular';   		    
	 'Regular';   		    
	 'Regular';   		    
	 'Regular';   		    
	 'Regular';   		    
	 'Regular';   		    
	 'Regular';   		    
	 'Regular';   		    
	 'Irregular: Jonswap (3.3)';   
	 'Irregular: Jonswap (3.3)';   		    
	 'Irregular: Jonswap (3.3)'}; 








periods	 = {6.5;	
	    9.5;  	
	    9.5;  	
	    9.5;  	
	    12.5;  	
	    15.5;  	
	    18.5;  	
	    12.5;  	
	    9.5;
	    14.0;
	    8.0;  	
	    14.0}; 	








wave_heights   = {2;
                  3;
		  3;
                  3;
                  4;
                  4;
                  4;
                  8;
                  6;
                  8;
                  2;
                  4};








dirnames    = {'25071204.a13';
               '25071232.a13';
               '25071513.a13';
               '25071532.a13';
               '25071413.a13';
               '25071432.a13';
               '25071453.a13';
               '25071558.a13';
               '25071620.a13';
               '25071643.a13';
               '25071717.a13';
               '25071743.a13'};



c39_prams   = struct('dirname',dirnames,...
		     'wave_height',wave_heights,...
		     'period',periods,...
		     'name',names,...
		     'type',types);
