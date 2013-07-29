%% C112083_Wave_calibration+Irregular_Spectrum.m
%% Author: Timothy Williams
%% Date:   20130724, 09:39:18 CEST
%% from C112083_Wave_calibration+Irregular_Spectrum.xlsx
function c79_prams = conc79_testspecs();

names = {'6'; 
	 '7'; 
	 '8'; 
	 '9'; 
	 '10'; 
	 '11'; 
	 '11a'; 
	 '11b'; 
	 '12'; 
	 '13'; 
	 '14'; 
	 '15'; 
	 '16'; 
	 '17'; 
	 '18'; 
	 '19'; 
	 '20'; 
	 '21'; 
	 '22'};

types = {'Regular';	    
	 'Regular';   		    
	 'Regular';   		    
	 'Regular';   		    
	 'Regular';   		    
	 'Regular';   		    
	 'Regular';   		    
	 'Regular';   		    
	 'Regular';   		    
	 'Regular';   		    
	 'Regular';   		    
	 'Regular';   		    
	 'Regular';   		    
	 'Regular';   		    
	 'Regular';   		    
	 'Irregular';   		    
	 'Irregular';   		    
	 'Irregular';   		    
	 'Irregular'}; 

periods	 = {6.5;	
	    8.0;  	
	    8.0;  	
	    9.5;  	
	    11.0;  	
	    12.5;  	
	    12.5;  	
	    12.5;  	
	    14.0;  	
	    14.0;  	
	    15.5;  	
	    17.0;  	
	    18.5;  	
	    20.0;  	
	    20.0;  	
	    8.0;  	
	    14.0;
	    14.0;
	    20.0}; 	

wave_heights   = {2;
		  2;
		  4;
		  3;
		  4;
		  4;
		  4;
		  4;
		  4;
		  8;
		  4;
		  4;
		  4;
		  4;
		  10;
		  2;
		  4;
		  8;
		  4};

dirnames    = {'23071525.a13';
	       '23071605.a13';
	       '23071625.a13';
	       '24071158.a13';
	       '24071121.a13';
	       '24071018.a13';
	       '24071037.a13';
	       '24071056.a13';
	       '24070958.a13';
	       '24071410.a13';
	       '24070941.a13';
	       '24070926.a13';
	       '23071748.a13';
	       '23071730.a13';
	       '24071427.a13';
	       '24071656.a13';
	       '24071611.a13';
	       '24071520.a13';
	       '24071452.a13'};

c79_prams   = struct('dirname',dirnames,...
		     'wave_height',wave_heights,...
		     'period',periods,...
		     'name',names,...
		     'type',types);
