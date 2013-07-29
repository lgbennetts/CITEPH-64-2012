%% C112083_Wave_calibration+Irregular_Spectrum.m
%% Author: Timothy Williams
%% Date:   20130724, 09:39:18 CEST
%% from C112083_Wave_calibration+Irregular_Spectrum.xlsx
function calib_prams = C112083_Wave_calibration_inc_Irregular_Spectrum();

names = {'Cal1'; 
	 'Cal2'; 
	 'Cal3'; 
	 'Cal4'; 
	 'Cal5'; 
	 'Cal6'; 
	 'Cal7'; 
	 'Cal8'; 
	 'Cal9'; 
	 'Cal10';
	 'Cal11';
	 'Cal12';
	 'Cal13';
	 'Cal14';
	 'Cal15';
	 'Cal16';
	 'Cal17'};

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
	 'Irregular : Joswap (3.3)'; 
	 'Irregular : Joswap (3.3)'; 
	 'Irregular : Joswap (3.3)'; 
	 'Irregular : Joswap (3.3)'};

periods	 = {6.5;	
	    8.0;  	
	    8.0;  	
	    9.5;  	
	    11.0; 	
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

wave_heights   = {2.0;
		  2.0;
		  4.0;
		  3.0;
		  4.0;
		  4.0;
		  4.0;
		  8.0;
		  4.0;
		  4.0;
		  4.0;
		  4.0;
		  10.0;
		  2.0;
		  4.0;
		  8.0;
		  4.0};

dirnames    = {'18071805.a13';
	       '18071410.a13';
	       '18071719.a13';
	       '18071644.a13';
	       '18071706.a13';
	       '18071425.a13';
	       '18071458.a13';
	       '18071731.a13';
	       '18071523.a13';
	       '18071548.a13';
	       '18071617.a13';
	       '18070905.a13';
	       '18070828.a13';
	       '19071216.a13';
	       '19071109.a13';
	       '19070915.a13';
	       '19070851.a13'};

calib_prams = struct('dirname',dirnames,...
		     'wave_height',wave_heights,...
		     'period',periods,...
		     'name',names,...
		     'type',types);
