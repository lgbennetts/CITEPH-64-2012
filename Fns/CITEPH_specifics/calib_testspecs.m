% function calib_prams = calib_testpecs();
%
% L Bennetts Oct 2013 / Adelaide


function calib_prams = calib_testpecs();

names = {'1'; 
	 '2'; 
	 '3'; 
	 '4'; 
	 '5'; 
	 '6'; 
	 '7'; 
	 '8'; 
	 '9'; 
	 '10';
     '11';
     '12';
     '13';
     '14';
     '15';
     '16';
     '17'}; 

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
	 'Irregular: Jonswap (3.3)';   
	 'Irregular: Jonswap (3.3)';   		    
	 'Irregular: Jonswap (3.3)'   
	 'Irregular: Jonswap (3.3)';}; 








periods	 = {6.5;
        8.0;
        8.0;
	    9.5;  
        11.0;
        12.5;  	
        14.0;
        14.0
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



calib_prams   = struct('dirname',dirnames,...
		     'wave_height',wave_heights,...
		     'period',periods,...
		     'name',names,...
		     'type',types);

return
