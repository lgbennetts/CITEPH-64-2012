% function prams = singlefloe_testspecs()
% 
% test parameters and file names for the single floe tests
% 
% L Bennetts Aug 2013 / Adelaide

function prams = singlefloe_testspecs();

names = {'37'; 
	 '38'; 
	 '39'; 
	 '40'; 
	 '41'; 
	 '42'; 
	 '43' };

types = {'Regular';	    
	 'Regular';   		    
	 'Regular';   		    
	 'Regular';   		    
	 'Regular';   		    
	 'Regular';   		    
	 'Regular' }; 

periods	 = {18.5;	
	    15.5;  	
	    12.5;  	
	    9.5;  	
	    6.5;  	
	    12.5;  	
	    9.5  }; 	

wave_heights   = {4;
		  4;
		  4;
		  3;
		  2;
		  8;
		  6};

dirnames    = {'02081148.a13';
	       '02081217.a13';
	       '02081407.a13';
	       '02081428.a13';
	       '02081503.a13';
	       '02081523.a13';
	       '02081543.a13' };

prams   = struct('dirname',dirnames,...
		     'wave_height',wave_heights,...
		     'period',periods,...
		     'name',names,...
		     'type',types);

return
            