function E=spec2d_exp(data,ns,plt)
% % Matlab program to determine the directional spectrum
% % from an array of wave gauges
% %
% % The Maximum Liklehood Method is used.
% %
% % ? This version automatically loops through a group of files
% 
% % Clear old variables
%   clear
%   close
% 
% % Loading data
%   
% fname = 'c010055.no7';
%   %fname = input('the filename is ');
%   %eval(['load e:\card\',fname]);
% %   eval(['load d:\data\matlab\lg\card\',fname]);
% eval(['load ',fname]);
%   N = length(fname)-4;
%   data = eval(fname(1:N));
%   disp('Data input complete')
if nargin<2|isempty(ns);ns=2048;end;
if nargin<3|isempty(plt);plt=0;end;

% Set variables
     shift = 0;    % Orientation of the array (zero, if West) 
  np=8;                         % 
  %ns=2048;                       % Number of points per block
%  nb=175;                       % Number of blocks
  nb = floor(length(data)/ns);
  nf=ns/2;                      % Number of frequencies in spectrum
  df=80.;                       % Sampling frequency (Hz)


% Set up matrix to force CPSD to be positive definate
  eps=ones(np,np)*0.99999; 
  for l=1:np
    eps(l,l)=1;
  end

  depth=3.0;                      % Water depth
  
% Define the orientation of array sensors of the external array
% It is assumed that Sensor 1 is at 0 degrees and radius R1
  R1=0.5;                      % Outer radius of sensors
  anginc=360/(np-1);
  r=0:anginc:359;               % Angles of array sensors
  r=r*pi/180;                   % Convert to radians
  r(np)=0;

% Define the orientation of array sensors of the internal array
% The central sensor is the same
%   R3(1)=sqrt(0.0032);
%   R3(2)=0.04;                   % Outer radius of sensors
%   r3(1)=-173;                 % Angles of array sensors
%   r3(2)=-128;
%   r3=r3*pi/180;                   % Convert to radians  

% Define the orientation of array sensors of the external array
% It is assumed that Sensor 1 is at 0 degrees and radius R1
  anginc=360/(np-1);
  r=0:anginc:359;               % Angles of array sensors
  r=r*pi/180;                   % Convert to radians
  r(np)=0;


% Set directional resolution for directional spectrum
  res=2;                        % Resolution in degrees;
  nang=360/res;                 % Number of angles in directional spectrum
  ang=0:res:359;                % Angles in degrees
  ang=ang*pi/180;               % Angles in radians

% Compute the window 
  wt=hanning(ns);
  factor=norm(wt)^2;            % Normalization factor
  wt=[wt wt wt wt wt wt wt wt];    % Set up weight vector for each sensor

% Compute the frequencies for each spectral bin
  f=(0:nf-1);
  f=f*df/ns;
  f=f';

% Compute the wavenumbers
%   temp1=cg(2*pi*f,depth);
temp1=wavek(f,depth);
  k=temp1(:,1);

% Define sqrt(-1)
  ii=sqrt(-1);

% Define complex phase lags
  x1=zeros(nf,nang);            % Zero the complex phase lag arrays
  x2=zeros(nf,nang);
  x3=zeros(nf,nang);
  x4=zeros(nf,nang);
  x5=zeros(nf,nang);
  x6=zeros(nf,nang);
  x7=zeros(nf,nang);
  x8=zeros(nf,nang);
  
  for m=1:nang
    x1(:,m)=exp(-ii*k(:)*R1*cos(r(1)-ang(m)));          % Sensor 1
    x2(:,m)=exp(-ii*k(:)*R1*cos(r(2)-ang(m)));          % Sensor 2
    x3(:,m)=exp(-ii*k(:)*R1*cos(r(3)-ang(m)));          % Sensor 3
    x4(:,m)=exp(-ii*k(:)*R1*cos(r(4)-ang(m)));          % Sensor 4
    x5(:,m)=exp(-ii*k(:)*R1*cos(r(5)-ang(m)));          % Sensor 5
    x6(:,m)=exp(-ii*k(:)*R1*cos(r(6)-ang(m)));          % Sensor 6
    x7(:,m)=exp(-ii*k(:)*R1*cos(r(7)-ang(m)));          % Sensor 7
    x8(:,m)=exp(-ii*k(:)*0);                                % Sensor 8 - Centre sensor
%     x6(:,m)=exp(-ii*k(:)*0);                            % Sensor 6 - Centre sensor
%     x7(:,m)=exp(-ii*k(:)*R3(1)*cos(r3(1)-ang(m)));      % Sensor 7
%     x8(:,m)=exp(-ii*k(:)*R3(2)*cos(r3(2)-ang(m)));      % Sensor 8 
  end
  

% Initialize the data array to increase speed
%  data=zeros(ns*nb,np+ni);

% Input the names and wind speed and direction
%  load /data/IRY/LGEORGE/DIRSPC/DATA/wind.dir
%  for i=1:length(wind)
%    names(i,1:11)=[int2str(wind(i,1)), '.a6'];          % Convert names to strings
%  end
 

% Loop through the file names
%  for k=1:length(names)
%    disp(['Loading ',names(k,1:11)])
%    eval(['!cp /data/IRY/LGEORGE/DIRSPC/DATA/',names(k,1:11),' data.dat'])
%    load data.dat;
%    disp('Data input complete')
%
%   Extract wind speed and direction
%    wspeedn=wind(k,2)
%    wdirn=wind(k,3)
%    wspeed=num2str(wspeedn);
%    wdir=num2str(wdirn);

%   Calculate the wave direction equivalent
%    wavdir=wdirn+180;
%    if wavdir > 360
%      wavdir=wavdir-360;
%    end

%   Rearrange the pole numbering
    ws(:,1)=data(:,1);
    ws(:,2)=data(:,2);
    ws(:,3)=data(:,3);
    ws(:,4)=data(:,4);
    ws(:,5)=data(:,5);
    ws(:,6)=data(:,6);
    ws(:,7)=data(:,7);                                 % Rearrange data
    ws(:,8)=data(:,8);                                 % Rearrange data

%     wsi(:,1)=data(:,1);
%     wsi(:,2)=data(:,2);
%     wsi(:,3)=data(:,3);
%     wsi(:,4)=data(:,4);
%     wsi(:,5)=data(:,5);
%     wsi(:,6)=data(:,6);
%     wsi(:,7)=data(:,7);                                 % Rearrange data
%     wsi(:,8)=data(:,8);                                 % Rearrange data

%   Detrend the data
    ws=detrend(ws,0);

%   Apply calibration factor
%     ws(:,1)=-0.00023*ws(:,1);
%     ws(:,2)=-0.00023*ws(:,2);
%     ws(:,3)=-0.00023*ws(:,3);
%     ws(:,4)=-0.00023*ws(:,4);
%     ws(:,5)=-0.00023*ws(:,5);
%     ws(:,6)=-0.00023*ws(:,6);
%     ws(:,7)=-0.00023*ws(:,7);
%     ws(:,8)=-0.00023*ws(:,8);

	Hs = 4*mean(std(ws));	% significant wave height	
%   Plot time series to check for spikes etc.

%   Adjust records so that they have the same variance
    sws=std(ws);
    asws=sqrt(mean(sws.^2));
    for m=1:np
      ws(:,m)=ws(:,m)*asws/sws(m);
    end

%   Initialize arrays
    ct=zeros(nf,(np)*(np));                         % Temporary store for cross spectra
    c=zeros(nf,(np)*(np));                          % The cross spectra
    S=zeros(nf,1);                              % 1-D spectrum
    s=zeros(nf,nang);                           % 2-D spectrum

%   Determine the cross spectra
    for i=1:nb
      is=1+(i-1)*ns;
      wss=ws(is:is+ns-1,:);                     % Extract ns points
      wss=wt.*wss;                              % Include window function

%     Calculate the fft for all np+ns sensors
      ft=fft(wss);
      ft1=ft(2:nf+1,:);                 % Exclude f=0 Hz
      ft1(nf,:)=ft1(nf,:)/2;            % The last frequency includes folded energy

%     Calculate and store cross spectra
      ij=0;
      for l=1:np
        for m=1:np
          ij=ij+1;
          ct(:,ij)=ft1(:,l).*conj(ft1(:,m))./(abs(ft1(:,l)).*abs(ft1(:,m)));
        end
        S(:)=S(:)+ft1(:,l).*conj(ft1(:,l));     % The 1-D spectrum
      end

%     Collect the block averages
      c=c+ct;
    end

%   Include number of blocks in average
    c=c/nb;

%   Normalize the 1-D spectrum and accounts for number of blocks and sensors
    S=S/(nb*(np));                        % Accounts for number of blocks and sensors
    S=S/factor;                         % Correction for window
    S=2*S/df;                           % Account for neg. f and sampling rate

%   Reconstruct matrices and invert
%   Take each frequency and each angle in turn
    for w=1:nf
      csub=c(w,:);                      % Extract all combinations at this freq
      cm=reshape(csub,(np),(np));           % Make into square matrix
      cm=cm.';                          % Transpose the matrix
      cm=cm.*eps;                       % Ensure diag. terms dominate
      cinv=inv(cm);                        % Invert the CSPD matrix
      for th=1:nang
        xm=[x1(w,th);x2(w,th);x3(w,th);x4(w,th);x5(w,th);x6(w,th);x7(w,th);x8(w,th)];
        s(w,th)=1/(xm'*cinv*xm);           % The spreading function
      end
    end

%   Ensure s and S are real
    s=real(s);
    S=real(S);

%   Normalize the directional spectrum to agree with the 1-D spectrum
    for i=1:nf
      area=sum(s(i,:))*res*pi/180; 
      s(i,:)=s(i,:)*S(i)/area;
    end

%   Rotate values to allow for array geometry and convert NS directions
    ang1=ang*180/pi;
   ang2=210-ang1;
%     ang2=90-ang1;
%    ang2=180-ang1;
    for th=1:nang
%      if ang2(th) >= 360
%        ang2(th)=ang2(th)-360;
      if ang2(th) < 0
        ang2(th)=ang2(th)+360;
      end
    end

    for th=1:nang
      idxold=round(ang1(th)/res)+1;
      idxnew=round(ang2(th)/res)+1;
      s1(:,idxnew)=s(:,idxold);
    end

%
%   Produce output plots
%
%   Plot one dimensional spectrum
%    S = S*4*pi;
%     subplot(2,2,1) 
%     loglog(f(2:nf-1),S(2:nf-1),'g')
%     xlabel('f [Hz]')
%     ylabel('E')
%
 %       Shifting the angle
          ang1 = ang1-shift;
          le = length(ang1);
          if ang1(1) < 0,
                  i = 1;
                  while ang1(i) < 0,
                          ang1(i) = ang1(i)+360;
                          i = i+1;
                  end
                  M = i-1;
                  for i = 1:le-M,
                          si(:,i) = s1(:,i+M);
                  end
                  for i = 1:M,
                          si(:,le-M+i) = s1(:,i);
                  end
          else
                  j = 1;
                  for i = 1:le,
                          if ang1(i) >= 360,
                                  ang1(i) = ang1(i)-360;
                                  j = j+1;
                          end
                  end
                  M = j-1;
          end
          if M ~= 0,
                  	for i = 1:M,
                       	   si(:,i) = s1(:,le-M+i);
                  	end
                  	for i = 1:le-M,
                       	   si(:,le-M+i) = s1(:,i);
                  	end
		  else
			si = s1;
		  end

          s1 = si;
          ang1 = sort(ang1);
          kj=find(ang1>180);
          ang1(kj)=ang1(kj)-360;
          [o,p]=max(max(s1));
          ang1(p);
          ang1=ang1-ang1(p);
          AA=find(ang1>180);
          BB=find(ang1<-180);
          B3=isempty(BB);A3=isempty(AA);
          if A3==0;ang1(find(ang1>180))=ang1(find(ang1>180))-360;end
          if B3==0;ang1(find(ang1<-180))=ang1(find(ang1<-180))+360;end;
          [ang1,P]=sort(ang1);
          nrj=s1(:,P);
          
          E=surf_framespec(nrj,f,ang1,[],depth);
          E.note='Freq in Hz!!!!';

%   Contour plot of spectrum
if plt > 0
    contour(ang1,f(1:50),nrj(1:50,:),25)
    %contour(ang1,f(11:30),s1(11:30,:),10)
    xlabel('\theta')
    ylabel('f [Hz]')
    grid
end

