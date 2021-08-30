function [Periods,Heights,Times]=fn_ModifiedMZeroCross(D)
% 
%   Function for to calculate wave heights and wave periods
%   using upcrossing or downcrossing method.
%  
%   This method require will analyze the sign change. The height is define
%   like the diference between the maximum and the minimum 
%   between the portion of record between two successive zero
%   upcrossings or downcrossing.
%   Any height that not comply with this restriction, it
%   is considerate no existed.
%
%  Input: 
%           D = Matrix with two colums.
%               First colum is Time.
%               Second column is the elevation.
%  Output:
%            H = Wave heights
%            T = upcrossing or downcrossing period 
% Example:
%         t=[0:pi/16:6.5*pi]';   % Note that is a column
%         y=sin(t);
%         D=[t,y],
%         [T,H]=MZerocross(D)
%       
% Answer:
%        T =
%            6.2832 6.2832 6.2832
%
%        H =
%
%           2     2     2
%
%
%   Created by R. Hernandez-Walls 
%             Facultad de Ciencias Marinas
%             Universidad Autonoma de Baja California
%             Ensenada, Baja California
%             Mexico.
%             rwalls@uabc.edu.mx
%
%   Modified by Eusebius Marcel
%               Fakultas Teknik Sipil dan Lingkungan
%               Institut Teknologi Bandung
%               Bandung, Jawa Barat
%               Indonesia.
%               eusebius.em@gmail.com
%   
%   The modification is done so the program can recognise different H and T
%   in random waves
%
%   To cite this file, this would be an appropriate format:
%
%   Marcel, Eusebius (2015).
%   ModifiedMZerocross:Zerocrossing Method. A MATLAB file.
%   http://www.mathworks.com/matlabcentral/fileexchange/
%
%   References:
%
%
%
%
% D(:,2)=D(:,2)-mean(D(:,2));
iD=find(diff(sign(D(:,2)))==2);%zero-upcrossing
%iD=find(diff(sign(D(:,2)))==-2);%zero-downcrossing
%put % in front of not wanted choice in line 64 or 65
for k=1:length(iD)-1
Heights(k)=max(D(iD(k):iD(k+1),2))-min(D(iD(k):iD(k+1),2));
Periods(k)=max(D(iD(k):iD(k+1),1))-min(D(iD(k):iD(k+1),1));
Times(k,1) = min(D(iD(k):iD(k+1),1));
Times(k,2) = max(D(iD(k):iD(k+1),1)); 
end
%return