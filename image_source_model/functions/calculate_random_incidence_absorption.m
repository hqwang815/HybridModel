%% function name
% Function to calculate random incidence absorption according to Kuttruff Room Acoustics 
%% Syntax
% # [output argument(s)] = functionname(input arguments)
%% Description
% Description
% Input arguments
% Output arguments
%% Example(s)
% Related functions
function [auni] = calculate_random_incidence_absorption(specificImpedance)


for iSurface=1:1:size(specificImpedance,1)


    
    
    
    
Z=specificImpedance(iSurface,:);

check = sum(imag(Z)~=0);


if check ~=0

 u=atan(imag(Z)./real(Z));
 
auni(iSurface,:)=8./(abs(Z).^2).*cos(u) .* (abs(Z)+(cos(2.*u)./(sin(u))).*atan((abs(Z)...
    .*sin(u))./(1+abs(Z).*cos(u)))-cos(u).*log(1+2.*abs(Z).*cos(u)+abs(Z).^2));


%% Discrete integration (Can be used for real-valued impedance)

elseif check ==0
    
s=50;
phi=(0:pi/s:.5*pi);
nSteps=length(phi);

for i=1:1:nSteps
    
    R= ((Z*cos(phi(i)))-1)./((Z*cos(phi(i)))+1);
    a= (1-abs(R).^2);  
    aTotal(i,:)=a*sin(2*phi(i));
    
    
end

auni(iSurface,:)=sum(aTotal)*(pi/s);

end


end