%-----------------------------------------------
function [lon,lat,Bval]=read_Bmaps_CS(fileIn)

%
ff=load(fileIn);
lon=ff(:,1);
lat=ff(:,2);
Bval=ff(:,3);

end
%-----------------------------------------------
