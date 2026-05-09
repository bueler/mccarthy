THESE ARE NOTES ON CREATING green_transect.pdf

#!/bin/bash

rm -f foo.nc
ncks -v topg,usrf,lat,lon -d y1,275,275 Greenland_5km_v1.1.nc foo.nc

exit

S = netcdf('foo.nc')
S(1).VarArray.Str

x1 = reshape(S(1).VarArray(6).Data,1,301);
usrf = reshape(S(1).VarArray(5).Data,1,301);
topg = reshape(S(1).VarArray(4).Data,1,301);

in = (x1>-450000) & (x1<350000);

>> max(usrf(in)-topg(in))
ans =  3208.9
>> max(x1(in))-min(x1(in))
ans =  790000
>> 3208/790000
ans =  0.0040608

plot(x1(in),topg(in),x1(in),usrf(in))
hold on
plot(x1(in),0.004 * (usrf(in)-topg(in)),'r',x1(in),zeros(size(x1(in))),'r')
xlabel x, hold off
print -dpdf green_transect.pdf

