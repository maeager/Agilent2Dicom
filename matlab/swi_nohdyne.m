function [swi_n,swi_p,phasemask_p]=swi_nohdyne(mag,pha,order,thr)
%%SWI_NOHDYNE suseptibility weighted imaging, without homedyne filter
%
%
%


narginchk(2,4);
if nargin < 3
    order=4;
end
if nargin <4
    thr=0.05;
end

thd = (thr)*max(mag(:)); 
if (thr ~= 0)
   x = find(mag(:) <= thd);
   pha(x)=pha(x)/10;
end
phasemask_n = (pha+pi)./pi;

phasemask_p = (pi-pha)./pi;

[x] = find(pha<=0);
phasemask_p(x) = 1;



[x] = find(pha>=0);
phasemask_n(x) = 1;

if mod(order,2)==1 || order <1
    phasemask_n=abs(phasemask_n);
end
swi_n = mag.*(phasemask_n.^order);
swi_p = phasemask_p.^order.*mag;