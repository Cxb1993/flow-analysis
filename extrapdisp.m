function out = extrapdisp(u)
% out  = extrapdisp(u)
%  extrapdisp crudely extrapolates the values of a scalar field, u, to 
%  replace NaN values.  it is designed to work when all the NaNs are
%  connected to the edge through other NaNs
%  INPUT:
%    u: a 2D matrix if values with NaNs 'framing' the data
%  OUTPUT:
%    out: a 2D matrix with the same dimension as u with values from the
%    edge of the data set extrapolated toward the edge of the matrix
%  Created by Eric R. Dufresne based on earlier code by Ye Xu.  10/20/2009

[nr,nc]=size(u);

for i=1:nr
   ind1 = find(~isnan(u(:,i)),1,'first');
   ind2 = find(~isnan(u(:,i)),1,'last');
   
   if ~isempty(ind1) & ind1~=1
       %extrapolate to the left
       u(1:ind1-1,i)=u(ind1,i);
   end
   if ~isempty(ind2) & ind2~=nc;
       %extrapolate to the right
       u((ind2+1):end,i)=u(ind2,i);
   end
end

for i=1:nc
   ind1 = find(~isnan(u(i,:)),1,'first');
   ind2 = find(~isnan(u(i,:)),1,'last');
   
   if ~isempty(ind1) & ind1~=1
       %extrapolate to the top
       u(i,1:ind1-1)=u(i,ind1);
   end
   if ~isempty(ind2) & ind2~=nc;
       %extrapolate to the bottom
       u(i,(ind2+1):end)=u(i,ind2);
   end
end

out=u;
end