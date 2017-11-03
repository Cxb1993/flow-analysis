function u = stress2disp(s,Q,D)
% u = stress2disp(s,Q,D)
% inverse of disp2stress (which has better documentation)
%calculate force on surface in each pixel from displacements sampled in a
%place below
%F(sigma) = Q.*F(u)
% Set D=2 if traction forces in plane.  Leave blank or set D=3 if traction forces are in plane and out of plane.
%MODIFICATION HISTORY
%   YX 05/2009
%Fourier transfer displacement u to U and shift the low frequence component
%to the center
%QEUSTION
%   - output is the comlex number
% modified jan 2010 by ERD to make compatible with 2D stuff

if nargin<3
    D=3;
end

for i=1:D
 s(i).S=fftshift(fft2(s(i).s));

end

    
    for i=1:D
        u(i).U=0;
         for j=1:D
            u(i).U = u(i).U+(Q(i,j).Qinv) .* s(j).S;     
         end
            u(i).u = ifft2(ifftshift(u(i).U),'symmetric');
           % u(i).u =ifft2(fftshift(circshift(u(i).U,[1,1])),'symmetric');
             
             
    end
    

   
end