function s = disp2stress(u,Q)
    % s = disp2stress(u,Q)
    %
    % Calculate traction stresses at z=h from displacements measured at z=z0*h.
    % For each component, FFT(sigma)_i = (Q_ij)*FFT(u)_j. See Supporting Information
    % of Xu et al. PNAS 2010 for more details (the theory).
    %
    % CREATED
    %   Ye Xu 05/2009
    % MODIFIED
    %   Ye Xu and Eric Dufresne 10/20/2009
    %   Ye Xu and Eric Dufresne 10/22/2009
    %     - to handle normal force = 0 assumption
    %   CWM 2013-05-09
    %     - Reformatted code
    %     - Revised and added comments
    %     - Removed unnecessary shft=0; circshift(*,[shft,shft]); in calc of FFT(s).
    %     - Replaced fftshift(circshift(*,[1,1])) with ifftshift(*) in calc of IFFT(s).
    %       The former breaks when N is odd, whereas the latter always
    %       does the right thing (see doc ifftshift)
    %
    %INPUTS
    %   u: a structure array such that u(i).u is a scalar field giving the
    %      displacement in the ith direction, i=3 is orthogonal
    %      to the substrate surface. u should be interpolated to a square grid.
    %      If u is 2D, the stress in the z-direction is assumed to be zero.
    %
    %   Q: the output of calcQ. 
    %
    %OUTPUTS
    %   s: a structure array such that s(i).s is a scalar field giving the
    %      traction stress in the ith direction, where i=3 is orthogonal
    %      to the substrate surface
    
    
    % Get the number of spatial dimensions (see defintion of u above)
    D = length(u);
    
    % Take the FFT of each component of the displacement field
    for i=1:D
        u(i).U = fftshift(fft2(u(i).u));
    end
    
    % Calculate the FFT of each component of the stress field
    for i=1:D
        s(i).S = 0;
        for j=1:D
            s(i).S = s(i).S + (Q(i,j).Q).*(u(j).U);
        end
    end
    
    % IFFT the stresses
    for i=1:D
        s(i).s =ifft2(ifftshift(s(i).S),'symmetric');
    end
    
end