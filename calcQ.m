function Q = calcQ(z,h,EM,nu,nr,dx,D)
    % function Q = calcQ(xy,z,h,EM,nu)
    % 
    % Calculate the matrix Q = P(h)/M(z) in Fourier space. See Supporting Information
    % of Xu et al. PNAS 2010 for more details (the theory).
    % 
    % CREATED
    %   YX 05/2009
    %     - based on calculation of Larry Wilen
    % MODIFIED
    %   Ye Xu, Eric Dufresne 10/20/2009
    %     - to include 2D stuff
    %   Eric Dufresne 10/11/2010
    %     - to pull out the exponentials from sinh and cosh to stabilize numerics
    %       for thick films and stuff very near a surface
    %   CWM 2013-05-10
    %     - Reformatted code
    %     - Revised and added comments
    %     - Removed unnecessary variables nrow and ncol
    %     - Changed determinant variable names from dt* to det* throughout
    %     - Renamed det to detQ2inv in calculation of Q2
    %     - Fixed calculation of ncenter, a, and b for even N.
    %   CWM 2013-05-20
    %     - Added Rob's "simpler" method for calculating Q2
    %
    % INPUTS
    %   z: distance of the observation plane above the zero displacment surface
    %   h: thickness of the film
    %   EM: Youngs Modulus
    %   nu: Poisson's ratio
    %   nr: number of rows and columns across field (must be square)
    %   dx: distance between succesive rows and columns
    %   D(optional): number of dimensions you have displacements for.  If D=2,
    %                the stress in the z-direction is assumed to be zero.
    %
    % OUTPUTS
    %   Q: a structure array containing fields .Q and .Qinv.
    %     - Q(i,j).Q: a scalar field showing the values of Qij in Fourier space.
    %     - Q(i,j).Qinv: a scalar field showing the values of the inverse of Qij in Fourier space.
    
    
    % Create grid in Fourier space
    fov = dx*nr;
    L = fov/2;
    N = nr;
    n = floor(N/2);
    ncenter = n+1; 
    
    if mod(N,2)==1
        a = repmat(linspace(-n*pi/L,n*pi/L,N),[N,1]);
        b = repmat((linspace(-n*pi/L,n*pi/L,N))',[1,N]); 
    else
        % FFT puts the zero-frequency mode at N/2+1 for even N,
        % so you lose the largest positive wavenumber.
        a = repmat(linspace(-n*pi/L,(n-1)*pi/L,N),[N,1]);
        b = repmat((linspace(-n*pi/L,(n-1)*pi/L,N))',[1,N]); 
    end
    k = sqrt(a.*a+b.*b);
    
    % Construct the matrix M(z) and M(h) (in Fourier space)
    %  NOTE: We have factored exp(kh)/2 and exp(kz)/2 out of sinh and cosh.
    %        We account for this below when constructing Q.
    Mz(1,1).M = (((3-4*nu)*a.^2+4*(1-nu)*b.^2).*(1 - exp(-2*k*z))./(4*k.^3*(1-nu)))+((a.^2*z).*(1 + exp(-2*k*z))./(4*k.^2*(1-nu)));
    Mh(1,1).M = (((3-4*nu)*a.^2+4*(1-nu)*b.^2).*(1 - exp(-2*k*h))./(4*k.^3*(1-nu)))+((a.^2*h).*(1 + exp(-2*k*h))./(4*k.^2*(1-nu)));
    Mz(1,2).M = (-b.*a./(4*k.^3*(1-nu))).*(1 - exp(-2*k*z))+(b.*a*z./(4*k.^2*(1-nu))).*(1 + exp(-2*k*z));
    Mh(1,2).M = (-b.*a./(4*k.^3*(1-nu))).*(1 - exp(-2*k*h))+(b.*a*h./(4*k.^2*(1-nu))).*(1 + exp(-2*k*h));
    Mz(1,3).M = complex(0,-z*a.*(1 - exp(-2*k*z))./(2*k*(1-2*nu)));
    Mh(1,3).M = complex(0,-h*a.*(1 - exp(-2*k*h))./(2*k*(1-2*nu)));
    
    Mz(2,1).M = -b.*a.*(1 - exp(-2*k*z))./(4*k.^3*(1-nu))+b.*a.*z.*(1 + exp(-2*k*z))./(4*k.^2*(1-nu));
    Mh(2,1).M = -b.*a.*(1 - exp(-2*k*h))./(4*k.^3*(1-nu))+b.*a.*h.*(1 + exp(-2*k*h))./(4*k.^2*(1-nu));
    Mz(2,2).M = ((3-4*nu)*b.^2+4*(1-nu)*a.^2).*(1 - exp(-2*k*z))./(4*k.^3*(1-nu))+b.^2*z.*(1 + exp(-2*k*z))./(4*k.^2*(1-nu));
    Mh(2,2).M = ((3-4*nu)*b.^2+4*(1-nu)*a.^2).*(1 - exp(-2*k*h))./(4*k.^3*(1-nu))+b.^2*h.*(1 + exp(-2*k*h))./(4*k.^2*(1-nu));
    Mz(2,3).M = complex(0,-z*b.*(1 - exp(-2*k*z))./(2*k*(1-2*nu)));
    Mh(2,3).M = complex(0,-h*b.*(1 - exp(-2*k*h))./(2*k*(1-2*nu)));
    
    Mz(3,1).M = complex(0,-z*a.*(1 - exp(-2*k*z))./(4*k*(1-nu)));
    Mh(3,1).M = complex(0,-h*a.*(1 - exp(-2*k*h))./(4*k*(1-nu)));
    Mz(3,2).M = complex(0,-z*b.*(1 - exp(-2*k*z))./(4*k*(1-nu)));
    Mh(3,2).M = complex(0,-h*b.*(1 - exp(-2*k*h))./(4*k*(1-nu)));
    Mz(3,3).M = (3-4*nu)*(1 - exp(-2*k*z))./(2*k*(1-2*nu))-z*(1 + exp(-2*k*z))/(2*(1-2*nu));
    Mh(3,3).M = (3-4*nu)*(1 - exp(-2*k*h))./(2*k*(1-2*nu))-h*(1 + exp(-2*k*h))/(2*(1-2*nu));
    
    % Construct the matrix N(h) = dM/dz(h) (in Fourier space)
    Nh(1,1).N = a.^2*h.*(1 - exp(-2 *k * h) )./(4*k*(1-nu))+(1 + exp(-2*k*h));
    Nh(1,2).N = b.*a*h.*(1 - exp(-2 *k * h) )./(4*k*(1-nu));
    Nh(1,3).N = complex(0,-a.*(1 - exp(-2 *k * h) )./(2*k*(1-2*nu))-h*a.*(1 + exp(-2*k*h))./(2*(1-2*nu)));
    
    Nh(2,1).N = b.*a*h.*(1 - exp(-2 *k * h) )./(4*k*(1-nu));
    Nh(2,2).N = b.^2*h.*(1 - exp(-2 *k * h) )./(4*k*(1-nu))+(1 + exp(-2*k*h));
    Nh(2,3).N = complex(0,-b.*(1 - exp(-2 *k * h) )./(2*k*(1-2*nu))-h*b.*(1 + exp(-2*k*h))./(2*(1-2*nu)));
    
    Nh(3,1).N = complex(0,-a.*(1 - exp(-2 *k * h) )./(4*k*(1-nu))-h*a.*(1 + exp(-2*k*h))./(4*(1-nu)));
    Nh(3,2).N = complex(0,-b.*(1 - exp(-2 *k * h) )./(4*k*(1-nu))-h*b.*(1 + exp(-2*k*h))./(4*(1-nu)));
    Nh(3,3).N = -k*h.*(1 - exp(-2 *k * h) )./(2*(1-2*nu))+(1 + exp(-2*k*h));
    
    % Construct the matrix P(h) = A1.*M(h) + A2.*N(h) (in Fourier space)
    % where A1 = [0 0 i*a; 0 0 i*b; i*nu*a i*nu*b 0]
    %       A2 = [1 0 0 ; 0 1 0 ; 0 0 1-nu]
    Ph(1,1).P = complex(0,a).*Mh(3,1).M + Nh(1,1).N;
    Ph(1,2).P = complex(0,a).*Mh(3,2).M + Nh(1,2).N;
    Ph(1,3).P = complex(0,a).*Mh(3,3).M + Nh(1,3).N;
    
    Ph(2,1).P = complex(0,b).*Mh(3,1).M + Nh(2,1).N;
    Ph(2,2).P = complex(0,b).*Mh(3,2).M + Nh(2,2).N;
    Ph(2,3).P = complex(0,b).*Mh(3,3).M + Nh(2,3).N;
    
    Ph(3,1).P = complex(0,nu*a).*Mh(1,1).M + complex(0,nu*b).*Mh(2,1).M + (1-nu)*Nh(3,1).N;
    Ph(3,2).P = complex(0,nu*a).*Mh(1,2).M + complex(0,nu*b).*Mh(2,2).M + (1-nu)*Nh(3,2).N;
    Ph(3,3).P = complex(0,nu*a).*Mh(1,3).M + complex(0,nu*b).*Mh(2,3).M + (1-nu)*Nh(3,3).N;
    
    % Construct the inverse of M
    %   Determinant of M
    detMz = Mz(1,1).M.* (Mz(2,2).M.*Mz(3,3).M - Mz(3,2).M.*Mz(2,3).M) ...
          - Mz(1,2).M.* (Mz(2,1).M.*Mz(3,3).M - Mz(3,1).M.*Mz(2,3).M) ...
          + Mz(1,3).M.* (Mz(2,1).M.*Mz(3,2).M - Mz(3,1).M.*Mz(2,2).M);
    %   Inverse of the determinant
    invdetMz = 1./detMz;
    %   Inverse of M
    Mzinv(1,1).M = invdetMz .* (Mz(2,2).M.*Mz(3,3).M - Mz(2,3).M.*Mz(3,2).M);
    Mzinv(1,2).M = invdetMz .* (Mz(1,3).M.*Mz(3,2).M - Mz(1,2).M.*Mz(3,3).M);
    Mzinv(1,3).M = invdetMz .* (Mz(1,2).M.*Mz(2,3).M - Mz(1,3).M.*Mz(2,2).M);
    Mzinv(2,1).M = invdetMz .* (Mz(2,3).M.*Mz(3,1).M - Mz(2,1).M.*Mz(3,3).M);
    Mzinv(2,2).M = invdetMz .* (Mz(1,1).M.*Mz(3,3).M - Mz(1,3).M.*Mz(3,1).M);
    Mzinv(2,3).M = invdetMz .* (Mz(1,3).M.*Mz(2,1).M - Mz(1,1).M.*Mz(2,3).M);
    Mzinv(3,1).M = invdetMz .* (Mz(2,1).M.*Mz(3,2).M - Mz(2,2).M.*Mz(3,1).M);
    Mzinv(3,2).M = invdetMz .* (Mz(1,2).M.*Mz(3,1).M - Mz(1,1).M.*Mz(3,2).M);
    Mzinv(3,3).M = invdetMz .* (Mz(1,1).M.*Mz(2,2).M - Mz(1,2).M.*Mz(2,1).M);
    
    % Construct matrix Q = <stuff> * Ph * Mzinv
    for i=1:3
        for j=1:3
            Q(i,j).Q = Ph(i,1).P.*Mzinv(1,j).M + Ph(i,2).P.*Mzinv(2,j).M + Ph(i,3).P.*Mzinv(3,j).M;
            
            % Account for the exp(kh)/2 and exp(kz)/2 we factored out of M above.
            % exp(kh)/2 comes from P(h), exp(-kz)/2 comes from Minv(z) -- the 2s cancel.
            Q(i,j).Q = exp( k * (h-z) ).* Q(i,j).Q;
        end
    end
    
    %   Replace the elements of Q for k = 0;
    Q(1,1).Q(ncenter,ncenter) = 1/z;
    Q(1,2).Q(ncenter,ncenter) = 0;
    Q(1,3).Q(ncenter,ncenter) = 0;
    
    Q(2,1).Q(ncenter,ncenter) = 0;
    Q(2,2).Q(ncenter,ncenter) = 1/z;
    Q(2,3).Q(ncenter,ncenter) = 0;
    
    Q(3,1).Q(ncenter,ncenter) = 0;
    Q(3,2).Q(ncenter,ncenter) = 0;
    Q(3,3).Q(ncenter,ncenter) = (1-nu)/z;
    
    %   Multiply by constant pre-factors -- the <stuff> from above.
    for i=1:3
        Q(1,i).Q = Q(1,i).Q * EM/(2*(1+nu));
        Q(2,i).Q = Q(2,i).Q * EM/(2*(1+nu));
        Q(3,i).Q = Q(3,i).Q * EM/((1+nu)*(1-2*nu));
    end
    
    
    % Construct the inverse of Q: Qgreen
    
    % Construct the inverse of P(h)
    %   Determinant of P(h)
    detPh = Ph(1,1).P.* (Ph(2,2).P.*Ph(3,3).P - Ph(3,2).P.*Ph(2,3).P) ...
          - Ph(1,2).P.* (Ph(2,1).P.*Ph(3,3).P - Ph(3,1).P.*Ph(2,3).P) ...
          + Ph(1,3).P.* (Ph(2,1).P.*Ph(3,2).P - Ph(3,1).P.*Ph(2,2).P);
    %   Inverse of determinant
    invdetPh = 1./detPh;
    % invdetMz=conj(detMz)./( detMz.*conj(detMz) + 1e-6*mean(mean(detMz.*conj(detMz))));
    %   Inverse of P(h)
    Phinv(1,1).P = invdetPh .* (Ph(2,2).P.*Ph(3,3).P - Ph(2,3).P.*Ph(3,2).P);
    Phinv(1,2).P = invdetPh .* (Ph(1,3).P.*Ph(3,2).P - Ph(1,2).P.*Ph(3,3).P);
    Phinv(1,3).P = invdetPh .* (Ph(1,2).P.*Ph(2,3).P - Ph(1,3).P.*Ph(2,2).P);
    
    Phinv(2,1).P = invdetPh .* (Ph(2,3).P.*Ph(3,1).P - Ph(2,1).P.*Ph(3,3).P);
    Phinv(2,2).P = invdetPh .* (Ph(1,1).P.*Ph(3,3).P - Ph(1,3).P.*Ph(3,1).P);
    Phinv(2,3).P = invdetPh .* (Ph(1,3).P.*Ph(2,1).P - Ph(1,1).P.*Ph(2,3).P);
    
    Phinv(3,1).P = invdetPh .* (Ph(2,1).P.*Ph(3,2).P - Ph(2,2).P.*Ph(3,1).P);
    Phinv(3,2).P = invdetPh .* (Ph(1,2).P.*Ph(3,1).P - Ph(1,1).P.*Ph(3,2).P);
    Phinv(3,3).P = invdetPh .* (Ph(1,1).P.*Ph(2,2).P - Ph(1,2).P.*Ph(2,1).P);
    
    %   Construct the inverse of Q, Qinv = (1/stuff)*M(z)*Pinv(h)
    for i=1:3
        for j=1:3
            Q(i,j).Qinv = Mz(i,1).M.*Phinv(1,j).P + Mz(i,2).M.*Phinv(2,j).P + Mz(i,3).M.*Phinv(3,j).P;
            Q(i,j).Qinv = exp(k * (z - h) ) .* Q(i,j).Qinv;
        end
    end
    Q(1,1).Qinv(ncenter,ncenter) = z;
    Q(1,2).Qinv(ncenter,ncenter) = 0;
    Q(1,3).Qinv(ncenter,ncenter) = 0;
    
    Q(2,1).Qinv(ncenter,ncenter) = 0;
    Q(2,2).Qinv(ncenter,ncenter) = z;
    Q(2,3).Qinv(ncenter,ncenter) = 0;
    
    Q(3,1).Qinv(ncenter,ncenter) = 0;
    Q(3,2).Qinv(ncenter,ncenter) = 0;
    Q(3,3).Qinv(ncenter,ncenter) = z/(1-nu);
    
    %   This is the (1/stuff) part.
    for i=1:3
        Q(i,1).Qinv = Q(i,1).Qinv / (EM/(2*(1+nu)));
        Q(i,2).Qinv = Q(i,2).Qinv / (EM/(2*(1+nu)));
        Q(i,3).Qinv = Q(i,3).Qinv / (EM/((1+nu)*(1-2*nu)));
    end
    
    
    % Special case: D=2
    if nargin==7 & D==2
        
        % Two different ways of calculating Q2
        %   They are algebraically equivalent. We retain both, but use the latter.
        
        % Method 1 (Eric's original method)
        % % Keep only the first 2x2 of Qinv
        % for i=1:2
        %     for j=1:2
        %         Q2(i,j).Qinv = Q(i,j).Qinv;
        %     end
        % end
        % % Q2.Q is the inverse of Q2.Qinv
        % detQ2inv = Q2(1,1).Qinv.*Q2(2,2).Qinv - Q2(1,2).Qinv.*Q2(2,1).Qinv;
        % Q2(1,1).Q =  Q2(2,2).Qinv./detQ2inv;
        % Q2(1,2).Q = -Q2(1,2).Qinv./detQ2inv;
        % Q2(2,1).Q = -Q2(2,1).Qinv./detQ2inv;
        % Q2(2,2).Q =  Q2(1,1).Qinv./detQ2inv;
        
        % Method 2 (Rob's method)
        % Calculate Q2 directly based on setting sigma33=0 and eliminating u33
        Q2(1,1).Q = Q(1,1).Q - Q(1,3).Q.*Q(3,1).Q./Q(3,3).Q;
        Q2(1,2).Q = Q(1,2).Q - Q(1,3).Q.*Q(3,2).Q./Q(3,3).Q;
        Q2(2,1).Q = Q(2,1).Q - Q(2,3).Q.*Q(3,1).Q./Q(3,3).Q;
        Q2(2,2).Q = Q(2,2).Q - Q(2,3).Q.*Q(3,2).Q./Q(3,3).Q;
        % Calculate Q2inv directly
        detQ2 = Q2(1,1).Q.*Q2(2,2).Q - Q2(1,2).Q.*Q2(2,1).Q;
        Q2(1,1).Qinv =  Q2(2,2).Q./detQ2;
        Q2(1,2).Qinv = -Q2(1,2).Q./detQ2;
        Q2(2,1).Qinv = -Q2(2,1).Q./detQ2;
        Q2(2,2).Qinv =  Q2(1,1).Q./detQ2;
        
        Q3 = Q; % Preserve 3D Q...
        Q = Q2;
        
    end
    
end