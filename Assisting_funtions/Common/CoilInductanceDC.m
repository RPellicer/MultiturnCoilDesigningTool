function [L,Lo,M]= CoilInductanceDC(s,rc,zc)
%% function [L,M]=L_M_Coaxial(R,nv,nl,do)
    % This function computes the self and mutual inductance of a coaxial coil
    % configuration having
    %R -> for internal radius [m]
    %nv-> vertical turns
    %nl-> number of layers
    %s-> distance in between wires [m]
    %dr -> radial filling factor [m]
    %dz -> axial filling factor  [m]
    % The calculation follows the formulation discussed in Jose et al
    % "Evaluation of the Inductance, the Resitance and the Capacitance of
    % Coaxial Inductos at Low Frequencies-A Review
    % 
    %---- Lo ---- Self inductance according to Eq.(2)
    % Lo= mu*L1n*Fn1;
    % L1n -> vector row having 
    % L1n =[2r1-a,2r2-a,...,2rn-a]
    % and
    % Fn1-> colum vector having
    % Fn1 = (Kn1-En1)
    % Kn1= Trampose[1-1/2*k1,1-1/2*k2,...1-1/2*kn]
    % En1= Trampose[E(k1),E(k2),...E(kn)]
    % and
    % ki =
    % 4[(r1*(r1-a),r2*(r2-a),...rn*(rn-a)]./[(2*r1-a)^2,(2*r2-a)^2,...(2*rn-a)^2]
    % with
    % ri = R +(2*r1-a)*a-(i-1)dr
    %---- M ---- Mutual inductance according to Eq.(10)
    % M = 2*mu *Transpose[(Mnn.*Fnn)*Transpose(1)]*Transpose(1)
    % Mnn = [0,sqrt(r2*r1)/k2,1,sqrt(r3*r1)/k3,1....; sqrt(r1*r2)/k1,2, 0, sqrt(r3*r2)/k3,2]
    %---- L =Lo + M---- Total inductance Eq.(14)
    % Call internal function XY_wire(nv,nl,Rm,do,dr,dz)
    % returning vectors for both radial, r and axial positions, zi of the wires
    % Constant "mu", the mangetic permeability of the free space
    % mu = mu=4*pi*10^(-10);      [H/mm]
    % Ruben Pellicer @ CAI 2015
    % v2. Feb 15
    n = length(rc);                                      % total number of turns 
    mu=4*pi*1e-7;                               % magnetic permeability in H/m
    r=rc(:)';z=zc(:)';                              % r(1xn),z(1xn) are row vectors having wires' center of coordinates
    a=s/2;                                               
    % ---- Self Inductance -------
    L1n = (2*r-a);                                  % row vector with dimensions (1xn)
    k = 4*(r.*(r-a))./(2*r-a).^2;  % k(1xn)       -> row vector having the coefficient of the eliptic integral   
%     ellipkeTol = eps;          
%     ellipkeTol = eps*n*1e10; % if the number of loops is to high reduce the tolerance of the elliptic integrals
%     if n>5000
    ellipkeTol = eps*5000*1e10;
%     end
    if (min(min(k))<0) || (max(max(k))>1)
        disp('error, inner diameter of the magnetometer is negative')
    end   
    [Kn1,En1]=ellipke(k,ellipkeTol);                             % E(1xn),K(1xn)-> row vectors for the first and second kind elliptical functions
    Kn1_ex = ((1-(1/2)*(k.^2)).*Kn1) ;                          % Kn1(nx1)-> colum vector; first term of Eq.(3)
    Fn1 = Kn1_ex-En1;                                  % Fn1(nx1) ->colum vector; Eq.(3)  
    Lo= mu*L1n*Fn1';                            % Self inductance in uH
    % ---- Mutual inductance ----
    uno =ones(1,n);
    rmr =(r'*uno+uno'*r).^2;                               % rmr(nxn) -> containing the sumation of all radial positions
    zmz =(z'*uno-uno'*z).^2;                               % zmz (nxn)-> containing the differences on axial positions 
    rr = 4*(r'*r);                                 % rr(nxn)  -> numerator of the k-coefficient 
    knt = rr./(rmr+zmz);                              % K(nxn)   -> matrix for kij 

    if (min(min(knt))<0) || (max(max(knt))>1)
        disp('error, wire insulation is negative')
    end   
    [Knt,Ent]=ellipke(knt,ellipkeTol);   
    Knt_ex = ((1-1/2*(knt.^2)).*Knt); 
    Fnt = Knt_ex-Ent;  
    E = eye(n);                                   % generates a matrix zero everywhere except on diagonals elements
    Mnn = sqrt(r'*r)./knt;  % Matrix of Eq.11. matrix Mnn but diagonal zeros in the diagonal
    for i = 1:n
        Mnn(i,i)=0;
        Fnt(i,i)=0;
    end
    % Calculation of Fnn
    M = 2*mu*(((Mnn.*Fnt)*uno'))'*uno';      % M1x1 of Equation (10)
    L = (M +Lo);   % Total inductance in uH      
end
