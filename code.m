clear all;
clc;

%Input of Basic Parameters
n = 8;						%numbers of lamina
ThetaAll = [0 45 -45 90 90 -45 45 0];		% ply angles in degrees, from top
t   = 0.000125;                      			% lamina thickness
h   = n * t;                         				% thickness of laminate         
distinct = 4;                        			% no of different values of theta 

%Calculation for Reference plane distances
zk(1) = 0;                            			% top layer as reference
for i = 2:n+1
    zk(i) = zk(i-1) + t;
end
zk = zk - h/2;                        % changing to datum plane
for i = 1:n
  zkbar(i) = - (h + t)/2 + i*t;
end

% Ply engineering properties 
E1   = 38.6e9; 					%Pa
E2   = 8.27e9; 					%Pa
G12  = 4.14e9;					%Pa
nu12 = 0.28;
nu21 = nu12 * E2 / E1;

% Ultimate strength of ply
S = [1062e6  -620e6 31e6 -118e6 72e6];

%More Input Parameters
a1 = 8.6e-6;                   				 % coefficients of thermal expansion
a2 = 22.1e-6;
b1 = 0;                           				 % coefficient of moisture expansion
b2 = 0;  
delT = 50;
delC = 0; 

% Thermal, Moisture coefficients in Matrix form
a = [a1 a2 0]';
b = [b1 b2 0]';

% Reduced stiffness matrix(Q) calculation
den = 1 - nu12 * nu21;
Q11 = E1 / den;
Q12 = nu12 * E2 / den;
Q22 = E2 / den;
Q66 = G12;

Q = [Q11 Q12 0; 
    Q12 Q22 0; 
    0 0 Q66];

% Q_bar Matrices (laminate coordinates) and contributions to ABBD matrices
A = zeros(3,3);
B = zeros(3,3);
D = zeros(3,3);

NT = zeros(3,1); 
NH = zeros(3,1);
MT = zeros(3,1); 
MH = zeros(3,1);
N = [1e5 0 0];                    			% N/m   force per unit width of laminate
M = [0 0 0];                  				% N-m/m torque per unit width of laminate

degradation = input('Enter value of degradation');  	%1 = complete, 2 = partial.

CD = 370;         % just to enter into the loop as we dont use 370 degrees
                          % at any time
 
%With degradation
for k =1:distinct           
    for i = 1:n                   
      theta  = ThetaAll(i); 
      c = cosd(theta);
      s = sind(theta);
      T = [ c^2 s^2 2*c*s; 
          s^2 c^2 -2*c*s; 
          -c*s c*s (c^2 - s^2)];
         
      Qbar = inv(T)*Q*(inv(T))';
      
      if theta == CD
          if degradation == 1
              Qbar = 0;
          else 
              if degradation == 2
              E2 =0; 
              G12 = 0;
              Q = [ Q11 Q12 0; 
                  Q12 Q22 0; 
                  0 0 Q66];
              Qbar = inv(T)*Q*(inv(T))';
             end
          end
      end
      
      abar = T' * a;
      bbar = T' * b;
      
      
      A = A + Qbar*t;
      B = B + Qbar*t*zkbar(i); 
      D = D + Qbar * (t * zkbar(i)^2 +t^3 / 12);
      
      NT = NT + Qbar * abar * t * delT;
      MT = MT + Qbar * abar * t * zkbar(i) * delT;
       
      NH = NH + Qbar * bbar * t * delC;
      MH = MH + Qbar * bbar * t * zkbar(i) * delC;
       
    end
        
    Qbar;
    A;
    B;
    D;
    
    ABD = [A B; 
          B D];
    
    %Engineering properties in global axes
    invA = inv(A);
    Ex = 1/(h*invA(1,1));
    Ey = 1/(h*invA(2,2));
    nuxy = -(invA(1,2)/invA(1,1));
    nuyx = -(invA(1,2)/invA(2,2));
    Gxy  = 1/(h*invA(3,3));
    
    midssdelT = ABD \ [NT' MT']';
    midssdelM = ABD \ [NH' MH']';
    midstrain = ABD \ [ N';M']  ;
    
    
    for i = 1:n                              % finding stresses and strains 
    
      strainG(:,i) = midstrain(1:3,1) + zk(i)*midstrain(4:6,1);
      strainTG(:,i) = midssdelT(1:3,1) + zk(i)*midssdelT(4:6,1);
      strainMG(:,i) = midssdelM(1:3,1) + zk(i)*midssdelM(4:6,1);
      
      theta  = ThetaAll(i) ;
      c = cosd(theta) ;
      s = sind(theta) ;
      T = [ c^2 s^2 2*c*s; s^2 c^2 -2*c*s;-c*s c*s (c^2 - s^2)];
      R = [1 0 0;0 1 0;0 0 2];
      Qbar = inv(T) * Q * (inv(T))';
      
      if theta == CD
          if degradation == 1
              Qbar = 0;
          else 
              if degradation == 2
              E2 =0; G12 = 0;
              Q = [ Q11 Q12 0; Q12 Q22 0; 0 0 Q66];
              Qbar = inv(T) * Q * (inv(T))';
             end
          end
      end
      
      abar = T' * a;
      FTS(:,i) = delT*abar;
      
      bbar = T' * b ;
      FME(:,i) = delC*bbar;
      
      resstrainT(:,i) = strainTG(:,i) - FTS(:,i);
      resstrainM(:,i) = strainMG(:,i) - FME(:,i);
      resstressG(:,i) = Qbar * (resstrainT(:,i)+resstrainM(:,i));
      resstressL(:,i) = T * resstressG(:,i);
      
      strainL(:,i) = R*T*inv(R)*strainG(:,i); 
      strainML(:,i) = R*T*inv(R)*strainMG(:,i); 
      strainTL(:,i) = R*T*inv(R)*strainTG(:,i); 
      stressL(:,i) = Q*strainL(:,i);
      
      totstress(:,i) = stressL(:,i) + resstressL(:,i);
    end
    
    stressL;
    resstressG;
    resstressL;
    if k == 1
        totstress
        totstrain = strainL + strainML + strainTL
    end
    
    % finding strenghth ratio of all plies
    for j = 1:n
        if ThetaAll(j) ~= CD 
            
            if(totstress(1,j)>0)
               SRL(1,j) = totstress(1,j)/S(1);
            else
               SRL(1,j) = totstress(1,j)/S(2);
            end
      
            if(totstress(2,j)>0)
               SRT(1,j) = totstress(2,j)/S(3);
            else
                SRT(1,j) = totstress(2,j)/S(4);
            end
      
            SRS(1,j) = abs(totstress(3,j)/S(5));
            
        end
    end
    
    SR = [SRL;SRT;SRS];                     	%assembly
    
    [m,I] = max(SR(:));               %I is the index maximun Here tu can change the function to max or min
    [I_row, I_col] = ind2sub(size(SR),I) ;  % I_row is the row index and I_col is the column index
    
    CD(k) = ThetaAll(I_col);                 %failure ply angle 
    
    if(I_row ~= 3)                           %sigma due to mechanical load(Nx)
        if (totstress(I_row, I_col)>0)
            SM = S(2*I_row-1) - resstressL(I_row, I_col);
        else
            SM = S(2*I_row  ) - resstressL(I_row, I_col);
        end
    else
        SM = S(5) - resstressL(I_row, I_col);
    end
    
    SM;
    
    PF(k,1) = N(1)*SM/stressL(I_row, I_col)        % kth ply failure in N/m
    
    N(1) = PF(k,1);
    A = 0;
    B = 0;
    D = 0;

end
