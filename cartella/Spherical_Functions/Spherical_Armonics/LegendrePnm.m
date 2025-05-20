%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%% LEGENDRE ASSOCIATED FUNCTIONS %%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
% This code implements Legendre associated functions as a file function to
% be called from other analysis

%INPUT: -n mode under test
%       -m index sum for harmonics 
%       -Angle for function dependence Pnm(cos(theta)

%OUTPUT: -Pnm (cos(theta) for selected m and n

% Reference: On the Computation af Derivatives of Legendre Functions (W.Bosh 2000)
%author: Vincenzo Miranda
%Date: 20/05/24
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


function Pnm = LegendrePnm(n,m,angle)

Pnn=legendre(n,cos(angle));

if m>=0
Pnm(:,:,:)=Pnn(m+1,:,:,:);

elseif m<0
Pnm(:,:,:)=c2/c1*Pnn(abs(m)+1,:,:,:);
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%