function [derPnm] = Legendre_derPnm(n,m,angle)

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Function for correct evaluation of Legendre associated function
% derivative fixed degree n and order m.

% Reference: On the Computation af Derivatives of Legendre Functions W.
%            Bosch 2000
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


c2=factorial(n-m)/factorial(n+m); 

% All legendre associated functions for m= 0,1..n
Pnm=legendre(n,cos(angle));
% Definition of Legendre Associated Function for m<0
Pnm_negative=c2*legendre(n,cos(angle));

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%Special case for m=0
if m==0
  Pn1(:,:,:)=Pnm(m+2,:,:,:);
  derPnm(:,:,:)=(-Pn1(:,:,:));

%Evaluation of Derivative Pnm for m>0 
elseif m>0 

  if m==n
   P_nm_minus(:,:,:)=Pnm(m,:,:,:);  %Pn m-1   
   derPnm=n*P_nm_minus ; %Recursive relation for Der Pnm     

  elseif m<n 
  P_nm_plus(:,:,:)=Pnm(m+2,:,:,:); %Pn m+1 
  P_nm_minus(:,:,:)=Pnm(m,:,:,:);  %Pn m-1   
  derPnm=0.5*( (n+m)*(n-m+1)*P_nm_minus - P_nm_plus ); %Recursive relation for Der Pnm 
  end

%Evaluation of Derivative Pnm for m<0 
elseif m<0 

  if abs(m)==n
   P_nm_minus(:,:,:)=Pnm_negative(abs(m),:,:,:);  %Pn m-1   
   derPnm=n*P_nm_minus ; %Recursive relation for Der Pnm     

  elseif abs(m)<n 
  P_nm_plus(:,:,:)=Pnm_negative(abs(m)+2,:,:,:); %Pn m+1 
  P_nm_minus(:,:,:)=Pnm_negative(abs(m),:,:,:);  %Pn m-1   
  derPnm=0.5*( (n+m)*(n-m+1)*P_nm_minus - P_nm_plus ); %Recursive relation for Der Pnm 
  end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
 end %end function
 %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%