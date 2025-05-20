function T=Tau(n,m,angle)

    ang = reshape(angle, [], 1)';
    Pm=legendre(n,cos(ang));

%Special case for m=0
if m==0
  Pn1=Pm(m+2,:);      %Pn m+1 
  T=(-Pn1);

%Evaluation of Derivative Pnm for m>0 
elseif m==n
P_nm_minus=Pm(m,:);  %Pn m-1   
T=n*P_nm_minus ; %Recursive relation for Der Pnm     

  elseif m<n 
  P_nm_plus=Pm(m+2,:); %Pn m+1 
  P_nm_minus=Pm(m,:);  %Pn m-1   
  T=0.5*( (n+m)*(n-m+1)*P_nm_minus - P_nm_plus ); %Recursive relation for Der Pnm 
end
    T=reshape(T,size(angle));
end





