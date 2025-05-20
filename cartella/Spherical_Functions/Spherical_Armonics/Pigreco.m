function Pig=Pigreco(n,m,angle)

   
    ang = reshape(angle, [], 1)';
    pig0=0*ones(1,length(ang));
    pig1=1*ones(1,length(ang));

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% 
        if m==0
            Pig=0*ones(1,length(ang));
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        elseif m==1

        Pig_aux=zeros(n+1,length(ang));
        Pig_aux(2,:)=pig1;
        Pig_aux(1,:)=pig0;

        for q=2:n
        Pig_aux(q+1,:)= (2*q-1)/(q-1)*cos(ang).*Pig_aux(q,:) - q/(q-1)*Pig_aux(q-1,:);
        end
        
        Pig=Pig_aux(n+1,:);
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        else 
      Pm=legendre(n,cos(ang));
      P=squeeze(Pm(m+1,:));
      Pig =m*P./sin(ang); % Funzione settoriale   
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

        end
 Pig(isnan(Pig)) = 0;       
 Pig=reshape(Pig,size(angle));

end