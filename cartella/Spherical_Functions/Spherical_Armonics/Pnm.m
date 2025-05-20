function P=Pnm(n,m,angle)

    ang = reshape(angle, [], 1)';
    Pm=legendre(n,cos(ang));
    P=squeeze(Pm(m+1,:));
    P=reshape(P,size(angle));

end