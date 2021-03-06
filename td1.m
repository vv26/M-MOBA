function pj=td1(xab,yab,xb,sig,k,n)
kp=n(k)*(n(k)+1)./(1*sig(k,1:2).^2);
mu=xb(k,1:2);
nu=n(k)-1;
txab=(xab-mu(1))*sqrt(kp(1));
tyab=(yab-mu(2))*sqrt(kp(2));
pj1=tcdf(txab(2),nu)-tcdf(txab(1),nu);
pj2=tcdf(tyab(2),nu)-tcdf(tyab(1),nu);
pj=pj1*pj2;