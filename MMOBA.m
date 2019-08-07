clear;%% use same random number and switch to equal allocation when pj=0
n0=5;
Mu=[0.5 5.5;1.9 4.2;2.8 3.3;3 3;3.9 2.1;4.3 1.8;4.6 1.5;3.8 6.3;4.8 5.5;5.2 5;5.9 4.1;6.3 3.8;6.7 7.2;7 7;7.9 6.1;9 9];
[r,cl]=size(Mu);tsig=2*ones(r,2);
for k=1:r 
    a=Mu(k,1)<Mu([1:k-1,k+1:r],1);
    b=Mu(k,2)<Mu([1:k-1,k+1:r],2);
    c=a|b;
    cn=sum(c);
    f0(k)=(cn==(r-1));%% if point k is on the pareto front, f0(k)=1,0 otherwise
end

jn=4000;   %%how many different values of budget will be tested
ct=zeros(3,jn);%%the times of the right choice
T1=zeros(1,jn);%%how many times algorithm used with each budget
T2=zeros(1,jn);
budgets=0;%%the smallest budget
budgeti=1;%%budget step
mre=1000; %% the maximum repitation
record=zeros(r,mre,jn);%% record in each iteration which alternative gets the sampling and the pj accordingly
for re=1:mre
    for i=1:r%++++++++++++++++++++++++++++++++++++++
        for j=1:cl
            sps(i,j,1:jn)=normrnd(Mu(i,j),tsig(i,j),jn,1);
        end
    end

%%%%%%%%%%%%%%%the initial sample mean value and variance%%%%%%%%%%%%%%%%%%
    xb0=zeros(r,2);sig0=zeros(r,2);
    for i=1:n0
        X(1:r,2*i-1:2*i)= (normrnd(Mu,tsig));
        xb0= (xb0+X(1:r,2*i-1:2*i));
    end
    xb0 = (xb0/n0);
    for i=1:n0
        sig0 = (sig0+(X(1:r,2*i-1:2*i)-xb0).^2);
    end
    sig0 = (sqrt(sig0/(n0-1)));
%%%%%%%%%%%%%%%%%%%%%%%%%%%%euqal allocation%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    xb=xb0;n=n0*ones(r,1);budget=0;
    
    for j=1:jn
        tbudget=budgets+(j-1)*budgeti;
        while(budget<tbudget)
            mn=mod(budget,r)+1;
            budget=budget+1;
            X=sps(mn,:,n(mn)-4);%++++++++++++++++++4 need to change when n0!=5
            xb(mn,1:2)= ((n(mn)*xb(mn,1:2)+X)/(n(mn)+1));
            n(mn)=n(mn)+1;
        end
        for k=1:r
            a=xb(k,1)<xb([1:k-1,k+1:r],1);
            b=xb(k,2)<xb([1:k-1,k+1:r],2);
            c=a|b;
            cn=sum(c);
            f1(k)=(cn==(r-1));
        end
        ct(1,j)=ct(1,j)+(sum(f1==f0)==r);
    end

%%%%%%%%%%%%%%%%allocation by probability t-distribution and tau=10 and to equal%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    xb=xb0; sig=sig0;n=n0*ones(r,1);budget=0;
    
    for k=1:r
        a=xb(k,1)<xb([1:k-1,k+1:r],1);
        b=xb(k,2)<xb([1:k-1,k+1:r],2);
        c=a|b;
        cn=sum(c);
        f1(k)=(cn==(r-1));%%if point k is on the pareto front based on five observations, f1(k)=1,0 otherwise
    end
    for j=1:jn
        tbudget=budgets+(j-1)*budgeti;
        while(budget<tbudget)
            pj= (paretot1(xb,sig,f1,r,n));
            record(:,re,j)=pj;%%record in each iteration which alternative gets the sampling and the pj accordingly
            if sum(pj)~=0
                [m,mn]=max(pj);%%m is the largest number of pj and mn is the correspoding alternative's sequence number 
                %[m,mn]= min(pj); 
                T1(1,j)=T1(1,j)+1;%%how many times use algorithm
            else
                pj=(paretot(xb,sig,f1,r,n));%% switch to tau=10
                if sum(pj)~=0
                    [m,mn]=max(pj);
                else
                    T2(1,j)=T2(1,j)+1;
                    mn=mod(budget,r)+1;%% switch twice to equal
                end
            end
            
            X=sps(mn,:,n(mn)-4);%++++++++++++++++++++++++++++++++++++++
            sig(mn,1:2)=(sqrt((n(mn)-1)/n(mn)*sig(mn,1:2).^2+1/(n(mn)+1)*(X-xb(mn,1:2)).^2));
            xb(mn,1:2)=((n(mn)*xb(mn,1:2)+X)/(n(mn)+1));
            n(mn)=n(mn)+1;
            budget=budget+1;
            
            for k=1:r 
                a=xb(k,1)<xb([1:k-1,k+1:r],1);
                b=xb(k,2)<xb([1:k-1,k+1:r],2);
                c=a|b;
                cn=sum(c);
                f1(k)=(cn==(r-1));
            end 
        end
        ct(2,j)=ct(2,j)+(sum(f1==f0)==r);
    end
    fprintf('re=%d\n',re);
end
x=budgets+[0:(jn-1)]*budgeti;
ct1=ct(1,1:jn);
ct2=ct(2,1:jn);

cp1=ct1/mre;
cp2=ct2/mre;
figure
plot(x,cp1,'-r',x,cp2,'-.b')
legend('Equal','M-MOBA')
xlabel('budget')
ylabel('P{CS}')
title('Probability of Correct Selection');
