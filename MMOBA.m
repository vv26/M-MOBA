clear;%% use same random number and switch to equal allocation when pj=0
n0=5;
Mu=[0.5 5.5;1.9 4.2;2.8 3.3;3 3;3.9 2.1;4.3 1.8;4.6 1.5;3.8 6.3;4.8 5.5;5.2 5;5.9 4.1;6.3 3.8;6.7 7.2;7 7;7.9 6.1;9 9];
%Mu=[1 2;3 1;5 5];
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
cp3=[0.473	0.489	0.488	0.484	0.509	0.512	0.539	0.561	0.553	0.561	0.56	0.562	0.554	0.551	0.563	0.567	0.57	0.58	0.582	0.599	0.595	0.597	0.604	0.605	0.602	0.615	0.611	0.625	0.63	0.648	0.653	0.66	0.649	0.67	0.676	0.669	0.675	0.684	0.689	0.691	0.7	0.693	0.702	0.695	0.705	0.71	0.708	0.713	0.726	0.721	0.728	0.743	0.738	0.745	0.749	0.76	0.76	0.765	0.756	0.755	0.757	0.767	0.766	0.771	0.768	0.773	0.772	0.777	0.775	0.773	0.779	0.776	0.778	0.792	0.784	0.779	0.788	0.797	0.797	0.802	0.806	0.802	0.807	0.81	0.807	0.814	0.823	0.829	0.829	0.83	0.837	0.834	0.838	0.835	0.824	0.825	0.824	0.826	0.831	0.83	0.83	0.832	0.83	0.831	0.837	0.841	0.844	0.844	0.842	0.844	0.838	0.843	0.841	0.851	0.852	0.85	0.853	0.853	0.864	0.859	0.86	0.858	0.862	0.865	0.865	0.858	0.86	0.864	0.865	0.865	0.871	0.871	0.87	0.869	0.877	0.875	0.88	0.884	0.887	0.887	0.892	0.888	0.89	0.889	0.886	0.891	0.897	0.9	0.896	0.898	0.904	0.904	0.903	0.907	0.906	0.911	0.91	0.905	0.905	0.904	0.901	0.904	0.904	0.906	0.906	0.907	0.905	0.908	0.907	0.909	0.914	0.917	0.912	0.915	0.915	0.917	0.919	0.921	0.92	0.914	0.916	0.912	0.918	0.916	0.913	0.911	0.914	0.915	0.915	0.92	0.919	0.918	0.921	0.921	0.924	0.927	0.92	0.92	0.923	0.925	0.925	0.924	0.927	0.922	0.928	0.928	0.928	0.932	0.934	0.936	0.932	0.935	0.937	0.934	0.933	0.931	0.932	0.935	0.938	0.934	0.939	0.941	0.937	0.939	0.94	0.941	0.94	0.94	0.942	0.941	0.937	0.936	0.939	0.939	0.937	0.939	0.941	0.94	0.941	0.941	0.939	0.938	0.938	0.939	0.939	0.943	0.941	0.941	0.942	0.941	0.939	0.937	0.938	0.937	0.94	0.94	0.938	0.938	0.94	0.939	0.941	0.943	0.943	0.944	0.947	0.945	0.946	0.942	0.941	0.943	0.944	0.944	0.946	0.946	0.954	0.95	0.947	0.948	0.949	0.948	0.95	0.951	0.954	0.952	0.953	0.951	0.951	0.952	0.953	0.951	0.951	0.953	0.956	0.954	0.956	0.958	0.959	0.959	0.96	0.957	0.959	0.959	0.961	0.962	0.963	0.961	0.962	0.962	0.962	0.964	0.964	0.965	0.967	0.968	0.968	0.97	0.966	0.965	0.96	0.963	0.964	0.962	0.965	0.965	0.962	0.964	0.964	0.963	0.965	0.966	0.966	0.967	0.964	0.966	0.965	0.964	0.964	0.965	0.967	0.966	0.967	0.966	0.963	0.967	0.966	0.967	0.967	0.965	0.966	0.967	0.967	0.965	0.967	0.965	0.968	0.967	0.967	0.967	0.968	0.971	0.971	0.968	0.97	0.97	0.969	0.969	0.971	0.971	0.97	0.967	0.968	0.97	0.971	0.969	0.969	0.97	0.969	0.97	0.972	0.972	0.971	0.97	0.972	0.971	0.973	0.974	0.973	0.974	0.975	0.975	0.974	0.973	0.976	0.976	0.977	0.976	0.973	0.973	0.974	0.975	0.973	0.973	0.972	0.973	0.976	0.974	0.976	0.976	0.976	0.975	0.978	0.977	0.978	0.979	0.98	0.979	0.98	0.979	0.98	0.982	0.98	0.98	0.979	0.979	0.979	0.98	0.979	0.98	0.979	0.98	0.979	0.978	0.98	0.981	0.982	0.981	0.98	0.979	0.98	0.978	0.981	0.98	0.979	0.98	0.977	0.977	0.978	0.98	0.981	0.982	0.983	0.982	0.984	0.984	0.985	0.984	0.986	0.985	0.986	0.984	0.985	0.985	0.986	0.986	0.986	0.985	0.987	0.986	0.986	0.986	0.986	0.987	0.988	0.987	0.988	0.988	0.989	0.989	0.989	0.986	0.987	0.99	0.99	0.989	0.989	0.99	0.987	0.986	0.986	0.985	0.985	0.985	0.988	0.986	0.985	0.987	0.988	0.986	0.987	0.986]
figure
plot(x,cp1,'-r',x,cp2,'-.b')
legend('Equal','M-MOBA')
xlabel('budget')
ylabel('P{CS}')
title('Probability of Correct Selection');