function [LL] = objfun(c,K,X)
%%
Beta=c(1:6,1);%Betas for the intercenpt,Medu,Famincome, Neutral, Rural and numsib
n(1,1)=c(7,1);%Unobserved heterogeneity terms for two types
n(2,1)=0;
pt1=c(8,1);%pi-tilde
%% Defining the probability of people belonging to type 1 (2 types)
p1=exp(pt1)/(1+exp(pt1));
%% Evaluating the probability of moving from grade x to y
Lg=zeros(1661,14,2);
Psi=c(9:22,1);
for j=1:1661
    for k=1:14
        for i=1:2
            Lg(j,k,i)=exp(X(j,:)*Beta+(k+6)*Psi(k,1)+n(i,1))/(1+exp(X(j,:)*Beta+(k+6)*Psi(k,1)+n(i,1)));
        end
    end
end
%% Evaluating the probablities of stopping after a certain grade
L1=ones(1661,2);
for j=1:1661
    if K(j,1)==20
        for i=1:2
        L1(j,i)=Lg(j,1,i)*Lg(j,2,i)*Lg(j,3,i)*Lg(j,4,i)*Lg(j,5,i)*Lg(j,6,i)*Lg(j,7,i)*Lg(j,8,i)*Lg(j,9,i)*Lg(j,10,i)*Lg(j,11,i)*Lg(j,12,i)*Lg(j,13,i)*Lg(j,14,i);
        end
    elseif K(j,1)==19
        for i=1:2
        L1(j,i)=Lg(j,1,i)*Lg(j,2,i)*Lg(j,3,i)*Lg(j,4,i)*Lg(j,5,i)*Lg(j,6,i)*Lg(j,7,i)*Lg(j,8,i)*Lg(j,9,i)*Lg(j,10,i)*Lg(j,11,i)*Lg(j,12,i)*Lg(j,13,i)*(1-Lg(j,14,i));
        end
    elseif K(j,1)==18
        for i=1:2
        L1(j,i)=Lg(j,1,i)*Lg(j,2,i)*Lg(j,3,i)*Lg(j,4,i)*Lg(j,5,i)*Lg(j,6,i)*Lg(j,7,i)*Lg(j,8,i)*Lg(j,9,i)*Lg(j,10,i)*Lg(j,11,i)*Lg(j,12,i)*(1-Lg(j,13,i));
        end
    elseif K(j,1)==17
        for i=1:2
        L1(j,i)=Lg(j,1,i)*Lg(j,2,i)*Lg(j,3,i)*Lg(j,4,i)*Lg(j,5,i)*Lg(j,6,i)*Lg(j,7,i)*Lg(j,8,i)*Lg(j,9,i)*Lg(j,10,i)*Lg(j,11,i)*(1-Lg(j,12,i));
        end
    elseif K(j,1)==16
        for i=1:2
        L1(j,i)=Lg(j,1,i)*Lg(j,2,i)*Lg(j,3,i)*Lg(j,4,i)*Lg(j,5,i)*Lg(j,6,i)*Lg(j,7,i)*Lg(j,8,i)*Lg(j,9,i)*Lg(j,10,i)*(1-Lg(j,11,i));
        end
    elseif K(j,1)==15
        for i=1:2
        L1(j,i)=Lg(j,1,i)*Lg(j,2,i)*Lg(j,3,i)*Lg(j,4,i)*Lg(j,5,i)*Lg(j,6,i)*Lg(j,7,i)*Lg(j,8,i)*Lg(j,9,i)*(1-Lg(j,10,i));
        end
    elseif K(j,1)==14
        for i=1:2
        L1(j,i)=Lg(j,1,i)*Lg(j,2,i)*Lg(j,3,i)*Lg(j,4,i)*Lg(j,5,i)*Lg(j,6,i)*Lg(j,7,i)*Lg(j,8,i)*(1-Lg(j,9,i));
        end
    elseif K(j,1)==13
        for i=1:2
        L1(j,i)=Lg(j,1,i)*Lg(j,2,i)*Lg(j,3,i)*Lg(j,4,i)*Lg(j,5,i)*Lg(j,6,i)*Lg(j,7,i)*(1-Lg(j,8,i));
        end
    elseif K(j,1)==12
        for i=1:2
        L1(j,i)=Lg(j,1,i)*Lg(j,2,i)*Lg(j,3,i)*Lg(j,4,i)*Lg(j,5,i)*Lg(j,6,i)*(1-Lg(j,7,i));
        end
    elseif K(j,1)==11
        for i=1:2
        L1(j,i)=Lg(j,1,i)*Lg(j,2,i)*Lg(j,3,i)*Lg(j,4,i)*Lg(j,5,i)*(1-Lg(j,6,i));
        end
    elseif K(j,1)==10
        for i=1:2
        L1(j,i)=Lg(j,1,i)*Lg(j,2,i)*Lg(j,3,i)*Lg(j,4,i)*(1-Lg(j,5,i));
        end
    elseif K(j,1)==9
        for i=1:2
        L1(j,i)=Lg(j,1,i)*Lg(j,2,i)*Lg(j,3,i)*(1-Lg(j,4,i));
        end
    elseif K(j,1)==8
        for i=1:2
        L1(j,i)=Lg(j,1,i)*Lg(j,2,i)*(1-Lg(j,3,i));
        end
    elseif K(j,1)==7
        for i=1:2
        L1(j,i)=Lg(j,1,i)*(1-Lg(j,2,i));
        end
    end
end
%% Calculating the contribution of each individual to the likelihood as a weighted average of the two types of people 
L1_type=zeros(1661,1);
L1_type(:,1)=p1*L1(:,1)+(1-p1)*L1(:,2);
%% Calculating the Logarithm of the likelihood for the sample of 1661 people
LL=0;
for i=1:1661
      LL=log(L1_type(i,1))+LL;
 %     LL=log(L1(i,1))+LL;
end
LL=-LL;

end