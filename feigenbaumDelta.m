% Compute the Feigenbaum delta
% Store approximate values in the row vector delta for assessment, where length(delta)= num_doublings and 
% delta(2:num_doublings) are computed from the algorithm described in Lectures 21-23.
num_doublings=11; delta=zeros(1,num_doublings); delta(1)=5;
% Write your code here
x_0=0.5;
xp_0=0.0;
m=zeros(1,num_doublings+1);
m(1)=2;
m(2)=1+sqrt(5);
for n=2:num_doublings
    mu_0=m(n)+(m(n)-m(n-1))/delta(n-1);
    mu=mu_0;
    %while mu/mu_0 < 1.0e-26
    for i=1:10 %Newton's iteration
        x=x_0;
        xp=xp_0;
        for j=1:2^n
            xp=x*(1-x)+mu*xp*(1-2*x);
            x=mu*x*(1-x);
        end
        mu = mu - (x-1/2)/xp;
    end
    m(n+1)=mu;
    delta(n)=(m(n)-m(n-1))/(m(n+1)-m(n));
end


% Output your results
fprintf('n        delta(n)\n');
for n=1:num_doublings
    fprintf('%2g %18.15f\n',n,delta(n));
end

