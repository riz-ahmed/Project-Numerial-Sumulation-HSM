clear all
close all
clc
i = 0;

for alpha = -5:0.001:1/3

   i=i+1; 
    
fun = @(x) x./(sqrt(2/3*(x.^3-1)-2*alpha*(x-1)));

q = integral(fun,1,0);

a(i) = alpha;

lambda(i) = 4*q^2;

end

i = 0;

for alpha = 1+0.001:0.001:5

   i=i+1; 
    
fun = @(x) x./(sqrt(2/3*(x.^3-1)-2*alpha*(x-1)));

q2 = integral(fun,1,0);

a2(i) = alpha;

lambda2(i) = 4*q2^2;

end


plot(a,lambda)
grid on
hold on
plot(a2,lambda2)
axis([-5 5 -10 10])
xlabel('\alpha')
ylabel('\lambda')


alpha_lambda = lambda.*a;
alpha_lambda2 = lambda2.*a2;


figure
plot(alpha_lambda,lambda)
grid on
hold on 
plot(alpha_lambda2,lambda2)
axis([-5 30 -10 30])
xlabel('\alpha lambda')
ylabel('\lambda')
