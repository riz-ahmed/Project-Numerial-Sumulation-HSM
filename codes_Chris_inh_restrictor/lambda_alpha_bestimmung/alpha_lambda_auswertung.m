clear all
close all
clc
i = 0;

for alpha = -5:0.00001:0.24999
 
    i=i+1; 
    
 fun = @(x) x./(sqrt(2/3*(x.^3-1)-2*alpha*(x-1)));
 
 q1(i) = integral(fun,1,0);
 
 a1(i) = alpha;
 
 lambda1(i) = q1(i)^2;
end
 i = 0;

 for alpha = 1/4:0.000001:1/3
 
    i=i+1; 
     
 fun = @(x) x./(sqrt(2/3*(x.^3-1)-2*alpha*(x-1)));
 
 q2(i) = integral(fun,1,0);
 
 a2(i) = alpha;
 
 lambda(i) = q2(i)^2;
 
 end
 i = 0;

for alpha = 1:0.0001:1.000001

   i=i+1; 
    
fun = @(x) x./(sqrt(2/3*(x.^3-1)-2*alpha*(x-1)));

q(i) = integral(fun,1,0);

a(i) = alpha;

lambda(i) = q(i)^2;

end
 i = 0;
 
for alpha = 1.00101:0.000001:1.1
 
    i=i+1; 
     
 fun = @(x) x./(sqrt(2/3*(x.^3-1)-2*alpha*(x-1)));
 
 q4(i) = integral(fun,1,0);
 
 a4(i) = alpha;
 
 lambda4(i) = q2(i)^2;
 
 end

i = 0;

 for alpha = 1.10001:0.00001:5
 
    i=i+1; 
     
 fun = @(x) x./(sqrt(2/3*(x.^3-1)-2*alpha*(x-1)));
 
 q5(i) = integral(fun,1,0);
 
 a5(i) = alpha;
 
 lambda5(i) = q2(i)^2;
 
 end


for i = 1:length(a)

    
fun = @(x) 1./(sqrt(2/3*(x.^3-1)-2*a(i)*(x-1))*sqrt(lambda(i)));

q(i) = integral(fun,1,0);



end

 
 for i = 1:length(a2)
     
 fun = @(x) 1./(sqrt(2/3*(x.^3-1)-2*a2(i)*(x-1))*sqrt(lambda2(i)));
 
 q2(i) = integral(fun,1,0);
 
end
 
 for i = 1:length(a1)
 
     
 fun = @(x) 1./(sqrt(2/3*(x.^3-1)-2*a(i)*(x-1))*sqrt(lambda(i)));
 
 q1(i) = integral(fun,1,0);
 
 

 end
 
 
 for i = 1:length(a2)
     
 fun = @(x) 1./(sqrt(2/3*(x.^3-1)-2*a2(i)*(x-1))*sqrt(lambda2(i)));
 
 q2(i) = integral(fun,1,0);
 
 end
for i = 1:length(a1)
 
     
 fun = @(x) 1./(sqrt(2/3*(x.^3-1)-2*a(i)*(x-1))*sqrt(lambda(i)));
 
 q1(i) = integral(fun,1,0);
 
 
 
 end
 
 
 for i = 1:length(a2)
     
 fun = @(x) 1./(sqrt(2/3*(x.^3-1)-2*a2(i)*(x-1))*sqrt(lambda2(i)));
 
 q2(i) = integral(fun,1,0);
 
 end

% plot(a,lambda)
% grid on
% hold on
% plot(a2,lambda2)
% axis([-5 5 -10 10])
% xlabel('\alpha')
% ylabel('\lambda')


alpha_lambda = lambda.*a;
% alpha_lambda2 = lambda2.*a2;
% 
% 
% figure
% plot(alpha_lambda,lambda)
% grid on
% hold on 
% plot(alpha_lambda2,lambda2)
% axis([-5 30 -10 30])
% xlabel('\alpha \lambda')
% ylabel('\lambda')
% 
% 
% figure 
% plot(a.*lambda,1./q)
% hold on
% plot(a2.*lambda2,-1./q2)

I = 1./q;

% I2 = -1./q2;


for i = 1:length(a)

    n(i) = 4*a(i)*lambda(i)/I(i)^2;


end


 for i = 1:length(a2)
 
     n2(i) = 4*a2(i)*lambda2(i)/I2(i)^2;%

 end

% figure
% plot(n,I)
% hold on
% plot(n2,I2)
% grid on
% xlabel('n')
% ylabel('I')
% N = [n';n2'];
% Iges = [I';I2'];
M = [n',I'];
dlmwrite('I über N 1308.txt',[M],'delimiter',' ','newline','pc');


for i = 1:length(a)

    m(i) = 4*lambda(i)/I(i)^2;

end
M2 = [m',I'];
dlmwrite('I über m 1308.txt',[M2],'delimiter',' ','newline','pc');

 for i = 1:length(a2)
 
     m2(i) = 4*lambda2(i)/I2(i)^2;
 
 end

% figure
% plot(n,m,'x-')
% hold on
% plot(n2,m2)
% grid on
% xlabel('n')
% ylabel('m')
% 
% 
% figure 
% subplot(1,2,1)
% plot(m,I,'r')
% hold on
% plot(m2,I2,'b')
% grid on
% subplot(1,2,2)
% plot(n,I,'r')
% hold on
% plot(n2,I2,'bx-')
% grid on


