clear all
close all
clc
i = 0;


 
sw = 0.0001;
l = -5;
r = 1/3;
schritte = round((abs(l-r))/sw+1);
q = zeros(2,schritte);
 fun = @(x,alpha) 2*x*alpha;     %x./(sqrt(2/3*(x.^3-1)-2*alpha*(x-1)));
 for alpha =l:sw:r
     i=i+1;
     q(1,i) = alpha;
     q(2,i) = integral(@(x)fun(x,alpha),1,0);
 end
     
 lambda = q(2,:).^2;
 
I1 = 1./q(2,:);





    n = 4*q(1,i)*lambda./I1.^2;
    m = 4*lambda./I1.^2;


 %lambda = q.^2;


% i = 0;
% 
% 
%  for alpha = 1:0.00001:1.1
%  
%     i=i+1; 
%      
%  fun = @(x) x./(sqrt(2/3*(x.^3-1)-2*alpha*(x-1)));
%  
%  q2(i) = integral(fun,1,0);
%  
%  a2(i) = alpha;
%  
%  lambda2(i) = q2(i)^2;
%  
%  end
% 
% 
% for i = 1:length(a1)
% 
%     
% fun = @(x) 1./(sqrt(2/3*(x.^3-1)-2*a1(i)*(x-1))*sqrt(lambda1(i)));
% 
% q1(i) = integral(fun,1,0);
% 
% 
% 
% end
% 
%  
%  for i = 1:length(a2)
%      
%  fun = @(x) 1./(sqrt(2/3*(x.^3-1)-2*a2(i)*(x-1))*sqrt(lambda2(i)));
%  
%  q2(i) = integral(fun,1,0);
%  
% end
%  
% I1 = 1./q1;
% 
% I2 = -1./q2;
% 
% 
% for i = 1:length(a1)
% 
%     n1(i) = 4*a1(i)*lambda1(i)/I1(i)^2;
%     m1(i) = 4*lambda1(i)/I1(i)^2;
% 
% end
% 
% 
%  for i = 1:length(a2)
%  
%     n2(i) = 4*a2(i)*lambda2(i)/I2(i)^2;%
%     m2(i) = 4*lambda2(i)/I2(i)^2;
%  end
%  
% m = [m2';m1'];
% n = [n2';n1'];
% Iges = [I2';I1'];
% M = [m,Iges];
% N = [n,Iges];
% 
% % Ausgabe
%     dlmwrite('I über M 1508.10.txt',[M],'delimiter',' ','newline','pc');
%     dlmwrite('I über N 1508.10.txt',[N],'delimiter',' ','newline','pc');
% 
% 
% % figure
% % plot(n,m,'x-')
% % hold on
% % plot(n2,m2)
% % grid on
% % xlabel('n')
% % ylabel('m')
% % 
% % 
% % figure 
% % subplot(1,2,1)
% % plot(m,I,'r')
% % hold on
% % plot(m2,I2,'b')
% % grid on
% % subplot(1,2,2)
% % plot(n,I,'r')
% % hold on
% % plot(n2,I2,'bx-')
% % grid on
% 

