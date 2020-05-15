clear all
close all
clc
i = 0;
sw1 = 0.000001;
a1 = zeros(1,int32((16/3)/sw1+1));
q1 = zeros(1,int32((16/3)/sw1+1));
for alpha = -5:sw1:1/3
 
    i=i+1; 
    
 fun = @(x) x./(sqrt(2/3*(x.^3-1)-2*alpha*(x-1)));
 
 q1(i) = integral(fun,1,0);
 
 a1(i) = alpha;
 
 lambda1(i) = q1(i)^2;
end
 i = 0;
sw2 = 0.1;
a2 = zeros(1,int32((16/3)/sw2+1));
q2 = zeros(1,int32((16/3)/sw2+1));


 for alpha = 1:sw2:5
 
    i=i+1; 
     
 fun = @(x) x./(sqrt(2/3*(x.^3-1)-2*alpha*(x-1)));
 
 q2(i) = integral(fun,1,0);
 
 a2(i) = alpha;
 
 lambda2(i) = q2(i)^2;
 
 end


for i = 1:length(a1)

    
fun = @(x) 1./(sqrt(2/3*(x.^3-1)-2*a1(i)*(x-1))*sqrt(lambda1(i)));

q1(i) = integral(fun,1,0);



end

 
 for i = 1:length(a2)
     
 fun = @(x) 1./(sqrt(2/3*(x.^3-1)-2*a2(i)*(x-1))*sqrt(lambda2(i)));
 
 q2(i) = integral(fun,1,0);
 
end
 
I1 = 1./q1;

I2 = -1./q2;


for i = 1:length(a1)

    n1(i) = 4*a1(i)*lambda1(i)/I1(i)^2;
    m1(i) = 4*lambda1(i)/I1(i)^2;

end


 for i = 1:length(a2)
 
    n2(i) = 4*a2(i)*lambda2(i)/I2(i)^2;%
    m2(i) = 4*lambda2(i)/I2(i)^2;
 end
 
m = [m1';m2'];
n = [n1';n2'];
Iges = [I1';I2'];
M = [m,Iges];
N = [m,Iges];

% Ausgabe
    dlmwrite('I über M 1508.2.txt',[M],'delimiter',' ','newline','pc');
    dlmwrite('I über N 1508.2.txt',[N],'delimiter',' ','newline','pc');


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

