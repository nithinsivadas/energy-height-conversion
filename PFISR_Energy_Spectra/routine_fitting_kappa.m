%%  Routine to use lsqcurvefit to fit data_thm.E to a kappa distribution
%  Not successfull 5/18 :(


m = (9.11*10^-31);

x =[70000,10^8,20];
lb=[5000,10^5,2];
ub=[100000,10^9,170];
Eb=x(1);
T =x(2);
k =x(3);
n =3*10^3;

xdata = data_thm.ebin(1,:);
E = xdata;
ydata = data_thm.E(1,:);


[x, resnorm, residual, existflag, output]=lsqcurvefit(@test_fn,x,xdata,ydata,lb,ub,optimoptions(@lsqcurvefit,'Display','iter'));

hold on;
j = kappa_j(Eb,T,k,n,m,E);

figure; plot(xdata,ydata); hold on; plot(xdata,j,'r');
legend('Measurement','Fit');

display(['E_b = ',num2str(Eb,'%10.2e\n'),' eV ; T = ',num2str(T,'%10.2e\n'),' K ; k = ',num2str(k,'%10.2e\n'),' ; n = ',num2str(n,'%10.2e\n'),' cm^-3 ; resnorm = ',num2str(resnorm,'%10.2e\n')]); 
