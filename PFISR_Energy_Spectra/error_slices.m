% Error slices!

alpha_0 = [5000,10000,10,1000];

x_0 = 100:100:20000;

y_0 = kappa_j(alpha_0,x_0);

a1=10:100:10000;      % Eb Mean Energy
a2=1000:1000:20000;  % T Temperature
a3=2:1:15;         % k Kappa
a4=100:100:2000;     % N Density of plasma

% Residue
for i=1:1:length(a1)
    for j=1:1:length(a2)
        for k=1:1:length(a3)
            for m=1:1:length(a4)
                R(i,j,k,m) = sum((kappa_j([a1(i),a2(j),a3(k),a4(m)],x_0)-y_0).^2);
            end
        end
    end
end

figure;
surf(squeeze(R(7,:,:,1)));