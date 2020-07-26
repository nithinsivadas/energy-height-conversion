function y=double_gaussian(x,a,m,s1,s2)
    % Plots two gaussians, with one-half that can specify a rise time
    % the other half a decay time. 
    y = a;
    y = y.*((x<m).*exp(-0.5.*((x-m)./s1).^2)+(x>=m).*exp(-0.5.*((x-m)./s2).^2));
end