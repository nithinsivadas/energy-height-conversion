function [fitresult, gof] = createFit_q(h, q_1)
%CREATEFIT1(H,Q_1)
%  Create a fit.
%
%  Data for 'untitled fit 1' fit:
%      X Input : h
%      Y Output: q_1
%  Output:
%      fitresult : a fit object representing the fit.
%      gof : structure with goodness-of fit info.
%
%  See also FIT, CFIT, SFIT.

%  Auto-generated by MATLAB on 17-May-2016 19:59:40
[xData, yData] = prepareCurveData( h, q_1 );

% Set up fittype and options.
ft = fittype( 'gauss1' );
opts = fitoptions( 'Method', 'NonlinearLeastSquares' );
opts.Display = 'Off';
opts.Lower = [-Inf -Inf 0];
opts.StartPoint = [123885.068304582 100.792799775062 14.7990129672142];

% Fit model to data.
[fitresult, gof] = fit( xData, yData, ft, opts );
end


