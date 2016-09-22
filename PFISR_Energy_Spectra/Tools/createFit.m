function [fitresult, gof] = createFit(log_ebin, eflux)
%CREATEFIT(LOG_ENERGY,TEST_ENERGY)
%  Create a fit.
%
%  Data for 'untitled fit 1' fit:
%      X Input : log_energy
%      Y Output: test_energy
%  Output:
%      fitresult : a fit object representing the fit.
%      gof : structure with goodness-of fit info.
%
%  See also FIT, CFIT, SFIT.

%  Auto-generated by MATLAB on 17-May-2016 19:42:44


%% Fit: 'untitled fit 1'.
[xData, yData] = prepareCurveData( log_ebin, eflux );

% Set up fittype and options.
ft = fittype( 'gauss1' );
opts = fitoptions( 'Method', 'NonlinearLeastSquares' );
% opts.Display = 'Off';
% opts.Lower = [-Inf -Inf 0];
% opts.StartPoint = [77282128.6 3.82411802045229 0.387844248373923];

% Fit model to data.
[fitresult, gof] = fit( xData, yData, ft, opts );

% Plot fit with data.
% figure( 'Name', 'untitled fit 1' );
% h = plot( fitresult, xData, yData );
% legend( h, 'test_energy vs. log_energy', 'untitled fit 1', 'Location', 'NorthEast' );
% % Label axes
% xlabel log_energy
% ylabel test_energy
% grid on

