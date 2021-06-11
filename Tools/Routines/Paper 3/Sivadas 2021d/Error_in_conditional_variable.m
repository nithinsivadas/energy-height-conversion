%% Error in conditional variable
% Y - is the dependent variable, X- is the control variable
% E(Y|X) is a good estimate of Y
%%
xNoisedB = -20;
yNoisedB = -1;

%% True signal
xTrue = 0:0.1:10;
yTrue = dirac(xTrue-5);
idx = yTrue == Inf;
yTrue(idx) = 2;

figure;
subplot(2,3,1);

plot(xTrue,yTrue);
xlabel('X');
ylabel('Y');
ylim([-1, 3]);
xlim([-1, 11]);
title('True Signal');

%% Noisy signal (Noise in Y)

samples = 1000;
yNoisey = repmat(yTrue,samples,1) + wgn(samples, size(yTrue,2),yNoisedB);
xNoisey = repmat(xTrue,samples,1);

yNoisey1 = yNoisey(:);
xNoisey1 = xNoisey(:);

subplot(2,3,2);
scatter(xNoisey1, yNoisey1, 10, 'filled','MarkerFaceAlpha',0.2);
ylim([-5, 5]);
xlim([-1, 11]);
xlabel('X');
ylabel('Y^*');

title('Noise in Y');
%% Conditional expectation of noisy Y given X
[xindx,E] = discretize(xNoisey1,length(xTrue));

for i = 1:max(xindx)
   YgX(i) = nanmean(yNoisey1(xindx==i));
   X(i)   = nanmean(xNoisey1(xindx==i));
end

subplot(2,3,3);
plot(X,YgX);
xlim([-1, 11]);
ylim([-1, 3]);
ylabel('<Y^*|X>');
xlabel('X');

title('Conditional expectation of Y');

%% Noisy singal in Y and X

xNoisey2 = xNoisey1  + wgn(length(xNoisey1),1,xNoisedB);

subplot(2,3,5);
scatter(xNoisey2,yNoisey1,10,'filled','MarkerFaceAlpha',0.2);
ylim([-5, 5]);
xlim([-1, 11]);
xlabel('X^*');
ylabel('Y^*');

title('Noise in Y and X');

%% Conditional expectation of noisy Y given noisy X

[xindx2,E2] = discretize(xNoisey2(:),E);
for i = 1:max(xindx2)
   YgX2(i) = nanmean(yNoisey1(xindx2==i));
   X2(i)   = nanmean(xNoisey2(xindx2==i));
end

subplot(2,3,6);
plot(X2,YgX2);
xlim([-1, 11]);
ylim([-1, 3]);
ylabel('<Y^*|X^*>');
xlabel('X^*');

title({'Conditional expectation of Y','given a noisy X'});

%% Comparing Y, conditional expectation of Y, and that conditional on noisy X
subplot(2,3,4);

plot(xTrue,yTrue,'k');
hold on;
plot(X,YgX,'.r');
hold on;
plot(X2,YgX2,'b');

legend('Y','<Y^*|X>','<Y^*|X^*>');
xlabel('X');
ylim([-1, 3]);
xlim([-1, 11]);
title('Combined');


