function [accuracy] = accuracy_of_star_calibration(calStarsAzEl,realStarsAzEl,...
                                                    azCal,elCal)
%ACCURACY_OF_STAR_CALIBRATION Calculates the accuracy of star calibration
%   Input
%   calStarsAzEl  - Calibrated [azimuth, elevation] in degrees of stars 
%                   extracted from image 
%   realStarsAzEl - [azimuth, elevation] in degrees of stars calculated 
%                   using a star chart 
%   azCal         - Calibrated azimuth values of all the pixels of the image
%   elCal         - Calibrated elevation values of all the pixels of the image

%   Output
%   accuracy
%          ->pixel.mean  - Mean of the distance between the real and calibrated star points
%          ->pixel.std   - its std deviation
%          ->angle.mean  - the mean of the angular distance
%          ->angle.std   - its std deviation

minDistance = min(pdist2(calStarsAzEl,realStarsAzEl),[],2);
sigma = std(minDistance);
medianDistance = median(minDistance);
matchedStars = minDistance < medianDistance + sigma;

accuracy.angle.mean = mean(minDistance(matchedStars));
accuracy.angle.std = std(minDistance(matchedStars));

imageSize = size(azCal,1);
matrixIndex = reshape((1:1:(imageSize.^2)),imageSize,imageSize);

nonNanIndx = ~isnan(azCal) & ~isnan(elCal); 
F = scatteredInterpolant(azCal(nonNanIndx),elCal(nonNanIndx),matrixIndex(nonNanIndx),'nearest','none');

calStarPixelNumber = F(calStarsAzEl(:,1),calStarsAzEl(:,2));
realStarPixelNumber = F(realStarsAzEl(:,1),realStarsAzEl(:,2));

[calStarPixelx,calStarPixely] = ind2sub(size(azCal),calStarPixelNumber);
[realStarPixelx,realStarPixely] = ind2sub(size(azCal),realStarPixelNumber);

minDistancePixels = min(pdist2([calStarPixelx,calStarPixely],...
[realStarPixelx,realStarPixely]),[],2);
sigmaP = std(minDistancePixels);
medianDistanceP = median(minDistancePixels);
matchedStarsP = minDistancePixels < medianDistanceP + sigmaP;

accuracy.pixel.mean = mean(minDistancePixels(matchedStarsP));
accuracy.pixel.std = std(minDistancePixels(matchedStarsP));

end

