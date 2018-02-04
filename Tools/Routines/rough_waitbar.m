hWait = waitbar(0);
steps = 1000;
for step = 1:steps
    % computations take place here
    custom_waitbar(h,step,steps,'Step 5/5: Generating Video');
end
close(hWait)