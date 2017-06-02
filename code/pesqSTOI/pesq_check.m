% run "pesq_check" to test pesqITU
[pp,ppm] = fileparts(mfilename('fullpath'));
x=pesqITU(8000,fullfile(pp,'refcheck.wav'),fullfile(pp,'degcheck.wav'));
fprintf('\nTest for PESQ\nPESQ obtained is %.3f, expected range [2.187, 2.287]\n',x);
if abs(x-2.237)<=0.05
    fprintf('*** PASSED ***\n\n');
else
    fprintf('*** FAILED ***\n\n');
end