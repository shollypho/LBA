sname = [19, 23, 24, 25, 26, 27, 28, 29, 30, 31, 32, 33, 527, 528, 529, 530, 533, 534];

for ss = 1:18%length(sname)
    clear D
    D = spm_eeg_load(sprintf('../s%d/devracbefffMs%d_FT_blk1.mat', sname(ss), sname(ss)));
    
    D = copy(D, sprintf('%s/s%d/%s',pwd, sname(ss), D.fname));
end