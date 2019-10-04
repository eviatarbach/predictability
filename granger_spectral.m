function [] = granger_spectral(offset, length)

nobs = 14235;
fres = 300;

cd mvgc_v1.0

startup

cd ..

sst = importdata('data/sst01_365.mat');
vort = importdata('data/vort01_365.mat');
div = importdata('data/div01_365.mat');
sp = importdata('data/sp01_365.mat');
temp = importdata('data/temp01_365.mat');
q = importdata('data/q01_365.mat');
times = importdata('data/times.mat');

f_atmos_to_sst = NaN(length, fres + 1);
f_sst_to_atmos = NaN(length, fres + 1);

'offset', offset

for cell = 1:length
    'progress', cell/length
    if offset+cell-1 > 88838
        continue
    end

    sst_ts = sst(:, offset+cell-1)';
    vort_ts = vort(:, offset+cell-1)';
    div_ts = div(:, offset+cell-1)';
    sp_ts = sp(:, offset+cell-1)';
    temp_ts = temp(:, offset+cell-1)';
    q_ts = q(:, offset+cell-1)';

    if (sst_ts(1) == 9999) || (vort_ts(1) == 9999) || (div_ts(1) == 9999) || (sp_ts(1) == 9999) || (temp_ts(1) == 9999) || (q_ts(1) == 9999)
        continue
    end

    moAIC = times(offset+cell-1);

    X = [sst_ts; vort_ts; div_ts; sp_ts; temp_ts; q_ts];

    [A, SIG] = tsdata_to_var(X, moAIC);
    assert(~isbad(A), 'VAR estimation failed');
    [G, info] = var_to_autocov(A, SIG);
    acerr = var_info(info, false);

    if acerr
        continue
    end

    f_atmos_cause = autocov_to_smvgc(G, [1], [2 3 4 5 6], fres);  % atmos -> sst

    f_sst_cause = autocov_to_smvgc(G, [2 3 4 5 6], [1], fres);  % sst -> atmos

    f_atmos_to_sst(cell, :) = f_atmos_cause;
    f_sst_to_atmos(cell, :) = f_sst_cause;
end

save(['data_atmos/f_atmos_to_sst_' num2str(offset) '.mat'], 'f_atmos_to_sst');
save(['data_atmos/f_sst_to_atmos_' num2str(offset) '.mat'], 'f_sst_to_atmos');
