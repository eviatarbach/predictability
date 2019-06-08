function [] = granger_seasonal_multiple(offset, length)

cd mvgc_v1.0

startup

cd ..

sst_winter = importdata('data/sst_winter.mat');
vort_winter = importdata('data/vort_winter.mat');
div_winter = importdata('data/div_winter.mat');
q_winter = importdata('data/q_winter.mat');
temp_winter = importdata('data/temp_winter.mat');
sp_winter = importdata('data/sp_winter.mat');

sst_spring = importdata('data/sst_spring.mat');
vort_spring = importdata('data/vort_spring.mat');
div_spring = importdata('data/div_spring.mat');
q_spring = importdata('data/q_spring.mat');
temp_spring = importdata('data/temp_spring.mat');
sp_spring = importdata('data/sp_spring.mat');

sst_summer = importdata('data/sst_summer.mat');
vort_summer = importdata('data/vort_summer.mat');
div_summer = importdata('data/div_summer.mat');
q_summer = importdata('data/q_summer.mat');
temp_summer = importdata('data/temp_summer.mat');
sp_summer = importdata('data/sp_summer.mat');

sst_fall = importdata('data/sst_fall.mat');
vort_fall = importdata('data/vort_fall.mat');
div_fall = importdata('data/div_fall.mat');
q_fall = importdata('data/q_fall.mat');
temp_fall = importdata('data/temp_fall.mat');
sp_fall = importdata('data/sp_fall.mat');

sst = {sst_winter, sst_spring, sst_summer, sst_fall};
vort = {vort_winter, vort_spring, vort_summer, vort_fall};
div = {div_winter, div_spring, div_summer, div_fall};
q = {q_winter, q_spring, q_summer, q_fall};
temp = {temp_winter, temp_spring, temp_summer, temp_fall};
sp = {sp_winter, sp_spring, sp_summer, sp_fall};

F_atmos_to_sst = {NaN(length), NaN(length), NaN(length), NaN(length)};
F_sst_to_atmos = {NaN(length), NaN(length), NaN(length), NaN(length)};

mspe_full_atmos_to_sst = {NaN(length), NaN(length), NaN(length), NaN(length)};
mspe_full_sst_to_atmos = {NaN(length), NaN(length), NaN(length), NaN(length)};

mspe_reduced_atmos_to_sst = {NaN(length), NaN(length), NaN(length), NaN(length)};
mspe_reduced_sst_to_atmos = {NaN(length), NaN(length), NaN(length), NaN(length)};

sig_atmos_to_sst = {NaN(length), NaN(length), NaN(length), NaN(length)};
sig_sst_to_atmos = {NaN(length), NaN(length), NaN(length), NaN(length)};

times = importdata('data/times.mat');

'offset', offset

for i = 1:4
    sst_i = sst{i};
    vort_i = vort{i};
    div_i = div{i};
    q_i = q{i};
    temp_i = temp{i};
    sp_i = sp{i};

    F_atmos_to_sst_i = F_atmos_to_sst{i};
    F_sst_to_atmos_i = F_sst_to_atmos{i};
    mspe_full_atmos_to_sst_i = mspe_full_atmos_to_sst{i};
    mspe_full_sst_to_atmos_i = mspe_full_sst_to_atmos{i};
    mspe_reduced_atmos_to_sst_i = mspe_reduced_atmos_to_sst{i};
    mspe_reduced_sst_to_atmos_i = mspe_reduced_sst_to_atmos{i};
    sig_sst_to_atmos_i = sig_sst_to_atmos{i};
    sig_atmos_to_sst_i = sig_atmos_to_sst{i};

    for cell = 1:length
        if offset+cell-1 > 88838
            continue
        end
        sst_ts = sst_i(:, :, offset+cell-1)';
        vort_ts = vort_i(:, :, offset+cell-1)';
        div_ts = div_i(:, :, offset+cell-1)';
        q_ts = q_i(:, :, offset+cell-1)';
        temp_ts = temp_i(:, :, offset+cell-1)';
        sp_ts = sp_i(:, :, offset+cell-1)';

        if (sst_ts(1) == 9999) || (vort_ts(1) == 9999) || (div_ts(1) == 9999) || (sp_ts(1) == 9999) || (temp_ts(1) == 9999) || (q_ts(1) == 9999)
            continue
        end

        [nobs, ntrials] = size(sst_ts);
        X = NaN(6, nobs, ntrials);
        X(1, :, :) = sst_ts;
        X(2, :, :) = vort_ts;
        X(3, :, :) = div_ts;
        X(4, :, :) = q_ts;
        X(5, :, :) = temp_ts;
        X(6, :, :) = sp_ts;

        moAIC = times(offset+cell-1);

        [A, SIG] = tsdata_to_var(X, moAIC);
        assert(~isbad(A), 'VAR estimation failed');
        [G, info] = var_to_autocov(A, SIG);
        acerr = var_info(info, false);

        if acerr
            continue
        end

        [F_atmos_cause, mspe_full_atmos_cause, mspe_reduced_atmos_cause] = autocov_to_mvgc(G, [1], [2 3 4 5 6]);

        [F_sst_cause, mspe_full_sst_cause, mspe_reduced_sst_cause] = autocov_to_mvgc(G, [2 3 4 5 6], [1]);

        F_atmos_to_sst_i(cell) = F_atmos_cause;
        F_sst_to_atmos_i(cell) = F_sst_cause;

        mspe_full_atmos_to_sst_i(cell) = mspe_full_atmos_cause;
        mspe_full_sst_to_atmos_i(cell) = mspe_full_sst_cause;

        mspe_reduced_atmos_to_sst_i(cell) = mspe_reduced_atmos_cause;
        mspe_reduced_sst_to_atmos_i(cell) = mspe_reduced_sst_cause;

        pval_atmos_cause = mvgc_pval(F_atmos_cause, moAIC, nobs, ntrials, 1, 5, 0);
        pval_sst_cause = mvgc_pval(F_sst_cause, moAIC, nobs, ntrials, 5, 1, 0);

        sig_atmos_to_sst_i(cell) = pval_atmos_cause;
        sig_sst_to_atmos_i(cell) = pval_sst_cause;
    end
    F_atmos_to_sst{i} = F_atmos_to_sst_i;
    F_sst_to_atmos{i} = F_sst_to_atmos_i;
    mspe_full_atmos_to_sst{i} = mspe_full_atmos_to_sst_i;
    mspe_full_sst_to_atmos{i} = mspe_full_sst_to_atmos_i;
    mspe_reduced_atmos_to_sst{i} = mspe_reduced_atmos_to_sst_i;
    mspe_reduced_sst_to_atmos{i} = mspe_reduced_sst_to_atmos_i;
    sig_sst_to_atmos{i} = sig_sst_to_atmos_i;
    sig_atmos_to_sst{i} = sig_atmos_to_sst_i;
end

save(['data_atmos/F_seasonal_atmos_to_sst_' num2str(offset) '.mat'], 'F_atmos_to_sst');
save(['data_atmos/F_seasonal_sst_to_atmos_' num2str(offset) '.mat'], 'F_sst_to_atmos');
save(['data_atmos/mspe_full_seasonal_atmos_to_sst_' num2str(offset) '.mat'], 'mspe_full_atmos_to_sst');
save(['data_atmos/mspe_full_seasonal_sst_to_atmos_' num2str(offset) '.mat'], 'mspe_full_sst_to_atmos');
save(['data_atmos/mspe_reduced_seasonal_atmos_to_sst_' num2str(offset) '.mat'], 'mspe_reduced_atmos_to_sst');
save(['data_atmos/mspe_reduced_seasonal_sst_to_atmos_' num2str(offset) '.mat'], 'mspe_reduced_sst_to_atmos');
save(['data_atmos/sig_seasonal_atmos_to_sst_' num2str(offset) '.mat'], 'sig_atmos_to_sst');
save(['data_atmos/sig_seasonal_sst_to_atmos_' num2str(offset) '.mat'], 'sig_sst_to_atmos');
