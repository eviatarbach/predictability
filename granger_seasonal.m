function [] = granger_seasonal(offset, length)
cd /lustre/ebach/causality

cd mvgc_v1.0

startup

cd ..

sst_winter = importdata('data/sst_winter.mat');
vort_winter = importdata('data/vort_winter.mat');
sst_spring = importdata('data/sst_spring.mat');
vort_spring = importdata('data/vort_spring.mat');
sst_summer = importdata('data/sst_summer.mat');
vort_summer = importdata('data/vort_summer.mat');
sst_fall = importdata('data/sst_fall.mat');
vort_fall = importdata('data/vort_fall.mat');

sst = {sst_winter, sst_spring, sst_summer, sst_fall};
vort = {vort_winter, vort_spring, vort_summer, vort_fall};

F_vort_to_sst = {NaN(length), NaN(length), NaN(length), NaN(length)};
F_sst_to_vort = {NaN(length), NaN(length), NaN(length), NaN(length)};

sig_vort_to_sst = {NaN(length), NaN(length), NaN(length), NaN(length)};
sig_sst_to_vort = {NaN(length), NaN(length), NaN(length), NaN(length)};

times = importdata(['data/times_' num2str(offset) '.mat']);
times_i = times{1};

'offset', offset

for i = 1:4
    sst_i = sst{i};
    vort_i = vort{i};

    F_vort_to_sst_i = F_vort_to_sst{i};
    F_sst_to_vort_i = F_sst_to_vort{i};
    sig_sst_to_vort_i = sig_sst_to_vort{i};
    sig_vort_to_sst_i = sig_vort_to_sst{i};

    for cell = 1:length
        if offset+cell-1 > 88838
            continue
        end
        sst_ts = sst_i(:, :, offset+cell-1)';
        vort_ts = vort_i(:, :, offset+cell-1)';

        if (sst_ts(1) == 9999) || (vort_ts(1) == 9999)
            continue
        end

        [nobs, ntrials] = size(sst_ts);
        X = NaN(2, nobs, ntrials);
        X(1, :, :) = sst_ts;
        X(2, :, :) = vort_ts;

        moAIC = times_i(cell);

        [A, SIG] = tsdata_to_var(X, moAIC);
        assert(~isbad(A), 'VAR estimation failed');
        [G, info] = var_to_autocov(A, SIG);
        acerr = var_info(info, false);

        if acerr
            continue
        end

        F = autocov_to_pwcgc(G);

        F_1 = F(1, 2);  % vort -> sst
        F_2 = F(2, 1);  % sst -> vort

        F_vort_to_sst_i(cell) = F_1;
        F_sst_to_vort_i(cell) = F_2;

        pval_vort_cause = mvgc_pval(F_1, moAIC, nobs*ntrials, 1, 1, 1, 0);
        pval_sst_cause = mvgc_pval(F_2, moAIC, nobs*ntrials, 1, 1, 1, 0);

        sig_vort_to_sst_i(cell) = pval_vort_cause;
        sig_sst_to_vort_i(cell) = pval_sst_cause;
    end
    F_vort_to_sst{i} = F_vort_to_sst_i;
    F_sst_to_vort{i} = F_sst_to_vort_i;
    sig_sst_to_vort{i} = sig_sst_to_vort_i;
    sig_vort_to_sst{i} = sig_vort_to_sst_i;
end

save(['data/F_seasonal_vort_to_sst_' num2str(offset) '.mat'], 'F_vort_to_sst');
save(['data/F_seasonal_sst_to_vort_' num2str(offset) '.mat'], 'F_sst_to_vort');
save(['data/sig_seasonal_vort_to_sst_' num2str(offset) '.mat'], 'sig_vort_to_sst');
save(['data/sig_seasonal_sst_to_vort_' num2str(offset) '.mat'], 'sig_sst_to_vort');
