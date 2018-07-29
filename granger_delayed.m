function [] = granger_delayed(offset)
cd /lustre/ebach/causality

nobs = {13870, 2774, 924};
max_delay = {30, 15, 10};
delays = {0:max_delay{1}, 0:max_delay{2}, 0:max_delay{3}};

cd mvgc_v1.0

startup

cd ..

sst_01 = importdata('data/sst_01.mat');
vort_01 = importdata('data/vort_01.mat');
sst_05 = importdata('data/sst_05.mat');
vort_05 = importdata('data/vort_05.mat');
sst_15 = importdata('data/sst_15.mat');
vort_15 = importdata('data/vort_15.mat');

sst = {sst_01, sst_05, sst_15};
vort = {vort_01, vort_05, vort_15};

F_vort_to_sst = {NaN(10, 241, max_delay{1} + 1), NaN(10, 241, max_delay{2} + 1), NaN(10, 241, max_delay{3} + 1)};
F_sst_to_vort = {NaN(10, 241, max_delay{1} + 1), NaN(10, 241, max_delay{2} + 1), NaN(10, 241, max_delay{3} + 1)};

mspe_vort_to_sst = {NaN(10, 241, max_delay{1} + 1), NaN(10, 241, max_delay{2} + 1), NaN(10, 241, max_delay{3} + 1)};
mspe_sst_to_vort = {NaN(10, 241, max_delay{1} + 1), NaN(10, 241, max_delay{2} + 1), NaN(10, 241, max_delay{3} + 1)};

sig_vort_to_sst = {NaN(10, 241, max_delay{1} + 1), NaN(10, 241, max_delay{2} + 1), NaN(10, 241, max_delay{3} + 1)};
sig_sst_to_vort = {NaN(10, 241, max_delay{1} + 1), NaN(10, 241, max_delay{2} + 1), NaN(10, 241, max_delay{3} + 1)};

times = importdata('data/times.mat');

%parpool(10)

'offset', offset

for i = 1:3
    nobs_i = nobs{i};
    sst_i = sst{i};
    vort_i = vort{i};
    delays_i = delays{i};

    times_i = times{i}; % Need to do this for parfor to work
    F_vort_to_sst_i = F_vort_to_sst{i};
    F_sst_to_vort_i = F_sst_to_vort{i};
    mspe_vort_to_sst_i = mspe_vort_to_sst{i};
    mspe_sst_to_vort_i = mspe_sst_to_vort{i};
    sig_sst_to_vort_i = sig_sst_to_vort{i};
    sig_vort_to_sst_i = sig_vort_to_sst{i};

    j = 1;

    for delay = delays_i
        F_vort_to_sst_ij = F_vort_to_sst_i(:, :, j);
        F_sst_to_vort_ij = F_sst_to_vort_i(:, :, j);
        mspe_vort_to_sst_ij = mspe_vort_to_sst_i(:, :, j);
        mspe_sst_to_vort_ij = mspe_sst_to_vort_i(:, :, j);
        sig_vort_to_sst_ij = sig_vort_to_sst_i(:, :, j);
        sig_sst_to_vort_ij = sig_sst_to_vort_i(:, :, j);

        for ilon = 1:10
            lon = str2num(offset) + (ilon - 1)
            for lat = 1:241
                lat, lon
                sst_ts = reshape(sst_i(lon, lat, :), [1, nobs_i]);
                vort_ts = reshape(vort_i(lon, lat, :), [1, nobs_i]);

                if (sst_ts(1) == -9999) || (vort_ts(1) == -9999)
                    continue
                end

                sst_ts = detrend(sst_ts);  % remove global warming signal

                moAIC = times_i(lon, lat);

                sst_ts_effect = sst_ts(1+delay:end);
                sst_ts_cause = sst_ts(1:end-delay);

                vort_ts_effect = vort_ts(1+delay:end);
                vort_ts_cause = vort_ts(1:end-delay);

                X_sst_cause = [sst_ts_cause; vort_ts_effect];
                X_vort_cause = [sst_ts_effect; vort_ts_cause];

                [A_sst_cause, SIG_sst_cause] = tsdata_to_var(X_sst_cause, moAIC);
                assert(~isbad(A_sst_cause), 'VAR estimation failed');
                [G_sst_cause, info_sst_cause] = var_to_autocov(A_sst_cause, SIG_sst_cause);
                acerr_sst_cause = var_info(info_sst_cause, false);

                [A_vort_cause, SIG_vort_cause] = tsdata_to_var(X_vort_cause, moAIC);
                assert(~isbad(A_vort_cause), 'VAR estimation failed');
                [G_vort_cause, info_vort_cause] = var_to_autocov(A_vort_cause, SIG_vort_cause);
                acerr_vort_cause = var_info(info_vort_cause, false);

                if acerr_sst_cause || acerr_vort_cause
                    continue
                end

                [F_1, mspe_1] = autocov_to_pwcgc(G_vort_cause);  % vort -> sst
                F_1 = F_1(1, 2)
                mspe_1 = mspe_1(1, 2)

                [F_2, mspe_2] = autocov_to_pwcgc(G_sst_cause);  % sst -> vort
                F_2 = F_2(2, 1)
                mspe_2 = mspe_2(2, 1)

                F_vort_to_sst_ij(ilon, lat) = F_1;
                F_sst_to_vort_ij(ilon, lat) = F_2;

                mspe_vort_to_sst_ij(ilon, lat) = mspe_1;
                mspe_sst_to_vort_ij(ilon, lat) = mspe_2;

                pval_vort_cause = mvgc_pval(F_1, moAIC, nobs_i - delay, 1, 1, 1, 0);
                pval_sst_cause = mvgc_pval(F_2, moAIC, nobs_i - delay, 1, 1, 1, 0);

                sig_vort_to_sst_ij(ilon, lat) = pval_vort_cause;
                sig_sst_to_vort_ij(ilon, lat) = pval_sst_cause;
            end
        end
        F_vort_to_sst_i(:, :, j) = F_vort_to_sst_ij;
        F_sst_to_vort_i(:, :, j) = F_sst_to_vort_ij;
        mspe_vort_to_sst_i(:, :, j) = mspe_vort_to_sst_ij;
        mspe_sst_to_vort_i(:, :, j) = mspe_sst_to_vort_ij;
        sig_vort_to_sst_i(:, :, j) = sig_vort_to_sst_ij;
        sig_sst_to_vort_i(:, :, j) = sig_sst_to_vort_ij;

        j = j + 1;
    end
    F_vort_to_sst{i} = F_vort_to_sst_i;
    F_sst_to_vort{i} = F_sst_to_vort_i;
    mspe_vort_to_sst{i} = mspe_vort_to_sst_i;
    mspe_sst_to_vort{i} = mspe_sst_to_vort_i;
    sig_sst_to_vort{i} = sig_sst_to_vort_i;
    sig_vort_to_sst{i} = sig_vort_to_sst_i;
end

save(['data/F_vort_to_sst_' offset '.mat'], 'F_vort_to_sst');
save(['data/F_sst_to_vort_' offset '.mat'], 'F_sst_to_vort');
save(['data/mspe_vort_to_sst_' offset '.mat'], 'mspe_vort_to_sst');
save(['data/mspe_sst_to_vort_' offset '.mat'], 'mspe_sst_to_vort');
save(['data/sig_vort_to_sst_' offset '.mat'], 'sig_vort_to_sst');
save(['data/sig_sst_to_vort_' offset '.mat'], 'sig_sst_to_vort');
