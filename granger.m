nobs = {14245, 2849, 949};
max_order = {60, 12, 4};
mhtc = 'FDR';
alpha = 0.05;

cd mvgc_v1.0

startup

cd ..

sst_names = {'data/sst_01.mat', 'data/sst_05.mat', 'data/sst_15.mat'};
vort_names = {'data/vort_01.mat', 'data/vort_05.mat', 'data/vort_15.mat'};

F_vort_to_sst = {NaN(88838), NaN(88838), NaN(88838)};
F_sst_to_vort = {NaN(88838), NaN(88838), NaN(88838)};
sig_95_vort_to_sst = {NaN(88838), NaN(88838), NaN(88838)};
sig_95_sst_to_vort = {NaN(88838), NaN(88838), NaN(88838)};
times = {NaN(88838), NaN(88838), NaN(88838)};

parpool(20)

for i = 1:3
    nobs_i = nobs{i};
    sst_i = importdata(sst_names{i});
    vort_i = importdata(vort_names{i});

    times_i = times{i}; % Need to do this for parfor to work
    F_vort_to_sst_i = F_vort_to_sst{i};
    F_sst_to_vort_i = F_sst_to_vort{i};
    sig_95_sst_to_vort_i = sig_95_sst_to_vort{i};
    sig_95_vort_to_sst_i = sig_95_vort_to_sst{i};

    parfor cell = 1:88838
        sst_ts = reshape(sst_i(:, cell), [1, nobs_i]);
        vort_ts = reshape(vort_i(:, cell), [1, nobs_i]);

        if (sst_ts(1) == 9999) | (vort_ts(1) == 9999)
            continue
        end

        sst_ts = detrend(sst_ts);  % remove global warming signal

        X = [sst_ts; vort_ts];

        [~, ~, moAIC, ~] = tsdata_to_infocrit(X, max_order{i});
        times_i(cell) = moAIC;

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

        pval = mvgc_pval(F, moAIC, nobs_i, 1, 1, 1, 0);

        sig = significance(pval, alpha, mhtc);
        sig_1 = sig(1, 2);
        sig_2 = sig(2, 1);
        sig_95_vort_to_sst_i(cell) = sig_1;
        sig_95_sst_to_vort_i(cell) = sig_2;
    end
    times{i} = times_i;
    F_vort_to_sst{i} = F_vort_to_sst_i;
    F_sst_to_vort{i} = F_sst_to_vort_i;
    sig_95_sst_to_vort{i} = sig_95_sst_to_vort_i;
    sig_95_vort_to_sst{i} = sig_95_vort_to_sst_i;
end

save('data/F_vort_to_sst.mat', 'F_vort_to_sst');
save('data/F_sst_to_vort.mat', 'F_sst_to_vort');
save('data/sig_95_vort_to_sst.mat', 'sig_95_vort_to_sst');
save('data/sig_95_sst_to_vort.mat', 'sig_95_sst_to_vort');
save('data/times.mat', 'times');
