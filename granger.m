nobs = 13870;
morder = 10;
alpha = 0.1;
mhtc = 'FDR';

cd mvgc_v1.0

startup

cd ..

sst_01 = importdata('data/sst_01.mat');
vort_01 = importdata('data/vort_01.mat');
sst_05 = importdata('data/sst_05.mat');
vort_05 = importdata('data/vort_05.mat');
sst_15 = importdata('data/sst_15.mat');
vort_15 = importdata('data/vort_15.mat');

ratios = zeros(480, 241);
F_vort_to_sst = zeros(480, 241);
F_sst_to_vort = zeros(480, 241);
sig_90_vort_to_sst = zeros(480, 241);
sig_90_sst_to_vort = zeros(480, 241);
sig_95_vort_to_sst = zeros(480, 241);
sig_95_sst_to_vort = zeros(480, 241);

parpool(20)

parfor lat = 1:241
    for lon = 1:480
	lat, lon
        sst = sst_01(lon, lat, :);
        vort = vort_01(lon, lat, :);
        if sst(1) == -9999
            ratios(lon, lat) = 0;
	    continue
        end

	X = [reshape(sst, [1, nobs]); reshape(vort, [1, nobs])];
	[A, SIG] = tsdata_to_var(X, morder, 'OLS');
	assert(~isbad(A),'VAR estimation failed');
	[G,info] = var_to_autocov(A,SIG);
	var_info(info,true);

	F = autocov_to_pwcgc(G);
	pval = mvgc_pval(F,morder,nobs,1,1,1,0);
	sig = significance(pval,alpha,mhtc);

	F_1 = F(1, 2);  % vort -> sst
	F_2 = F(2, 1);  % sst -> vort

	if sig(1, 2) | sig(2, 1)
	    if F_1 > F_2
	        ratios(lon, lat) = F_1/F_2;
            else
		ratios(lon, lat) = -F_2/F_1;
	    end
	else
	    ratios(lon, lat) = 0;
	end
    end
end

save('granger.mat', 'ratios')
