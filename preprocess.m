nobs = 13870;

fid=fopen('data/sstah_eraint_0.75x0.75_daily365_1979-2016.dat','r');
data=fread(fid,inf,'real*4');
fclose(fid);
sst_daily=reshape(data,[480,241,nobs]);

fid=fopen('data/vort850ahocn_eraint_0.75x0.75_daily365_1979-2016.dat','r');
data=fread(fid,inf,'real*4');
fclose(fid);
vort_daily=reshape(data,[480,241,nobs]);

% Need to take out 10 days to be divisible by 15
sst_15 = mean(reshape(sst_daily(:, :, 11:end), [480, 241, 13860/15, 15]), 4);
vort_15 = mean(reshape(vort_daily(:, :, 11:end), [480, 241, 13860/15, 15]), 4);

save('data/sst_15.mat', 'sst_15', '-v7.3');
save('data/vort_15.mat', 'vort_15', '-v7.3');

fid=fopen('data/sstah_eraint_0.75x0.75_daily365_1979-2016.dat','r'); 
data=fread(fid,inf,'real*4');
fclose(fid);
sst_daily=reshape(data,[480,241,nobs]);

fid=fopen('data/vort850ahocn_eraint_0.75x0.75_daily365_1979-2016.dat','r'); 
data=fread(fid,inf,'real*4');
fclose(fid);
vort_daily=reshape(data,[480,241,nobs]);

save('data/sst_01.mat', 'sst_daily', '-v7.3');
save('data/vort_01.mat', 'vort_daily', '-v7.3');

fid=fopen('data/sstah_eraint_0.75x0.75_pentad_1979-2016.dat','r'); 
data=fread(fid,inf,'real*4');    
fclose(fid);
sst_pentad=reshape(data,[480,241,nobs/5]);

fid=fopen('data/vort850ahocn_eraint_0.75x0.75_pentad_1979-2016.dat','r'); 
data=fread(fid,inf,'real*4');    
fclose(fid);
vort_pentad=reshape(data,[480,241,nobs/5]); 

save('data/sst_05.mat', 'sst_pentad', '-v7.3');
save('data/vort_05.mat', 'vort_pentad', '-v7.3');
