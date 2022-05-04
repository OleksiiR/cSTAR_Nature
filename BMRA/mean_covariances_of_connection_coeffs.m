data=dlmread('C:\Users\TAPESH\Dropbox\codes\BVSA_MRA\Netinf\results\lm1215_tgf_cc_all.txt');
C=cov(data);
M=mean(data);
params.mean=M;
params.cov=C;
save('C:\Users\TAPESH\Dropbox\codes\BVSA_MRA\Netinf\results\lim1215_tgf','params');

data=dlmread('C:\Users\TAPESH\Dropbox\codes\BVSA_MRA\Netinf\results\lm1215_igf_cc_all.txt');
C=cov(data);
M=mean(data);
params.mean=M;
params.cov=C;
save('C:\Users\TAPESH\Dropbox\codes\BVSA_MRA\Netinf\results\lim1215_igf','params');

data=dlmread('C:\Users\TAPESH\Dropbox\codes\BVSA_MRA\Netinf\results\hct116_tgf_cc_all.txt');
C=cov(data);
M=mean(data);
params.mean=M;
params.cov=C;
save('C:\Users\TAPESH\Dropbox\codes\BVSA_MRA\Netinf\results\hct116_tgf','params');

data=dlmread('C:\Users\TAPESH\Dropbox\codes\BVSA_MRA\Netinf\results\hct116_igf_cc_all.txt');
C=cov(data);
M=mean(data);
params.mean=M;
params.cov=C;
save('C:\Users\TAPESH\Dropbox\codes\BVSA_MRA\Netinf\results\hct116_igf','params');

data=dlmread('C:\Users\TAPESH\Dropbox\codes\BVSA_MRA\Netinf\results\ht29_tgf_cc_all.txt');
C=cov(data);
M=mean(data);
params.mean=M;
params.cov=C;
save('C:\Users\TAPESH\Dropbox\codes\BVSA_MRA\Netinf\results\ht29_tgf','params');

data=dlmread('C:\Users\TAPESH\Dropbox\codes\BVSA_MRA\Netinf\results\ht29_igf_cc_all.txt');
C=cov(data);
M=mean(data);
params.mean=M;
params.cov=C;
save('C:\Users\TAPESH\Dropbox\codes\BVSA_MRA\Netinf\results\ht29_igf','params');

data=dlmread('C:\Users\TAPESH\Dropbox\codes\BVSA_MRA\Netinf\results\sw403_tgf_cc_all.txt');
C=cov(data);
M=mean(data);
params.mean=M;
params.cov=C;
save('C:\Users\TAPESH\Dropbox\codes\BVSA_MRA\Netinf\results\sw403_tgf','params');

data=dlmread('C:\Users\TAPESH\Dropbox\codes\BVSA_MRA\Netinf\results\sw403_igf_cc_all.txt');
C=cov(data);
M=mean(data);
params.mean=M;
params.cov=C;
save('C:\Users\TAPESH\Dropbox\codes\BVSA_MRA\Netinf\results\sw403_igf','params');

data=dlmread('C:\Users\TAPESH\Dropbox\codes\BVSA_MRA\Netinf\results\sw480_tgf_cc_all.txt');
C=cov(data);
M=mean(data);
params.mean=M;
params.cov=C;
save('C:\Users\TAPESH\Dropbox\codes\BVSA_MRA\Netinf\results\sw480_tgf','params');

data=dlmread('C:\Users\TAPESH\Dropbox\codes\BVSA_MRA\Netinf\results\sw480_igf_cc_all.txt');
C=cov(data);
M=mean(data);
params.mean=M;
params.cov=C;
save('C:\Users\TAPESH\Dropbox\codes\BVSA_MRA\Netinf\results\sw480_igf','params');