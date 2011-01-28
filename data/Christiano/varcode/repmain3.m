load datrep;
data=datrep;
%('Y','P','PCOM','FF','TOTR','NBR','M');

nlags		= 4;
hasconst = 1;

[betaz,sigma,residuals]=estimatevar(data,nlags,hasconst);

a0rec=inv(chol(sigma)');

errshk = 3;
nstep = 15;
 
 impzmat1 =mkimprep(betaz,a0rec,nlags,errshk,nstep);
 
%To change the ordering is simple. Here we switch around y and p. 
%[betaz,sigma,residuals]=estimatevar(data(:,[2 1 3:end]),nlags,hasconst);
%a0rec=inv(chol(sigma)');
%errshk = 3;
%nstep = 15;
%impzmat2=mkimprep(betaz,a0rec,nlags,errshk,nstep);

 
pctg=0.05; %the confidence interval will be (1-2*pctg)*100 wide.
ndraws=100; % the number of draws. 
nobs=120; % the number of observations per draw. 
 % these are the bootstrap confidence intervals. 
[cilb,ciub,cilvarb,ciuvarb,varb]=mkimpci(betaz,a0rec,nlags,errshk,nstep,ndraws,nobs,pctg,residuals);

% this does the variance decomposition. 
[mse,msedecomp]=vardecomp(betaz,a0rec,nlags,errshk,nstep);

vdpctg = msedecomp./mse;

% this calculates the confidence intervals for the variance decomposition.
[cilvd,ciuvd,cilvarbv,ciuvarbv,varbv]=mkvdci(betaz,a0rec,nlags,errshk,nstep,ndraws,nobs,pctg,residuals);
kstep = [2 4 8 12];

for zk=1:4;
   lowci(:,zk) = diag(squeeze(cilvd(:,:,kstep(zk))));
   lowvar(:,zk) = diag(squeeze(cilvarbv(:,:,kstep(zk))));

table(:,zk) = diag(squeeze(vdpctg(:,:,kstep(zk))));
highci(:,zk) = diag(squeeze(ciuvd(:,:,kstep(zk))));
highvar(:,zk) = diag(squeeze(ciuvarbv(:,:,kstep(zk))));

end;

disp('The first column of Table 3 from CEE')
[lowci(:,1) table(:,1) highci(:,1)]

varnames=char('Y','P','PCOM','FF','TOTR','NBR','M');

figure
for zr=1:4;
   subplot(2,2,zr)
   hold on
   plot([ 100*cilvarb(:,zr) 100*ciuvarb(:,zr) 100*impzmat(:,zr)]);
   plot([100*cilb(:,zr) 100*ciub(:,zr)],'x-');
   title(['Response by ' varnames(zr,:)]); 
   axis tight
end;

suptitle('Shock to Fed Funds')
figure
for zr=5:7;
   subplot(2,2,zr-4)
   hold on
      plot([ 100*cilvarb(:,zr) 100*ciuvarb(:,zr) 100*impzmat(:,zr)]);
   plot([100*cilb(:,zr) 100*ciub(:,zr)],'x-');
   title(['Response by ' varnames(zr,:)]); 
   axis tight
end;
suptitle('Shock to Fed Funds')

%gtext(char('Solid Lines 2 Standard Deviation Based', 'x- Percentile Based '));


