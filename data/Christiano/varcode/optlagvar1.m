function optlagvar=optlagvar(data,hasconstant,maxlag);
% OPTLAGVAR displays different criteria for the choice of the optimal lag length for a VAR
%
% optlagvar(DATA,M0,MAXLAG)
%
% DATA:		each column is a variable
% M0 =		1 if a constant is included
% MAXLAG:	is the maximum order of the VAR
% Written by Riccardo Di Cecio.

capt=length(data)-maxlag-1;
ms_d=size(data,2);
for i=1:maxlag;
   [beta,sigma,u]=estimatevar(data,i,hasconstant);
   eval(['sig' num2str(i) '= sigma;']);
   ln_d(i)=log(det(eval(['sig' num2str(i)])));
   c_t(1)=2/(capt);
   c_t(2)=2*log(log(capt))/(capt);
   c_t(3)=(log(capt))/(capt);
   phi(i)=i*(ms_d^2);
end;
n_lag=(1:1:maxlag');
rstr=num2str(n_lag);
for i=1:3;
   c_r(:,i)=ln_d'+c_t(i)*phi';
end;
dgf=ms_d^2;
for i=maxlag:-1:2;
   LR(i)=capt*(ln_d(i-1)-ln_d(i));
   LRC(i)=(capt-1-ms_d*i)*(ln_d(i-1)-ln_d(i));
end;
LR=LR(2:end);
LRC=LRC(2:end);
for i=1:length(LR)
   chi_2(i)=1-chi2cdf(LR(i),dgf);
	chi_2c(i)=1-chi2cdf(LRC(i),dgf);   
end;
rstr1=2:1:maxlag;
rstr1=num2str(rstr1);
dgff=ones(maxlag-1,1)*dgf;
lrmat=[LR; chi_2; LRC;chi_2c;dgff'];
for i=1:maxlag-1;
   LR1(i)=capt*(ln_d(i)-ln_d(maxlag));
   LRC1(i)=(capt-1-ms_d*maxlag)*(ln_d(i)-ln_d(maxlag));
   dgf1(i)=(maxlag-i)*dgf;
end;
for i=1:length(LR1)
   chi_21(i)=1-chi2cdf(LR1(i),dgf1(i));
   chi_21c(i)=1-chi2cdf(LRC1(i),dgf1(i));
end;
rstr2=rstr(1:end-1);
lrmat1=[LR1;chi_21;LRC1;chi_21c;dgf1];
printmat(c_r,'Information criteria',rstr,'AIC HQ SC');
printmat(lrmat','LR tests',rstr1,'LRk/(k-1) p-value corLRk/(k-1) p-value dgf');
printmat(lrmat1','LR tests',rstr2,'LRkmax/k p-value corLRkmax/k p-value dgf');