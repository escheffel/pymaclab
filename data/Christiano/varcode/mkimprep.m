function impz=mkimprep(betaz,a0,nlags,errshk,nstep);
% MKIMPREP Calculate the first nstep impulse response to a shock to the errshk error term
% to the system (1-betaz(L))y(t) = inv(a0)*error(t);

betaz = betaz';
a0inv = inv(a0);

nvars = size(a0,2);

if nvars ~= size(betaz,2);
   error('betaz and a0 are not the correct dimensions')
end;

hascon = size(betaz,1)-nvars*nlags;
if ~or(hascon==1,hascon==0);
   error('betaz is the wrong size')
end;

shock = [a0inv(:,errshk)';zeros(nstep-1,nvars)]

impz=zeros(nstep,nvars);  
zlag=zeros(1,nlags*nvars);

for jj=1:nstep;
   impz(jj,:)=shock(jj,:)+zlag*betaz(1+hascon:length(betaz),:);
   zlag=[impz(jj,:) zlag(1,1:(nlags-1)*nvars)];
end;
