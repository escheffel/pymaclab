load c:\robswork\identeich\szdat_dm.dat;
data=szdat_dm;
data(:,1)=data(:,1)./sqrt(10.0); 

nl 		= 4;

[bz,sig,rz]=estimatevar(data,nl);

a0=estnonreca(sig,20,'mksimzha','svsz')

eshk=3;
nstp = 15;
 
 
 impz=mkimprep(bz,a0,nl,eshk,nstp);

pc=0.05;
nd=10; %should be more like a 100
nb=120;


[lp,up,lv,uv]=mkimpci2(bz,a0,nl,eshk,nstp,nd,nb,pc,rz,sig,3,'mksimzha','svsz');


[mse,msed]=vardecomp(bz,a0,nl,eshk,nstp);

vdpctg = msed./mse;

[vlp,vup,vlv,vuv]=mkvdci(bz,a0,nl,eshk,nstp,nd,nb,pc,rz,sig,3,'mksimzha','svsz');

kstep = [2 4 8 12];

for zk=1:4;
   lowci(:,zk) = diag(squeeze(vlp(:,:,kstep(zk))));
   lowvar(:,zk) = diag(squeeze(vlv(:,:,kstep(zk))));

table(:,zk) = diag(squeeze(vdpctg(:,:,kstep(zk))));
highci(:,zk) = diag(squeeze(vup(:,:,kstep(zk))));
highvar(:,zk) = diag(squeeze(vuv(:,:,kstep(zk))));

end;

 
 
varnames=char('Y','P','PCOM','FF','TOTR','NBR','M');
for zr=1:7;
   subplot(4,2,zr)
   hold on
   plot([impz(:,zr) cilb(:,zr) ciub(:,zr)],'.');
   % title(['Response by ' varnames(zr,:)]); 
   axis tight
end;

