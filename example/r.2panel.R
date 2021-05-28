### see note.txt
indset=(1:10); 
nind=length(indset); 
afn=c('./output/test.ps21.txt'); 
print(afn); 

######
	xx=scan("admix-1cm.truth");
	dim(xx)=c(length(xx)/20,20); 
	xx=t(xx); 
	xs=read.table("admix-1cm.marker",1)
	ys=read.table("./output/test.snpinfo.txt",1);
	ind=c(1:dim(xs)[1])[xs$rs %in% ys$rs]
	xxi=as.matrix(xx[,ind]); 
	zz=xxi[1:19,]+xxi[2:20,];
	zz=zz[c(1:10)*2-1,]
	zz=as.matrix(zz); 
	ns=dim(zz)[2]; 
###### zz is the true ancestry. 
    
	myfun=function(vx) {
		dim(vx)=c(2,length(vx)/2); 
		vx=as.matrix(t(vx)); 
		mu=vx[,1]; 
	}
	sy=matrix(0, nrow=nind, ncol=ns); 
	nt = 0; 
	sumy = 0;
	for ( fn in afn) 
	{
		if(!file.exists(fn)) {next;}
		print(fn); 
		nt = nt+1; 
		yy=scan(fn,skip=0); 
		sumy = sumy+yy;
	}
	if(nt==0){print("nt == 0");}
	{
		yy=sumy/nt;
		dim(yy)=c(length(yy)/nind,nind);  
   		yy=as.matrix(t(yy)); 
		sy=t(apply(yy,1,myfun)); 

	}
	fn="";
####### sy is elai inferred ancestry. 


   dist=matrix(0, nrow=10, ncol=2); 
   pset=c(1:ns); 
   par(mfrow=c(5,2)); 
   par(mar=c(1,1,1,1));
   for(k in c(1:10))
{
	   mu=sy;
	   plot(zz[k,pset],ylim=c(0,2),type="n",lwd=2.5,cex.lab=1.8, cex.axis=1.5, xlab="position along the chromosome",ylab="reference allele count");
	   points(zz[k,pset],ylim=c(0,2),type="l",lwd=2,cex.lab=1.8, cex.axis=1.5); 
	   points(sy[k,pset],ylim=c(0,2), type="l",col=2,lwd=2.2);                     
	   d1=mean(abs(mu[k,]-zz[k,]));
	   d3=cor(mu[k,], zz[k,]);
	   dist[k,]=c(d1,d3);
	   print(c(d1,d3)); 
}
###### plot the truth as the inferred. 

