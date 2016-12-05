Args <- commandArgs()
filename=Args[6]
pvalue_cutoff=0.1
read.table(filename,header=F)->mydata
rownum=nrow(mydata)
ks_fp_pvalues<-rep(NA,rownum)
ks_tp_pvalues<-rep(NA,rownum)
y<-rep(0,76)
for(i in c(60:76))
{
	y[i]=exp(-16.3505+0.2198*(i))
}
y[76]=y[76]+4
for(i in c(1:rownum))
{

	as.character(mydata[i,14][1][1])->base
	as.character(mydata[i,15][1][1])->site
	strsplit(base,split="")[[1]]->base
	as.numeric(strsplit(site,split=",")[[1]])->site
	snp<-site[which(base=="A"|base=="a"|base=="T"|base=="t"|base=="C"|base=="c"|base=="G"|base=="g")]
	snp<-sort(snp)
	##all the above code is to get the snp site in increase order
	fpbg<-y
	ts<-rep(0,(length(snp)+length(which(snp==76))))
	num76=length(which(snp==76))
	loop=1
	for(j in c(1:length(snp)))##let the snp site 
	{
		if(snp[j]<76)
		ts[j]=77-snp[j]
		if(snp[j]>76)
		ts[j]=snp[j]-75
		if(snp[j]==76)
		{
			ts[j]=1
			ts[length(ts)-num76+loop]=1
			loop=loop+1
		}
	}
	fpbg<-length(ts)*fpbg/sum(fpbg)
	bgs<-c(1:sum(round(fpbg)))
	kbg=1
	for( j in c(1:length(fpbg)) )
	{
		bnum=round(fpbg[j])
		if(bnum>0)
		for(l in c(1:bnum))
		{
			bgs[kbg]=j
			kbg=kbg+1
		}
	}
	ks_fp_pvalues[i]<-ks.test(ts,bgs)$p.value
	len<-length(ts)
	minsite<-site[which.min(site)]
	maxsite<-site[which.max(site)]
	midsite<-(minsite+maxsite)/2
	start<-(midsite-minsite+2)/(2*len)-0.5
	tpbg<-seq(start,by=(midsite-minsite+2)/len,length.out=len)
	tpbg<-sort(tpbg)
	tpbg<-round(tpbg)
	tpbg[tpbg<minsite]=minsite
	tpbg[tpbg>midsite]=midsite
	tpbg<-round(tpbg)
	ks_tp_pvalues[i]<-ks.test(ts,tpbg)$p.value
	
}

result=ks_tp_pvalues/ks_fp_pvalues
which(result>1|ks_tp_pvalues>pvalue_cutoff)->final_seq
rep(0,rownum)->check
rep(0,rownum)->pvalues
check[final_seq]=1
for(i in c(1:rownum))
{
	pvalues[i]=round(0.000,3)
	if(result[i]!="NaN")
	if(ks_tp_pvalues[i]>=0.001&&result[i]<=1)
	pvalues[i]=round(1.000*ks_tp_pvalues[i],3)
}
pvalues[result>1]=round(1.000,3)
cbind(pvalues,check,mydata)->mydata_new
resultname=paste(filename,".lable",sep="")
write.table(mydata_new,file=resultname,col.names=F,row.names=F)

