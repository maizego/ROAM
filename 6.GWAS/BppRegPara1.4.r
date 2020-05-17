# FORWARD AND BACKWARD STEPWISE REGRESSION WITH BOOTSTRAP RESAMPLING STRATEGY CARRY OUT THE GENOWIDE ASSOCIATION ANALYSIS
# ON THE MULTIPLE PARENTAL RECOMBINATION POPULATIONS IN MAIZE. THIS GWAS SCRIPT FIRST HOLD THE 'POP' TERM IN PRIOR, AND THEM 
# PERFORM THE STEPWISE REGRESSION ALONG THE 'GENOME' SNP BY SNP WITH THE DEGREE OF FREE OF ONE FOR EACH SNP. 
# THE MUTLTIPLE-VARIABLE REGRESSION MODEL IS: y = Xb + Sa + e

# INPUT DATA INCLUDE FOLLOWING FOUR FILES:
# gen --- A GENOTYPE MATRIX WHITH 'N' ROWS AND 'M' COLUMN, REPRESENTING INDIVIDUALS AND MARKERS(SNP), RESPECTIVELY.
# phe --- A DATAFRAME WITH THREE COLUMNS: 'line','pop','trait.name'.
# inf --- A DATAFRAME OF GENOTYPE INFOMATION WITH THREE COLUMNS: 'snp','alleles','chr','pos'.
# qtl --- A DATAFRAME WITH AT LEAST FOWLLOWING COLUMNS: 'chr','pos.p','pos.l','pos.r',prior QTL information for controlling background effect.

# FIVE PARAMETERS OF GWAS SCANING AS FOLLOWING:
# cores --- CORE NUMBERS FOR PARALLEL CALCULATING PROCEDURES.
# nboot --- NUMBERS OF BOOTSTRAP SAMPLING.
# nchunk ---- NUMBERS OF CHUNKS SPLITED FROM THE TOTAL DATASET FOR ENHANCING THE CALCULATING SPEED OF REGRESSION IF THE COMPUTER HAS NOT ENOUGH MEMORY.
# nperm_reg --- PERMUTATION TIMES TO DETERMINE THE P VALUE CUTOFF WHEN THE SNP ENTERRING AND LEAVING THE REGRESSION MODEL.
# nperm_boot --- PERMUTATION TIMES TO DETERMINE THE BPP CUTOFF. INCREASING THE VALUE WILL SIGNIFICANTLY INCREASE THE OPERATION PERIOD.

# OUTPUT RESULTS INCLUDES FOLLOWING NINE FILES:
# Cutoffs: the cutoff values of chromosome based scanning, bootstrap and final backward procedures, respectively.
# SingleReg: the results of single stepwise regression analysis, "#NA"s indicate no available significance in specific chromosomes.
# Permbpp: the permutation results of determining bootstrap cutoff.
# SigSNPbpp: the raw results of bootstrap based stepwise regression analysis.
# SigSNPfilt: the filter results of bootstrap based stepwise regression analysis.
# SNP.finalmodel: the results of final backward regression analysis.
# Params.finalmodel: the significance parameters of the final backward regression.
# Fit.model: the fit of three models including only population term, only SNP term as well as both population and SNP terms, respectively.
# uniq.QTL: the prior QTLs detected by the whole GWAS procedure.

# REQUIRED R PACKAGES:
# parallel: optional for enhancing the calculation efficiency, following information shows the details.
# xlsx: output the results as Microsoft excel files.

# NOTES:
# THE SCRIPT IMPLEMENT THE PARALLEL ALGORITHM FROM THE R PACKAGE 'parallel' WHICH IS VERY SIMILAR TO 'multicores', 
# IF THE SCRIPT RUN ON THE WINDOWS SYSTEM, THE 'cores' SHOULD BE SET TO ONE, AND IF ON THE LINUX SYSTEM,
# THE 'cores' CAN SET TO MORE THAN ONE TO ACHIEVE THE PARALLELE. THE GENOTYPE SHOULD HAVE NO MISSING VALUE.
# UPDATAE V1.4:
# REVISE THE CODE STRUCTURE, AND ELIMINATE THE COFACTOR MATRIX FROM THE INPUT FILES, WHICH IS REPLACED WITH QTL DATAFRAME.
# THESE REVISION REDUCE THE COMPLAX LABOR DURING PREPARING THE INPUT FILES.

# ANY BUGS OR QUESTIONS WILL BE WELCOME AND PLEASE CONTACT WITH YINGJIE XIAO(Ph.D), EMAIL: shanren0179@gmail.com.

BppRegPara = function(gen,phe,inf,qtl,cores=1,nboot=100,nchunk=1,nperm_reg=500,nperm_boot=2){

if(Sys.info()[1] == 'Windows' & cores > 1) {
stop('---------OS \'WINDOWS\' CAN NOT IMPLEMENT MULTIPLE-CORES PARALLELIZATION BY \'PARALLEL\' PACKAGE--------\n')}
start_time = Sys.time()

gen = as.matrix(gen)
stopifnot(sum(is.na(gen)) == 0)
n0 = sum(gen == 0)
n1 = sum(gen == 1)
n2 = sum(gen == 2)
stopifnot(n0+n1+n2 == length(gen))

phe = na.omit(phe)
trait = names(phe)[3]
line.id = intersect(rownames(gen),phe$line)
phe = phe[phe$line %in% line.id,]
phe = phe[order(phe$pop),]
stopifnot(all(qtl$pop %in% phe$pop))
y = phe[,3]

pop = factor(seq(length(unique(phe$pop)))[match(phe$pop,unique(phe$pop))])
qtl$pop = factor(seq(length(unique(phe$pop)))[match(qtl$pop,unique(phe$pop))])
gen = gen[match(phe$line,rownames(gen)),]
stopifnot(identical(as.character(inf$snp),colnames(gen)))

# EXCLUDE NON-POLYMORPHIC SNP MARKERS IN THE WHOLE POPOLATION	
if(any(!c('chr','snp','pos','alleles') %in% colnames(inf))){
stop('------------\'inf\'SHOULD CONTAIN FOUR COLUMNS:\'chr\',\'pos\',\'snp\',\'alleles\'-------------\n')}
gen = gen[,apply(gen,2,function(x) any(x[1] != x))]
inf = inf[match(colnames(gen),inf$snp),]
inf = inf[order(inf$chr,inf$pos),]
gen = gen[,match(inf$snp,colnames(gen))]
stopifnot(length(y) == nrow(gen))
stopifnot(length(y) == length(pop))

if(!is.null(qtl) & any(!c('chr','pos.p','pos.l','pos.r') %in% colnames(qtl))){
stop('-----------QTL INFO MUST CONTAIN AT LEAST TWO COLUMNS,\'chr\',\'pos.p\',\'pos.l\',\'pos.r\'-------------\n')}
if(!is.null(qtl)){
qtl = qtl[order(qtl$chr,qtl$pos.l),]
maxp = c(0, cumsum(tapply(inf$pos, inf$chr, max) + sum(tapply(inf$pos, inf$chr, max)) / 200))	
mpos = inf$pos + maxp[inf$chr]
qtl.pos = qtl$pos.p + maxp[qtl$chr]
qtl.l = qtl$pos.l + maxp[qtl$chr]
qtl.r = qtl$pos.r + maxp[qtl$chr]

cat('-----------SELECT THE QTL PEAK MARKERS AS COVARIATES STRAT!-----------------------------------------------\n\n')
library(parallel)
mc1=getOption('mc.cores',ifelse(cores>20,10,cores/2))
pop_poly = do.call(cbind,mclapply(split(as.data.frame(gen,check.names=F),pop),
				function(gp)apply(gp,2,function(x)any(x[1] != x)),mc.cores=mc1))
if(ncol(pop_poly) != length(unique(pop))) {
stop('-----------THE COLUMN DOES NOT INDICATE POPULATIONS-------------------------------------------------------\n')}

Mk = vector()
for(j in 1:nrow(qtl)){
	dp = qtl.pos[j] - mpos[pop_poly[,colnames(pop_poly)%in%qtl$pop[j]]]
	Mk[j] = as.character(inf$snp[pop_poly[,colnames(pop_poly)%in%qtl$pop[j]]][which.min(abs(dp))])
	}

cof = as.matrix(gen[,match(Mk,colnames(gen))])
colnames(cof) = qtl$chr
qtl$snp = Mk
} else cof = NULL

# EXCLUDE THE COMPLETE POLYMORPHYSMS WITH COMPLETE LINKAGE DISEQUILIBRIUM
if(is.data.frame(gen)) gen=as.matrix(gen)
gen = unique(gen,MARGIN=2)	
inf.uni = inf[match(colnames(gen),inf$snp),]
chr = inf.uni$chr

cat('-----------ASSEMBLE THE GENOME BIN BASED ON LINKAGE DISEQUILIBRIUM STRAT!--------------------------------\n\n')
bin_inf = genome_bin(inf,inf.uni)

rm(pop_poly,dp,inf,inf.uni)
gc()

false_p_tol = length(unique(chr)) * nperm_boot		
			
mc = getOption('mc.cores',cores)

if(is.null(cof)) {
	x = data.frame(gen,check.names=F)
	rm(gen)
	cofc = as.matrix(rep(1,nrow(x)))

	cat('---------PERMUTATION TESTS FOR REGRESSION START!---------------------------------------------------------','\n')
	perm_plist = unlist(mclapply(1:nperm_reg,perm_reg,x,y,pop,cofc,nchunk,mc.cores = mc))
	p_cutoff = quantile(perm_plist,0.05,names = F,na.rm = T,type = 1)
	cutoff_reg = c('P_cutoff',p_cutoff)
	
	# @STEP 1: SINGLE STEPWISE REGRESSION SCANING IN THE FULL DATASET
	comm = 1
	res.reg = stepwise_reg(x,y,pop,cofc,comm,p_cutoff,nboot,nchunk)
	reg_single = data.frame(bin_inf[match(rownames(res.reg),bin_inf$Tag_snp),],res.reg)
	cat('---STEP 1: SINGLE STEPWISE SCANING FOR FULL DATA DONE!---------------','\n')
	
	# @STEP 2: BOOTSTRAP BASED STEPWISE REGRESSION
	cat('---STEP 2: BOOTSTRAP BASED STEPWISE SCANING START!------------------','\n')
	# PERMUTATION TESTS FOR BOOTSTRAP STARTS
	comm = 2
	res.perm = do.call('rbind',lapply(1:nperm_boot,perm_boot,x,y,pop,cofc,bin_inf,comm,p_cutoff,nboot,nperm_boot,nchunk,mc))
	perm_bpp = res.perm[!is.na(res.perm$Tag_snp),]
	
	# BOOTSTRAP START BY RESAMPLING THE EQUAL NUMBER OF INDIVIDUALS REPEATALLY, IN WHICH THE STEPWISE REGRESSION WAS PERFORMED
	cat('---------BOOTSTRAP REGRESSION ANALYSIS START!-----------------------------------------------------------------','\n')
	comm = 3
	res.sig = nbootf(x,y,pop,cofc,comm,p_cutoff,nboot,nchunk,mc)
	sig_bpp = boot_summary(res.sig,nboot,bin_inf)																																			
	} 
	else {
	xc = lapply(tapply(1:ncol(gen),chr,unique),function(m)gen[,m])
	reg_chr = list()
	perm_chr = list()
	sig_chr = list()
	cutoff_chr = list()
	rm(gen)
	
	for(c in 1:length(unique(chr))){
	x = data.frame(xc[[c]],check.names=F)
	cofc = cof[,colnames(cof) != c]

	cat('---------CHR',c,'.PERMUTATION TESTS FOR REGRESSION START!---------------------------------------------------------','\n')
	perm_plist = unlist(mclapply(1:nperm_reg,perm_reg,x,y,pop,cofc,nchunk,mc.cores = mc))
	p_cutoff = quantile(perm_plist,0.05,names = F,na.rm = T,type = 1)
	cutoff_chr[[c]] = c(paste('P_cutoff.chr',c,sep=''),p_cutoff)

	# @STEP 1: SINGLE STEPWISE REGRESSION SCANING IN THE FULL DATASET
	comm = 1
	res.reg = stepwise_reg(x,y,pop,cofc,comm,p_cutoff,nboot,nchunk)
	reg_chr[[c]] = data.frame(bin_inf[match(rownames(res.reg),bin_inf$Tag_snp),],res.reg)
	cat('---STEP.1:SINGLE STEPWISE SCANING IN FULL-DATA FOR CHR',c,'DONE!---------------','\n')
	
	# @STEP 2: BOOTSTRAP BASED STEPWISE REGRESSION
	cat('---STEP.2:BOOTSTRAP BASED STEPWISE SCANING FOR CHR',c,'START!------------------','\n')
	cat('---------CHR',c,'.PERMUTATION FOR BOOTSTRAP STRAT!-----------------------------------------------------------------','\n')
	comm = 2
	res.perm = do.call('rbind',lapply(1:nperm_boot,perm_boot,x,y,pop,cofc,bin_inf,comm,p_cutoff,nboot,nperm_boot,nchunk,mc))
	perm_chr[[c]] = res.perm[!is.na(res.perm$Tag_snp),]

	cat('---------CHR',c,'.BOOTSTRAP REGRESSION ANALYSIS START!-----------------------------------------------------------------','\n')
	comm = 3
	res.sig = nbootf(x,y,pop,cofc,comm,p_cutoff,nboot,nchunk,mc)
	sig_chr[[c]] = boot_summary(res.sig,nboot,bin_inf)
	
	cat('---------CHR',c,'.PERMUTATIONS BASED P-VALUE CUTOFF:',p_cutoff,'----------------------------------------------','\n')
	}
	
	reg_single = do.call(rbind,reg_chr)
	perm_bpp = do.call(rbind,perm_chr)
	sig_bpp = do.call(rbind,sig_chr)
	cutoff_reg = do.call(rbind,cutoff_chr)
	}
	
	bpp_cutoff = ifelse(nrow(perm_bpp) == 0,0.02,ifelse(with(perm_bpp,min(sort(BPP,decreasing=T)[1:false_p_tol],na.rm=T) < 0.02),
													0.02,min(sort(perm_bpp$BPP,decreasing=T)[1:false_p_tol],na.rm=T)))
	sig_filt = subset(sig_bpp,BPP >= bpp_cutoff)
	
# @STEP 3: FINAL MODEL SELECTION REGARDING FULL-DATA SINGLE REGRESSION AND BOOTSTRAP STEPWISE REGRESSION TOGETHER AS CANDIDATE SNPS
	cat('---STEP 3: FINAL MODEL SELECTION COMBINING FULL-DATA REGRESSION & BOOTSTRAP SCANING STRAT!----------','\n')
	combin.snp = c(as.character(reg_single$Tag_snp),as.character(sig_filt$Tag_snp))
	if(is.null(cof)) x.final = data.frame(x[,colnames(x) %in% combin.snp],check.names=F)
	else x.final = data.frame(do.call('cbind',xc)[,colnames(do.call('cbind',xc)) %in% combin.snp],check.names=F)
	
	bin.final = bin_inf[match(colnames(x.final),bin_inf$Tag_snp),]
	sig.pos = bin.final$Pos_sta + maxp[bin.final$Chr]
	
	# SUMMARY THE GWAS POINTS OVERRED BY QTL INTERVALS
	gwas.hit = sapply(1:nrow(bin.final), function(i) {
				hit = sig.pos[i] >= qtl.l & sig.pos[i] <= qtl.r
				hit.id = ifelse(any(hit),which(hit),0)
				return(hit.id)
				})
	
	# SPLIT FINAL SNPs INTO TWO GROUPS BASED ON WETHER LOCATED IN THE QTL REGIONS
	
	p.cut.final = ifelse(is.null(cof),p_cutoff,median(as.numeric(cutoff_reg[,2])))
	Cutoffs = rbind(cutoff_reg,c('bpp_cutoff',bpp_cutoff),c('p.cut.final',p.cut.final))
	
	# HIT SNP STEPWISE SELECTION PER QTL BY QTL
	hit.snp = NULL
	for(i in unique(gwas.hit[gwas.hit != 0])){
	x.hit = x.final[,colnames(x.final) %in% bin.final$Tag_snp[gwas.hit == i],drop=F]
	reg.hit = stepwise_reg(x.hit,y,pop,cof[,qtl$chr != qtl$chr[i]],comm=1,p.cut.final,nboot,nchunk=1,b=1)
	hit.snp = c(hit.snp,rownames(reg.hit)[!is.na(reg.hit[,1])])
	}
	
	# SUMMARY THE QTL INVERVALS COVERRED BY THE GWAS SNPs
	qtl.hit = sapply(1:nrow(qtl), function(i) any(qtl.l[i] <= sig.pos & qtl.r[i] >= sig.pos))
	
	qtl.name = paste(paste('QTL',qtl$chr,sep=''),unlist(tapply(qtl$chr,qtl$chr,function(x)seq(length(x))),use.names=F),sep='.')
	colnames(cof)=qtl.name
	qtl = cbind(qtl.num=qtl.name,qtl)
	# UNHIT SNP STEPWISE SELECTION BY FORCING THE POP, UNCOVERED QTLs, HITTED SNPs AS COVARIATE IN MODEL
	x.fix = as.matrix(x.final[,colnames(x.final) %in% hit.snp,drop=F])
	x.test = x.final[,colnames(x.final) %in% bin.final$Tag_snp[gwas.hit == 0],drop=F]
	reg.unhit = stepwise_reg(x.test,y,pop,x.fix,comm=1,p.cut.final,nboot,nchunk=1,b=1)
	unhit.snp = rownames(reg.unhit)[!is.na(reg.unhit[,1])]
	
	x.full = data.frame(x.final[,colnames(x.final) %in% c(hit.snp,unhit.snp),drop=F],check.names=F)
	name.full = colnames(x.full)
	
	# GOODNESS OF FIT ESTIMATES FOR VARIOUS MODELS
	RSS.h0 = sum((y-mean(y))^2)
	res.lm.full = lm(y ~ pop + ., data=x.full)
	RSS.full = sum(residuals(res.lm.full)^2)
	RSS.pop = sum(residuals(lm(y ~ pop))^2)
	RSS.snp = sum(residuals(lm(y ~ ., data=x.full))^2)

	
	fit.full = 1 - RSS.full/RSS.h0
	fit.pop = 1 - RSS.pop/RSS.h0
	fit.snp = 1 - RSS.snp/RSS.h0
	
	# R SQUARE ESTIMATES OF GWAS SNP BASED ON THE FULL MODEL
	Rsq.sinsnp = vector()
	for(i in 1:length(name.full)){
	mod.reduce = formula(paste('.~.',paste('`',name.full[i],'`',sep=''),sep='-'))
	reduce.lm = update(res.lm.full,mod.reduce,data=x.full)
	RSS.reduce = sum(reduce.lm$residuals^2)
	Rsq.sinsnp[i] = 1 - RSS.full/RSS.reduce
	}
	
	Predict.out = data.frame(pop=fit.pop, snp=fit.snp, full=fit.full)
	para.est = summary(res.lm.full)$coefficient
	plist = summary(res.lm.full)$coefficient[-c(1:length(unique(pop))),4]
	elist = summary(res.lm.full)$coefficient[-c(1:length(unique(pop))),1]
	names(plist) = gsub('`','',names(plist))
	snp.out = data.frame(bin_inf[match(names(plist),bin_inf$Tag_snp),],BPP=sig_bpp$BPP[match(names(plist),sig_bpp$Tag_snp)],
				pvalue=plist,effect=elist,Rsq=Rsq.sinsnp[match(names(plist),name.full)])
	
	reg_single = cbind(reg_single,hits=gwas.hit[match(reg_single$Tag_snp,colnames(x.final))]*1)
	sig_filt = cbind(sig_filt,hits=gwas.hit[match(sig_filt$Tag_snp,colnames(x.final))]*1)
	snp.out = cbind(snp.out,hits=gwas.hit[match(snp.out$Tag_snp,colnames(x.final))])
	qtl$pop = levels(phe$pop)[match(qtl$pop,seq(length(levels(phe$pop))))]
	qtl$hits = qtl.hit*1
	
	library(xlsx)
		
	write.xlsx(Cutoffs,paste(trait,'bpp.res.xlsx',sep = '.'),'Cutoffs',row.names = F,col.names=F)
	write.xlsx(reg_single,paste(trait,'bpp.res.xlsx',sep='.'),'SingleReg',row.names = F,append=T)
	write.xlsx(perm_bpp,paste(trait,'bpp.res.xlsx',sep = '.'),'Permbpp',row.names = F,append=T)
	write.xlsx(sig_bpp,paste(trait,'bpp.res.xlsx',sep = '.'),'SigSNPbpp',row.names = F,append=T)
	write.xlsx(sig_filt,paste(trait,'bpp.res.xlsx',sep= '.'),'SigSNPfilt',row.names = F,append=T)
	
	write.xlsx(snp.out,paste(trait,'bpp.res.xlsx',sep='.'),'SNP.finalmodel',row.names=F,append=T)
	write.xlsx(para.est,paste(trait,'bpp.res.xlsx',sep='.'),'Params.finalmodel',append=T)
	write.xlsx(Predict.out,paste(trait,'bpp.res.xlsx',sep='.'),'Fit.model',row.names=F,append=T)
	write.csv(bin_inf,paste(trait,'Genome.bin.info',sep='.'),row.names=F,quote=F)
	write.xlsx(qtl,paste(trait,'bpp.res.xlsx',sep='.'),'uniq.QTL',row.names=F,append=T)
	
cat('-----------------OUTPUT GWAS RESULTS INTO FOLLOWING SHEETS:--------------------------------------------------------','\n',
    '-----------------------------------------------------------SigSNPbpp------------------------------------------','\n',
    '-----------------------------------------------------------Permbpp--------------------------------------------','\n',
    '-----------------------------------------------------------SigSNPfilt-----------------------------------------','\n',
    '-----------------------------------------------------------SingleReg------------------------------------------','\n',
    '-----------------------------------------------------------SNP.finalmodel-------------------------------------','\n',
    '-----------------------------------------------------------Param.finalmodel-----------------------------------','\n',
    '-----------------------------------------------------------Fit.model------------------------------------------','\n',
    '-----------------------------------------------------------Cutoffs--------------------------------------------','\n',
    '-----------------------------------------------------------Genome.bin.info------------------------------------','\n',
    '-----------------------------------------------------------uniq.QTL-------------------------------------------','\n')
	
stop_time = Sys.time()
cat('\nTHE GWAS PROCEDURE START:\n')
print(start_time)
cat('THE GWAS PROCEDUURE END:\n')
print(stop_time)
}	

# SET THE FUNCTION FOR REGRESSION ANALYSIS USING VECTORIZATION METHOD
reg_vector = function(x,rss_h0,pop,j,ncof){
		RSS = apply(x,2,function(x){sum(lm.fit(as.matrix(x),rss_h0)$residuals^2)})
		RSS_H0 = sum(rss_h0^2)
		df_snp = 1
		df_res = nrow(x)-df_snp-ncof-j-length(unique(pop))+1
		p_value = pf((rep(RSS_H0,length(RSS))/RSS-1)*df_res/df_snp,df_snp,df_res,lower.tail=FALSE)
			
		return(p_value)
	}

# SET THE FUNCTION TO PERFORM THE FORWARD AND BACKWARD REGRESSION
stepwise_reg = function(x,y,pop,cofc,comm,p_cutoff,nboot,nchunk,b=1){

	x = data.frame(x,check.names=F)
	cofc = as.matrix(cofc)
	ncof = ifelse(ncol(cofc) == 1,0,sum(!duplicated(cofc[,apply(cofc,2,function(x)any(x != x[1]))],MARGIN=2)))
	stopifnot(length(y) == nrow(x))
	stopifnot(length(pop) == nrow(x))
	
	fwd_cof = 1
	p_min = 0
	p_fwd = NULL
	fwd_lm0 = lm(y ~ pop + cofc)
	
	# FORWARD REGRESSION STRAT TO OPTIMAZE THE MODEL VARIABLES
	while(p_min <= p_cutoff){
	fwd_cof = c(fwd_cof,ifelse(is.null(names(p_fwd)[which.min(p_fwd)]),1,paste('`',names(p_fwd)[which.min(p_fwd)[1]],'`',sep='')))
	fwd_mod = formula(paste('.~.',paste(fwd_cof,collapse='+'),sep='+'))
	fwd_lm = update(fwd_lm0,fwd_mod, data=x)
	rss_h0 = residuals(fwd_lm)
	
	if(nchunk == 1) p_fwd = reg_vector(x,rss_h0,pop,length(fwd_cof),ncof)
	else {
	snp_n = ncol(x)
	size = floor(snp_n / nchunk)
	nchunk_factor = c(rep(1:(nchunk-1), each = size), rep(nchunk, snp_n - size*(nchunk-1)))
	p_fwdn = lapply(tapply(1:snp_n, nchunk_factor, unique), function(n) reg_vector(x[,n],rss_h0,pop,length(fwd_cof),ncof))			
	p_fwd = unlist(p_fwdn)
	names(p_fwd)=colnames(x)
	rm(p_fwdn)}
	p_min = p_fwd[which.min(p_fwd)]	
	} 
	

	# BACKWARD REGRESSION START TO AVOID THE OVER-FIT OF THE MODEL	
	bwd_cof = 0
	p_max = 1
	p_bwd = NULL
	bwd_lm0 = fwd_lm
	
	while(p_max > p_cutoff){
	bwd_cof = c(bwd_cof,names(p_bwd)[which.max(p_bwd)])
	bwd_mod = formula(paste('.~.',paste(bwd_cof,collapse='-'),sep='-'))
	bwd_lm = update(bwd_lm0,bwd_mod,data=x)
	p_bwd = summary(bwd_lm)$coefficient[-c(1:(length(unique(pop))+ncof)),4]
	snp.id = rownames(summary(bwd_lm)$coefficient)[-c(1:(length(unique(pop))+ncof))]
	names(p_bwd) = snp.id
	if(length(p_bwd) > 0) p_max = p_bwd[which.max(p_bwd)]
	else break
	}
	
	if(p_max == 1){
	pvalue = c('NA' = NA)
	effect = c('NA' = NA)
	} else{
	pvalue = p_bwd
	effect = summary(bwd_lm)$coefficient[-c(1:(length(unique(pop))+ncof)),1]
	names(effect) = snp.id
	names(pvalue) = gsub('`','',names(pvalue))
	names(effect) = gsub('`','',names(effect))
	stopifnot(all.equal(names(pvalue),names(effect)))
	}
	
	if(comm == 3) cat('---------------BOOTSTRAP',b,'OF',nboot,'DONE!',sum(!is.na(pvalue)),'SNPS IN THE FINAL MODEL----------------------------','\n')
	else  if(comm == 2 & b %% 10 == 0) cat('------------BOOTSTRP FOR PERMUTATION',b,'OF',nboot,'DONE!------------------------------------------------','\n')
	return(cbind(pvalue,effect))
	}

# PERMUTATION TESTS FOR EVALUATING THE FALSE POSITIVES IN BOOTSTRAP PROGRESS
perm_reg = function(p,x,y,pop,cofc,nchunk){
			if(p %% 10 == 0) cat('-----------------PERMUTATION',p,'START!--------------------------------------------------------------------------','\n')
			set.seed(p)
			id = unlist(lapply(tapply(1:nrow(x),pop,unique),function(m)sample(m,length(m),replace=F)))
			pop_perm = pop
			y_perm = y[id]
			rss_h0 = residuals(lm(y_perm ~ pop_perm + cofc))
			ncof = ifelse(ncol(cofc) == 1,0,sum(!duplicated(cofc[,apply(cofc,2,function(x)any(x != x[1]))],MARGIN=2)))
			
			if(nchunk == 1) ps = reg_vector(x,rss_h0,pop,1,ncof)
			else {
				snp_n = ncol(x)
				size = floor(snp_n/nchunk)
				nchunk_factor = c(rep(1:(nchunk-1),each = size),rep(nchunk,snp_n-size*(nchunk-1)))
				
				ps_perm = lapply(tapply(1:snp_n,nchunk_factor,unique),function(n)reg_vector(x[,n],rss_h0,pop,1,ncof))				
				ps = unlist(ps_perm)
				rm(ps_perm)
				gc()
			}

		perm_p = min(ps,na.rm=T)
		return(perm_p)
	}

# SET THE FUNCTION FOR BOOTSTRAP ANALYSIS BY WHICH THE EQUAL NUMBER OF INDIVIDUALS WERE REPEATALLY SAMPLED
nbootf = function(x,y,pop,cofc,comm,p_cutoff,nboot,nchunk,mc){
	boot.fun = function(b){
				set.seed(b)
				id_pop = lapply(tapply(1:nrow(x),pop,unique),function(m)sample(m,0.8*length(m),replace=F))
				id = unlist(id_pop)
				x_b = x[id,]
				cofc_b = as.matrix(cofc[id,])
				
				pop_b = pop[id]
				y_b = y[id]
				res = stepwise_reg(x_b,y_b,pop_b,cofc_b,comm,p_cutoff,nboot,nchunk,b)
				
				rm(id,x_b,y_b,pop_b,cofc_b)
				gc()
				colnames(res) = paste(colnames(res),paste('boot',b,sep = ''),sep = '')
				cbind(rownames(res),res)
				}
		return(mclapply(1:nboot,boot.fun,mc.cores = mc))
			}

# PERMUTATION TESTS DETERMINE THE OPTIMAL 'BPP' CUTOFF FOR BOOTSTRAP ANALYSIS
perm_boot = function(p,x,y,pop,cofc,bin_inf,comm,p_cutoff,nboot,nperm_boot,nchunk,mc){
			cat('------------------PERMUTATION TEST',p,'OF',nperm_boot,'START!-----------------------------------------------','\n')
			set.seed(p)
			y_perm = y[unlist(lapply(tapply(1:nrow(x),pop,unique),function(m)sample(m,length(m),replace = F)))]
			boot_ = nbootf(x,y_perm,pop,cofc,comm,p_cutoff,nboot,nchunk,mc)
			perm_boot = boot_summary(boot_,nboot,bin_inf)
			
			return(perm_boot)
	}

# SUMMARIZE THE BOOTSTRAP RESULT LIST INTO FINAL OUTPUTS
boot_summary = function(boot.list,nboot,bin_inf){
	boot_t = as.matrix(Reduce(function(...)merge(...,all = T),boot.list))
	p_res = matrix(as.numeric(matrix(boot_t[,-1],nrow(boot_t),2*nboot)[,1:(2*nboot) %% 2 != 0]),nrow = nrow(boot_t))
	e_res = matrix(as.numeric(matrix(boot_t[,-1],nrow(boot_t),2*nboot)[,1:(2*nboot) %% 2 == 0]),nrow = nrow(boot_t))
	BPP = apply(p_res,1,function(x) {sum(!is.na(x))/nboot})
	p = apply(p_res,1,function(x) median(x,na.rm = T))
	effect = apply(e_res,1,function(x) median(x,na.rm = T))
	sig_boot = cbind(bin_inf[match(boot_t[,1],bin_inf$Tag_snp),],BPP,p,effect)
	sig_boot = sig_boot[!is.na(sig_boot$Tag_snp),]
	return(sig_boot)
}

# EXTRACT THE CONTINUOUS BIN INFORMATION USED IN GWAS SCANING
genome_bin = function(inf,inf.uni){
	
uniq.id = which(inf$snp %in% inf.uni$snp)
chr_max = tapply(1:nrow(inf),inf$chr,max)
id_sta = uniq.id
id_end = c(uniq.id[-1]-1,NA)
id_end[tapply(1:nrow(inf.uni),inf.uni$chr,max)]=chr_max
pos_sta = inf$pos[id_sta]
pos_end = inf$pos[id_end]
bin = data.frame(Bin_num=1:nrow(inf.uni),Tag_snp=inf.uni$snp,
	Tag_alle=inf.uni$alleles,Chr=inf.uni$chr,Pos_sta=pos_sta,Pos_end=pos_end)

library(plyr)

bin_inf = ddply(bin,'Chr',transform,Pos_sta_just=floor(c(0,(Pos_sta[-1]+Pos_end[-length(Pos_end)])/2)),
		Pos_end_just=ceiling(c((Pos_sta[-1]+Pos_end[-length(Pos_end)])/2,Pos_end[length(Pos_end)])))
bin_inf$Len=with(bin_inf,Pos_end_just-Pos_sta_just)
return(bin_inf)
}





