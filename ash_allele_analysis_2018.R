if (.Platform$OS.type=="windows"){
	setwd("C:/Users/seanmh/Dropbox/Projects/MSB_Kew_ash_Tree_Seed/")
	setwd("C:/Users/shoban/Dropbox/Projects/MSB_Kew_ash_Tree_Seed/") }
if (.Platform$OS.type=="unix") {
	setwd("/home/user/Dropbox/Projects/MSB_Kew_ash_Tree_Seed/")
	#setwd("/media/sean/Windows7_OS/Users/seanmh/Dropbox/Projects/MSB_Kew_ash_Tree_Seed/") 
}
	
library(adegenet); library(diveRsity)
library(parallel)
#library(foreach); 
#library(doMC); registerDoMC()	#OR..
#library(doParallel); cl<-makeCluster(8); registerDoParallel(cl)
source("Simulations_and_Code/src/sample_funcs.R")

colMax <- function(data) sapply(data, max, na.rm = TRUE)

thresh_freq_H<- 0.20; thresh_freq_L<- 0.05
	
thisgrid<-"large"

if (thisgrid=="small"){
scenarios_ran<-list.dirs(path = "./Scenarios_Nov_28_2016/10rep_10mark", full.names = TRUE,recursive=F)
SIZE_OF_GRID<-209;		NUM_GRID_ROWS<-19;		N_REGIONS<-15
#first_row_region<-c(1,3,5,7,9,11,13,15,17);			last_row_region<-c(2,4,6,8,10,12,14,16,19)
first_row_region<-	c(1,3,5,7,8,9,10,11,12,13,14,15,16,17,18)
last_row_region<-	c(2,4,6,7,8,9,10,11,12,13,14,15,16,17,19)
}
if (thisgrid=="medium"){
scenarios_ran<-list.dirs(path = "./Scenarios_Nov_22_medium", full.names = TRUE,recursive=F)
SIZE_OF_GRID<-1215;		NUM_GRID_ROWS<-45;		N_REGIONS<-15
first_row_region<-	c(1,3,6,9,12,15,18,21,24,27,30,33,36,39,42)
last_row_region<-	c(2,5,8,11,14,17,20,23,26,29,32,35,38,41,45)	
}
if (thisgrid=="large"){
scenarios_ran<-list.dirs(path = "./Scenarios_Nov_22_large", full.names = TRUE,recursive=F)
SIZE_OF_GRID<-4717;		NUM_GRID_ROWS<-89;		N_REGIONS<-18
first_row_region<-	c(1,6,11,16,21,26,31,36,41,46,51,56,61,66,71,76,81,86)		
last_row_region<-	c(5,10,15,20,25,30,35,40,45,50,55,60,65,70,75,80,85,89)	
}		

scenarios_ran<-list.dirs(path = "./Simulations_and_Code/test2018", full.names = TRUE,recursive=F)

#Set up file path where scenario is
num_scen<-length(scenarios_ran)
num_reps<-1000

#set up results
#Last 4 are MSB, NTSB, 10*MSB, 10*NTSB, TOTAL
summ_results<-array(dim=c(1485+4+1,8,num_scen,num_reps))
colnames(summ_results)<-c("total plants","total populations","G", "GLF", "GR", "RC", "LC", "LR")
type_samp<-c("random","each fr diff reg", "only N 2 rows", "center 2 rows", "only S 2 rows", 
					 "core", "edge", "focus S", "focus N")
rownames(summ_results)<-c(rep(type_samp,each=165),"MSB","NTSB","MSB10","NTSB10","TOTAL")

######################################
#---FILE CHECK AND FILE CONVERSION---#
######################################

#for (scen in 1:length(scenarios_ran)){	
#	#Check for and remove genind
#	gen_files<-dir(path=scenarios_ran[scen],pattern="gen")
#	if (length(gen_files)!=0)  file.remove(file.path(scenarios_ran[scen],gen_files))
#	#convert to genind	
#	reps_ran_arp<-list.files(scenarios_ran[scen], pattern="arp")
#	arp_file_list<-file.path(scenarios_ran[scen],reps_ran_arp,sep="")
#	if (.Platform$OS.type=="unix") arp_file_list<-substr(arp_file_list,1,nchar(arp_file_list)-1)
#	mclapply(arp_file_list,arp2gen,mc.cores=16)
#	#foreach (nrep = 1:length(reps_ran_arp)) %dopar% {	arp2gen(arp_file_list[nrep])	}		#alternative MC loop
#}


##############################
#---DEFINE REGIONAL MAKEUP---#
#---AND WHERE MSB SAMPLED----#
##############################
scen<-1
region_makeup<-set.regions(scenarios_ran[scen], SIZE_OF_GRID, NUM_GRID_ROWS, N_REGIONS, 
						   first_row_region, last_row_region)
		
file_MSB_all<-"Data_from_simon/sampled_so_far/2018/sampled_locations_all_2018.csv"
file_NTSB<-"Data_from_simon/sampled_so_far/2018/sampled_locations_NTSB_2018.csv"
MSB_locations<-MSB.samp(get.uk.grid(scenarios_ran[scen], SIZE_OF_GRID, NUM_GRID_ROWS), file_MSB_all)
MSB_locations<-MSB_locations[MSB_locations!=0]
MSB_plant_nums<-read.csv(file_MSB_all)[,4]*10
MSB_plant_nums[MSB_plant_nums>=250]<-240
NTSB_locations<-MSB.samp(get.uk.grid(scenarios_ran[scen], SIZE_OF_GRID, NUM_GRID_ROWS), file_NTSB)
NTSB_locations<-NTSB_locations[NTSB_locations!=0]
NTSB_plant_nums<-read.csv(file_NTSB)[,4]*10
NTSB_plant_nums[NTSB_plant_nums>=250]<-240
	
	
	
#####################
#---MAIN ANALYSIS---#
#####################

for (scen in 1:length(scenarios_ran)){	
	reps_ran_gen<-list.files(scenarios_ran[scen], pattern="gen")
	for (nrep in 1:num_reps){
		print(scenarios_ran[scen])
		
		#make a genind (by individual) and genpop (by population)
		temp_file_name<-file.path(scenarios_ran[scen],reps_ran_gen[nrep],sep="")
			if (.Platform$OS.type=="unix")  temp_file_name<-substr(temp_file_name,1,nchar(temp_file_name)-1)
		UK_genind<-read.genepop(temp_file_name,ncode=3)
		UK_genpop<-genind2genpop(UK_genind)
		

		#--NUMBER OF POPULATIONS, INDIVIDUALS, REGIONS, REGIONAL MAKEUP--#
		n_pops<-length(levels(UK_genind@pop))
		n_total_indivs<- length(UK_genind@tab[,1])
		n_ind_p_pop<-table(UK_genind@pop)
		allele_freqs<-colSums(UK_genpop@tab)/(n_total_indivs*2)	
		
		
		#######################
		#---#DETERMINE WHAT ALLELES FALL IN WHAT CATEGORIES---#
		#######################
		allele_cat<-get.allele.cat(UK_genpop, region_makeup, N_REGIONS, n_ind_p_pop)	
		glob_com<-allele_cat[[1]]; 		glob_lowfr<-allele_cat[[2]];		glob_rare<-allele_cat[[3]]
		reg_com_int<-allele_cat[[4]];	loc_com_int<-allele_cat[[5]]; 		loc_rare<-allele_cat[[6]]
		
	
		#######################
		#--SUBSAMPLE SECTION
		#######################
		
		### LIST OF POPULATIONS SAMPLED
		#sample certain number of populations, and individuals, and by region
		center_pops_vect<-read.csv("Grids/center_edge/center_pops.txt",sep=",",header=F)[[1]]
		edge_pops_vect<-read.csv("Grids/center_edge/edge_pops.txt",sep=",",header=F)[[1]]
		type_samp<-c("random","each fr diff reg", "only N 2 rows", "center 2 rows", "only S 2 rows", 
					 "core", "edge", "focus S", "focus N")
		l_plt_smp<-length(c(2,seq(5,50, by=5)))
		n_pops_to_samp<- rep(c(rep(2,l_plt_smp),rep(5,l_plt_smp),rep(10,l_plt_smp),rep(15,l_plt_smp),
							  rep(20,l_plt_smp),rep(25,l_plt_smp),rep(30,l_plt_smp),rep(35,l_plt_smp),
							  rep(40,l_plt_smp),rep(45,l_plt_smp),rep(50,l_plt_smp),rep(55,l_plt_smp),
							  rep(60,l_plt_smp),rep(65,l_plt_smp),rep(70,l_plt_smp)),length(type_samp))  #rep(n_pops,l_plt_smp),
		N_SAMPS_P_POP<-as.list(rep(c(2,seq(5,50, by=5)),(length(type_samp)*length(unique(n_pops_to_samp)))))
		POPS_TO_SAMP<-list();		this_slot<-1
		for (t in 1:length(type_samp)) for (p in 1:length(unique(n_pops_to_samp))) {
				if (p>6) REPLC_N=T; if (p<=6) REPLC_N=F
				if (t==1)	the_pops<-sample(1:n_pops, unique(n_pops_to_samp)[p])
						#t=2 makes list of 1 pop per region then samples n_pops from that 
				if (t==2)	the_pops<-sample(sapply(region_makeup,sample,3), unique(n_pops_to_samp)[p],replace=REPLC_N)
						#t=3,4,5 takes the top, middle, & end of the grid, unlists, then samples n_pops from that 
				if (t==3)	the_pops<-sample(unlist(region_makeup[1:2]),unique(n_pops_to_samp)[p],replace=REPLC_N)
				if (t==4)	{
						temp_mid<-length(region_makeup)/2
						the_pops<-sample(unlist(region_makeup[(temp_mid-1):(temp_mid+1)]),unique(n_pops_to_samp)[p])
					}
				if (t==5)	the_pops<-sample(unlist(tail(region_makeup)[4:5]),unique(n_pops_to_samp)[p])	
				if (t==6)	the_pops<-sample(center_pops_vect, unique(n_pops_to_samp)[p])
				if (t==7)	the_pops<-sample(edge_pops_vect, unique(n_pops_to_samp)[p])
				if (t==8)	{
						num_in_focus<-floor(0.75*unique(n_pops_to_samp)[p])
						the_pops<-c(sample(unlist(tail(region_makeup)[4:5]),num_in_focus),
									sample(1:n_pops, (unique(n_pops_to_samp)[p]-num_in_focus)))
					}
				if (t==9)	{
						num_in_focus<-floor(0.75*unique(n_pops_to_samp)[p])
						the_pops<-c(sample(unlist(region_makeup[1:2]),num_in_focus, replace=REPLC_N),
									sample(1:n_pops, (unique(n_pops_to_samp)[p]-num_in_focus)))
					}
				
				for (x in 1:l_plt_smp) { POPS_TO_SAMP[[this_slot]]<-the_pops; this_slot=this_slot+1}
			}
		for (i in 1:length(POPS_TO_SAMP)) 
			N_SAMPS_P_POP[[i]]<-unlist(rep(N_SAMPS_P_POP[i],length(POPS_TO_SAMP[[i]])))
				#This checks the above code- just swap out "34" for larger numbers
				#temp_ind<-vector(length=15)
				#for (i in 1:15)	temp_ind[i]<-(which(UK_grid==POPS_TO_SAMP[[34]][i],arr.ind=T)[1])
				#sort(temp_ind)
			
		#This is for the MSB samples
			POPS_TO_SAMP[[this_slot]]<-MSB_locations; N_SAMPS_P_POP[[this_slot]]<-MSB_plant_nums/10
			N_SAMPS_P_POP[[this_slot]][N_SAMPS_P_POP[[this_slot]][]==1]<-2; this_slot=this_slot+1
			POPS_TO_SAMP[[this_slot]]<-NTSB_locations; N_SAMPS_P_POP[[this_slot]]<-NTSB_plant_nums/10
			N_SAMPS_P_POP[[this_slot]][N_SAMPS_P_POP[[this_slot]][]==1]<-2; this_slot=this_slot+1
			POPS_TO_SAMP[[this_slot]]<-MSB_locations; N_SAMPS_P_POP[[this_slot]]<-MSB_plant_nums
			this_slot=this_slot+1
			POPS_TO_SAMP[[this_slot]]<-NTSB_locations; N_SAMPS_P_POP[[this_slot]]<-NTSB_plant_nums
		
	UK_genind_sep<-seppop(UK_genind)
		for (samp in 1:length(POPS_TO_SAMP)){
			alleles<-matrix(nrow=length(POPS_TO_SAMP[[samp]]),ncol=length(allele_freqs))
			alleles<-sample.pop(UK_genind_sep,POPS_TO_SAMP[[samp]],N_SAMPS_P_POP[[samp]])
			summ_results[samp,1,scen,nrep]<-sum(N_SAMPS_P_POP[[samp]]); summ_results[samp,2,scen,nrep]<-length(N_SAMPS_P_POP[[samp]])

  
			#NOW SEE WHAT IS CAUGHT
			#this will record each time an allele is in one of our sampled populations
			for (p in 1:n_pops){
				caught_glob<-colSums(alleles,na.rm=T)
				caught_glob_lowfr<-colSums(alleles,na.rm=T)[glob_lowfr]
				caught_glob_rare<-colSums(alleles,na.rm=T)[glob_rare]
				caught_reg_com<-colSums(alleles,na.rm=T)[reg_com_int]
				caught_loc_com<-colSums(alleles,na.rm=T)[loc_com_int]
				caught_loc_rare<-colSums(alleles,na.rm=T)[loc_rare]
				}
			#as long as its not missed (0) it counts as being caught (could try other minima)
			got_G<-sum(caught_glob>=1);  got_GLF<-sum(caught_glob_lowfr>=1);	got_GR<-sum(caught_glob_rare>=1)
			got_RC<-sum(caught_reg_com>=1);	got_LC<-sum(caught_loc_com>=1);   	got_LR<-sum(caught_loc_rare>=1)
			
			summ_results[samp,3,scen,nrep]<-got_G; summ_results[samp,4,scen,nrep]<-got_GLF
			summ_results[samp,5,scen,nrep]<-got_GR; summ_results[samp,6,scen,nrep]<-got_RC
			summ_results[samp,7,scen,nrep]<-got_LC; summ_results[samp,8,scen,nrep]<-got_LR
			}
		summ_results[1490,3:8,scen,nrep]<-c(length(allele_freqs),length(glob_lowfr),length(glob_rare),
										   length(reg_com_int),length(loc_com_int),length(loc_rare))
		
	}
			#difference between the scenarios_ran
			#apply(summ_results[,,1,1:10]-summ_results[,,2,1:10],c(1,2),mean)
			#apply(summ_results[1089:1092,,1,1:num_reps],c(1,2),mean)
}
			
save(summ_results,file="summ_results_mina_2_14_18_MSB2.Rdata")
					


#####################################
#			RESULTS					#
#####################################
	
#FOR THE CHOSEN SCENARIO
if (.Platform$OS.type=="windows"){
	setwd("C:/Users/shoban.DESKTOP-DLPV5IJ/Dropbox/Projects/MSB_Kew_ash_Tree_Seed/")
	setwd("C:/Users/shoban/Dropbox/Projects/MSB_Kew_ash_Tree_Seed/") }
if (.Platform$OS.type=="unix") 	setwd("/home/user/Dropbox/Projects/MSB_Kew_ash_Tree_Seed/")
load(file="summ_results_2_14_18_MSB2.Rdata")

load(file="summ_results_mina_2_16_18_MSB2.Rdata")

type_samp<-c("random","each fr diff reg", "only N 2 rows", "center 2 rows", "only S 2 rows", 
					 "core", "edge", "focus S", "focus N")
type_allele<-c("global", "global low frequency", "global rare", "regional", "local common", "local rare")
num_samp_strat<-length(summ_results[,1,1,1]); num_mom_pop_comb<-11*15
num_trees<-c(2,seq(5,50,5))
num_reps<-1000



###########################################################
#			COMPARE SPATIAL STRATEGIES GRAPH					#
###########################################################
		
#The first question is about ranking of the large scale spatial strategies relative to each other
#This uses 'matplot' to plot the different sampling strategies 

total_exists<-apply(summ_results[num_samp_strat,3:8,,1:num_reps],1,mean)
pdf(file="compare_spatial.pdf",width=14,height=11)
par(mfrow=c(3,2),oma=c(4,3,4,3),mar=c(3,3,2,2))
#this will show 2 populations (lines 1:11), 20 populations (line 45) and 45 populations (line 100)
#If I wanted 5 populations sampled I would use 12
#and we will focus on allele category 3 (global) and 7 (locally common)
#FOR SUPPLEMENTAL
for (start_pop in c(1,45,100)){
	#we need to get, for example, 2 populations for each of the sampling strategies (all nine of them)
	#each spatial sampling strategy has 121 (or 11 times 11) slots for all tree/pop combinations
	#so essentially need slots 1:11, then (1+121):(11+121) etc. etc.  The next three lines get that list
	start_each_samp<-seq(start_pop,num_samp_strat-5,by=num_mom_pop_comb)
	for (i in 1:(length(num_trees)-1)) unique(start_each_samp<-c(start_each_samp,start_each_samp+1))
	some_samps<-sort(unique(start_each_samp))
		#For the two allele categories of focus, get the results for that list, make a matrix with 
		#the spatial strategies as columns and adding more individuals as rows, then plot this
	for (i in c(3,7)) {
		select_samp<-apply(summ_results[some_samps,i,,1:num_reps],1,mean,na.rm=T)
		matplot(matrix(select_samp/total_exists[i-2],nrow=length(num_trees)),type="l",lwd=2,
			col=c(rep("black",2),rep("red",3),rep("blue",5)), lty=c(1,2,1,2,3,1,2,3,4,5),
			   ylab="",xaxt="n",xlab="",cex.axis=1.8)
		axis(4,labels=F); axis(1,labels=num_trees, at=1:11,cex.axis=1.8) }
mtext(side=4,line=2,paste(summ_results[start_pop,2,1,1]," populations sampled"),cex=1.4) }
mtext(side=1,line=1.3,"number of trees sampled per population",outer=T,cex=1.4)
mtext(side=2, line=1.3,"proportion of alleles preserved",outer=T,cex=1.4)
mtext(side=3,"    global alleles                                           locally common alleles",outer=T,cex=2.1)
legend(6,.8,legend=type_samp, col=c(rep("black",2),rep("red",3),rep("blue",5)), lty=c(1,2,1,2,3,1,2,3,4,5), cex=1.3,bg="white", lwd=2)
dev.off()

#this as above but just the global alleles and assumes 40 populations
#FOR MAIN MANUSCRIPT
pdf(file="compare_spatial_global_40.pdf",width=9,height=7)
select_samp<-apply(summ_results[some_samps,3,,1:num_reps],1,mean,na.rm=T)
matplot(matrix(select_samp/total_exists[1],nrow=length(num_trees)),type="l",lwd=2.5,	
		col=c(rep("black",2),rep("red",3),rep("blue",5)), lty=c(1,2,1,2,3,1,2,3,4,5), 
		ylab="proportion of alleles captured",xaxt="n",xlab="number of trees sampled", 
		cex=1.5,cex.axis=1.5,cex.lab=1.5)
legend(6,.71,legend=type_samp, col=c(rep("black",2),rep("red",3),rep("blue",5)), lty=c(1,2,1,2,3,1,2,3,4,5), cex=1.3,bg="white", lwd=2.5)
axis(1,at=1:11,labels=(c(2,seq(5,50, by=5))), cex.axis=1.5)
dev.off()




###########################################################
#			COMPARE SPATIAL STRATEGIES TABLE					#
###########################################################	
#This analysis is to rank the spatial strategies
#NOT REALLY USED IN PAPER- THE PLOT TELLS ALL THAT WE NEED
library("tidyr"); library("reshape2")

#melt and dcast will swap L and R, so have to fix order
total_exists<-apply(summ_results[num_samp_strat,3:8,,1:100],1,mean)
total_exists_ord<-total_exists; total_exists_ord[4:5]<-total_exists[5:6]; total_exists_ord[6]<-total_exists[4]

#first rbind the reps together and cbind a column to add sample strategies to each row
bound_results<-as.data.frame(summ_results[,,1,1])
for (i in 2:100)	bound_results<-rbind(bound_results,summ_results[,,1,i])
labels_one_rep<-as.factor(c(rep(type_samp,each=num_mom_pop_comb),"MSB","NTSB","MSB10","NTSB10","TOTAL"))
bound_results<-cbind(rep(labels_one_rep,100),bound_results)
colnames(bound_results)<-c("strat","plants","pops","G","GLF","GR","RC","LC","LR")

melted_results<-gather(bound_results,allcat,alleles,G,GLF,GR,RC,LC,LR,-pops,-plants,-strat)

#first get MSB results
MSB_results<-melted_results[which(melted_results$strat=="MSB"),];	means_MSB<-dcast(MSB_results,strat~allcat,mean)
means_MSB[2:7]/total_exists_ord

#remove MSB and total for now
melted_results<-melted_results[-which(melted_results$strat=="MSB"),]; 	melted_results<-melted_results[-which(melted_results$strat=="MSB10"),]
melted_results<-melted_results[-which(melted_results$strat=="NTSB"),]; 	melted_results<-melted_results[-which(melted_results$strat=="NTSB10"),]
melted_results<-melted_results[-which(melted_results$strat=="TOTAL"),]

means_all<-dcast(melted_results,strat~allcat,mean)

ss_ranking<-function(num_plants,num_pops){
	spec_mean<-dcast(melted_results[melted_results$plants==num_plants&melted_results$pops==num_pops,],strat~allcat,mean)	#50/ 50
	rownames(spec_mean)<-spec_mean[,1]; spec_mean<-spec_mean[,-1]
	#print(rownames(t(t(spec_mean)/total_exists_ord)[order(spec_mean[,4]),])[5:9])
		#identify the best by ranking
	write.table(cbind(rownames(t(t(spec_mean)/total_exists_ord)[order(spec_mean[,1]),])[5:9],
		t(t(spec_mean)/total_exists_ord)[order(spec_mean[,1]),][5:9]),
		file="ss_ranking2.csv", append=T)
	write.table( t(t(spec_mean)/total_exists_ord), file="ss_ranking.csv", append=T)
}
#do this for low, medium, high sampling, more pops, more individuals see if ranking changes
temp_plants<-c("2500","625","200","50","250","250"); temp_pops<-c("50","25","10","5","5","50")
for (i in 1:length(temp_plants)) ss_ranking(temp_plants[i],temp_pops[i])
	
#determine which are significantly different than each other
anov_res<-aov(alleles~strat+plants+pops+allcat,data=melted_results)
anov_res<-aov(alleles~strat,data=melted_results)
Tuk_res<-TukeyHSD(anov_res); 	Tuk_res$strat[order(Tuk_res$strat[,4]),]
#Interestingly the edge and core dont significantly differ from each other; nearly all the rest differ
barplot(means_all[,2]/total_exists_ord[1],names.arg=means_all[,1],las=2)
	
	
	
			
###########################################################
#			MSB vs. POTENTIAL SAMPLING					#
###########################################################
#Compare MSB to other samplings, see which of the strategies did better at proportion of alleles
#and/ or how many samples does it take to beat MSB, with a diff strategy (i.e. north)?
	
summ_results_ex<-summ_results
means_caught<-apply(summ_results_ex[1:1490,,1,1:1000],c(1,2),mean)
means_caught[1490,1:2]<-1
prop_caught<-t(t(means_caught)/means_caught[1490,])
write.csv(prop_caught,file="prop_caught.csv")
#which sampling strategies could have beaten the MSB in terms of choice of populations
which_beat_MSB<-prop_caught[which(prop_caught[,3]>prop_caught[1487,3]),]
write.csv(which_beat_MSB,"which_beat_MSB_G.csv")
which_beat_MSB<-prop_caught[which(prop_caught[,7]>prop_caught[1487,7]),]
write.csv(which_beat_MSB,"which_beat_MSB_LC.csv")
#to get as much as they did, if only sampling in the south, they would have had to sample 
#very hard to capture all the alleles!
	
	
	
pdf(file="local_vs_global_accum.pdf")
plot(prop_caught[1:1093,3],prop_caught[1:1093,7],xlim=c(0,1),ylim=c(0,1), xlab="locally common alleles", ylab="all alleles")
	abline(0,1,col="green",lwd=2)
dev.off()
prop_caught[as.numeric(which(prop_caught[,3]>prop_caught[,7])),1:2]
#Get global alleles "faster" for small collections, and local alleles faster with big collections- 
#this means it is easier to get all the local common alleles because the accumulation of global ones slows one




###########################################################
#			COMPARE TREES VS. POP'NS					#
###########################################################
#The question is it better to sample more trees or more populations can partly be answered with a graph
#This makes two plots, one for number of trees on the X axis and one for number of populations on the X axis
#With additional lines being the other variable (populations and trees)
	#The following loop will do this for every allele type
	for (A in 3:8){
pdf(paste(A,"trees_vs_popns.pdf"), height=5,width=8)
all_mean<-apply(summ_results[,,1,1:1000],c(1,2),mean)
n_pops_samp<- c(2,seq(5,70,by=5)); n_trees_samp<- c(2,seq(5,50,by=5))
all_mean_glob<-matrix(all_mean[1:num_mom_pop_comb,A],nrow=11,dimnames=list(n_trees_samp,n_pops_samp))/all_mean[num_samp_strat,A]
#plots of adding more trees and more populations
par(mfrow=c(1,2),mar=c(5,5,2,1))
plot(all_mean_glob[1,]~colnames(all_mean_glob),type="l",ylim=c(0,1),xlim=c(-4,50),xlab="number of populations",ylab="proportion of genetic variation")
for (t in 1:11) lines(all_mean_glob[t,]~colnames(all_mean_glob))
for (pop_index in c(1:4,6,11)) text(-3,all_mean_glob[pop_index,1],n_pops_samp[pop_index],cex=1.1)
axis(4, at= c(0,.2,.4,.6,.8,1), labels=F)
#plot(all_mean_glob[1,]~colnames(all_mean_glob),type="l",ylim=c(0.4,.9),xlab="number of populations",ylab="proportion of genetic variation")
#for (t in 2:10) lines(all_mean_glob[t,]~colnames(all_mean_glob))
par(mar=c(5,2,2,4))
plot(all_mean_glob[,1]~rownames(all_mean_glob),type="l",ylim=c(0,1),xlim=c(-5,50),xlab="number of trees",yaxt="n")
for (t in 1:11) lines(all_mean_glob[,t]~rownames(all_mean_glob))
for (mom_index in c(1:4,6,11)) text(-3,all_mean_glob[1,mom_index],n_trees_samp[mom_index])
axis(4, at= c(0,.2,.4,.6,.8,1), labels=T); axis(2, at= c(0,.2,.4,.6,.8,1), labels=F)
dev.off()
	}
#So the first black bar is 5 trees/ 35 populations, the third red bar is 10 populations/ 35 trees
#Second black bar is 10 trees/ 35 populations, the fourth red bar is 15 populations/ 35 trees

	
	
#########################################################################
#			GET THE DIMINISHING RETURNS POINT (PLATEAU 1%)					#
#########################################################################	
#ALSO CALL IT THE STOPPING POINT
#this looks at more trees or more pops
#it will take in the all_mean (all alleles means) matrix 
#first for trees it will go through the all_mean matrix and calculate the difference between row r+1 and r 
#then will divide this by the number of trees sampled from r+1 to r e.g. divide by 5 or 10 trees
#then this is stored in tree_diff
#then we determine for every sampling group (2 to 50 trees, 11 types) we find the first diff less then thresh
#we do skip the first minimum sampling (2) because it is the diff from the last sampling (50)
	#to look at other allele categories just replace the 3]<thresh below with other categories

thresh<-0.001

#first get the proportions
all_mean[,3:8]<-t(t(all_mean[,3:8])/(all_mean[1489,3:8]))
#add a column that is number of plants sampled
all_mean<-cbind(all_mean,all_mean[,1]/all_mean[,2])
diff<-all_mean[1:1489,]
#calculate the difference from the next closest sampling
for (i in 1:length(diff[,1])) for (c in 3:8) diff[i,c]<-(all_mean[i+1,c]-all_mean[i,c])/(all_mean[i+1,9]-all_mean[i,9])	
#mean(diff[diff[,3]<0.005,9])	#this is wrong

#Note that all_mean is sorted by number of populations so we can look at difference as we add more trees
	#The 99 is 9 sampling spatial strategies * 11 tree possibilities
b<-0;	plateau_tree<-vector(length=99)
#Go through each set of 11, identify which of those is less than thresh, take the minimum of that 
for (i in 0:98)	plateau_tree[i+1]<-diff[1+i*10+i+min(which(diff[(1+i*10+i+1):(11+i*10+i),3]<thresh)),9]
	mean(plateau_tree,na.rm=T)
	#0.001 tree = 28.1; 0.005 = 11.9; 0.01 = 7.22

#Reorder all_mean by number of trees
all_mean_temp<-all_mean
all_mean_temp<-all_mean_temp[order(rownames(all_mean_temp),all_mean_temp[,9]),]
diff<-all_mean_temp[1:1489,]
#calculate the difference from the next closest sampling
for (i in 1:length(diff[,1])) for (c in 3:8) diff[i,c]<-(all_mean_temp[i+1,c]-all_mean_temp[i,c])/(all_mean[i+1,1]-all_mean[i,1])	

#so we can look at difference as we add more populations
	#The 135 is 9 sampling spatial strategies * 15 tree possibilities
b<-0;	plateau_pop<-vector(length=135)
for (i in 0:134)	plateau_pop[i+1]<-diff[1+i*10+i+min(which(diff[(1+i*10+i+1):(11+i*10+i),3]<thresh)),2]
	mean(plateau_pop,na.rm=T)
	#0.005 = 17.1 populations, 0.001 = 20.4 populations



 
	
##########	
#JUNKYARD
##########
#separate data into single populations
	#UK_genind_sep<-lapply(seppop(UK_genind), function(x) x[sample(1:nrow(x$tab), 10)])
#this will look at what is captured in a given population
	#sum(as.vector(colSums(UK_genind_sep[[39]]@tab)[glob_com])==0)
	#sum(as.vector(colSums(UK_genind_sep[[39]]@tab)[glob_lowfr])==0)
	#sum(as.vector(colSums(UK_genind_sep[[39]]@tab)[glob_rare])==0)
#but really I want what is captured in all sampled populations
	#pooled_data<-repool(UK_genind_sep)
	#sum(as.vector(colSums(pooled_data_north@tab)[glob_com])==0)	
#this doesn't work to repool, because it makes new allele counts



#colSums lists all alleles and the number of populations (N) having 0 copies of that allele
#we then pull out the indices (which) of alleles for which N = number of populations - 1 
#(only one pop'n does not have 0 copies)
these<-as.vector(which(colSums(UK_genpop@tab<12)==(n_pops-1)))
#then take each of these alleles and sort the counts per population- how many counts are in that pop'n
sapply(these, function(x) sort(UK_genpop@tab[,x]))
big_enough<-table_counts[39,]>=30
these<-as.vector(which(colSums(UK_genpop@tab<3)==(n_pops-1)))
sapply(these, function(x) sort(UK_genpop@tab[,x]))
#CYCLE THROUGH THIS FROM 0 TO ABOUT 14, record the population and the allele, its count in that pop, and count in next pop

sum(as.vector(colSums(north_genpop@tab)[glob_lowfr])==0)
