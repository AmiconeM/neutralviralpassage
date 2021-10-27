#---
#  title: "Virus replication and migration"
#author: "Amicone Massimo"
#date: "12/04/2021"
#---

U=0.1 ##mutation rate
N=2000000
bottleneck=1/1000 ##dilution
growth=1/bottleneck ## average growth factor
L=30000 ##number of sites
migration=1/10 ## migrating fraction

## input output file names
input_file=paste("Virus_Evolution_Pool_N",N,"_L",L,".RData",sep="") ## to store pools of genotypes to be used in the migration simulations
output_file=paste("Virus_Evolution_M",migration,"_N",N,"_L",L,".RData",sep="")

generations=15 ##number of passages
simulations=100

Freq_all=(1:generations)*U
Counts_all=c()
Mutation_freq_all=list()
N_all=rep(N,generations)

N_pool=N*bottleneck*migration

### load pools of genomes (list over time)
load(file = input_file)

for(pop in 1:simulations){
  print(paste("Population",pop))
  M=1 #number of lineages/mutations (SGV)
  
  
  id=1:M ## label for each lineage
  
  Genotypes=matrix(0,nrow=M,ncol=L)
  
  frequencies=rep(1/M,M) ##sgv frequencies
  n=N*bottleneck*frequencies ##abundances
  
  
  freq=c()
  N_t=c()
  for(t in 1:generations){
    ##growth (expand by x)
    N_new=n*rpois(length(n),growth)
    N_t=c(N_t,sum(N_new))
    #print(sum(N_new))
    
    ## mutations
    parents=c()
    mut_id=c()
    M_t=M
    N_old=N_new
    tot_mutations=0
    
    for(m in 1:M){
      n_i=N_new[m]
      mutations= min(n_i,rpois(1,U*n_i))
      #mutations= U*n_i
      if(mutations>0){
        parents=c(parents,rep(m,mutations))
        mut_id=c(mut_id,((M_t+1):(M_t+mutations)))
        N_new[m]=n_i-mutations
        N_new=c(N_new,rep(1,mutations)) ##mutants strat with freq 1/N
        M_t=M_t+mutations
        
      }
      tot_mutations=tot_mutations+mutations
    }
    id=1:M_t
    
    #if(t==generations){
    Counts=c()
    for(j in 1:L){
      Counts=c(Counts,sum(Genotypes[,j]*N_old))
    }
    
    Counts=Counts[which(Counts>0)]
    freq=c(freq,(sum(Counts)+(M_t-M))/(sum(N_new)))
    #}
    
    ##sample genotype for next cycle  (bottleneck)
    N_sample=round(sum(N_new)*bottleneck,digits = 0) ## absolute number of individuals to pick
    sampled_pop=sample(id,N_sample,prob=N_new,replace = T) ## IDs of the sampled genotypes
    id_new=sort(unique(sampled_pop)) #unique IDs
    
    ##compute new abundances
    n=c()##new abundances
    for(x in id_new){
      n=c(n,length(which(sampled_pop==x)))
    }
    
    ##update mutant genotypes
    new_mutant=which(id_new>M) ##retrieve mutants positions that survive the sampling
        
    new_mut_id=id_new[which(id_new>M)]
    surv=setdiff(id_new,new_mut_id)
    
    if(length(new_mut_id)>0){
      G_new=matrix(0,nrow=1,ncol=L)
      for(i in 1:length(new_mut_id)){
        parent=parents[which(mut_id==new_mut_id[i])]
        G=Genotypes[parent,]
        site=sample(1:L,1)
        G[site]=(1+G[site])%%2 ##allow back mutations
        G_new=rbind(G_new,G)
      }
    }
    
    ##remove dead genotypes
    if(length(surv)==1){
      Genotypes=t(as.matrix(Genotypes[surv,]))
    }else{Genotypes=Genotypes[surv,]}
    
    
    ##add mutations
    if(length(new_mut_id)>0){
      Genotypes=rbind(Genotypes,G_new[2:(length(new_mut_id)+1),])
    }
    
    M=length(n)
    
    #### Migration bit ############
    ## move out/kill random geotypes to be replaced by migration
    id=1:M
    tokill=sample(rep(id,n),N_pool,replace = F) ## sample N_pool individuals
    id_kill=sort(unique(tokill)) ## get their IDs
    
    ##update abundances
    for(i in 1:length(id_kill)){
      j=id_kill[i]
      n[j]=n[j]-min(n[j],length(which(tokill==j)))
    }
    
    ## remove extinct
    surv=which(n>0)
    Genotypes=Genotypes[surv,]
    n=n[surv]
    M=length(n)
    
    ##import from pool
    n_pool=Pool_n[[t]]
    id=1:length(n_pool)
    pool=sample(rep(id,n_pool),N_pool,replace = F) #sample N_pool from the pool
    id_pool=sort(unique(pool))
    
    ##update abundances
    
    for(i in 1:length(id_pool)){
      n=c(n,length(which(pool==id_pool[i])))
    }
    ##update Genotypes
    Genotypes=rbind(Genotypes,Pool_G[[t]][id_pool,])
    M=length(n)
    
  }
  Counts_all=c(Counts_all,Counts)  ##store counts only at last generation
  Freq_all=rbind(Freq_all,freq)
  N_all=rbind(N_all,N_t)
  Mutation_freq_all[[pop]]=c(Counts,rep(1,tot_mutations))/N_t[generations]
  
}

# Save multiple objects
save(Counts_all,Freq_all,N_all,Mutation_freq_all, file = output_file)

#```
