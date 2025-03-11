#function to recursively create communities within larger communities
#inputs:
#   j1 - counter to keep track of the amount of reiterations (starts at wave 0 - entire graph communities)
#   c1 - the current list of communities to split
#   comslist - the list of communities already found (starts as an empty list which is built onto once the communities are found)
#   pracmat - the working prescriber tie matrix
comm.recurs <- function(j1, c1, comslist, maxJ, pracmat){
  
  cmore <- communities(c1)[unname(which(sizes(c1) > 1))] #this selects all communities >1, we can change this to a different number (>4), and test it
  
  #recursive end check, are we past the third iteration? or are the communities all of size 1?
  if(j1==maxJ | length(cmore)==0){
    #j1==3 is the iterative step, do we stop after the 3rd iteration?
    comslist2 <- append(comslist, unname(communities(c1)))
    return(comslist2)
  }
  
  #if communities over size 1 and not at the 3rd wave, then continue
  else {
    #increase counter, we are at wave +1
    j.new <- j1+1
    
    for(i in 1:length(cmore)){
      
      #since we are using the numbers as names in later recursive attempts, they will be given as characters
      #need to turn them back to numeric in order to get the correct matrix to model
      testc <- as.numeric(cmore[[i]])
      matrem <- pracmat[testc,testc]
      
      #name the columns of the new matric so we can use them and keep track of the vertice numbers
      colnames(matrem) <- testc
      
      #create the new graph, the colnames is used from the matrix, which keeps the vertix numbers the same instead of renumbering them from 1
      g.new <- graph_from_adjacency_matrix(matrem,mode="undirected", weighted = TRUE, add.colnames = "name")
      
      #create the new communities from the new graph
      c.new <- cluster_fast_greedy(g.new, weights=E(g.new)$weight)
      
      #recursive call
      new.coms <- comm.recurs(j.new, c.new, comslist, maxJ,pracmat)
      
      #after the recursive call returns, add the recieved communities to our community list
      comslist <- unique(append(comslist, new.coms)) 
      #append isnt working like I thought it would, it's adding multiples, so unique removes the repeats
      #look for a better function that actually works like the .add() function in JAVA
      #this may be due to the recursive nature of the function/loop
    }
    cnone <- communities(c1)[unname(which(sizes(c1) <= 1))] 
    complus0 <- unique(append(comslist,cnone))
    return(complus0) 
  }
}

#function to create the network under the given tie definition
#inputs:
#   docPat - data to pass to the method, the list of doctor patient pairings
#   type - the type of tie definition, should accept only 'absolute' or 'relative'
#   num - number that sets the minimum required of patients or relative % based on the type
#   itr.max - the number of iterations to run under the community algorithm 
    #itr.max also determine which community detection algorithm will run.
    #-5 is to call leading eigenvector - cluster_leading_eigen
    #-4 is to call label propogation - cluster_label_prop
    #-3 us to call infomap - cluster_infomap
    #-2 is to call Leiden - cluster_leiden
    #-1 is to call the walktrap algorithm instead - cluster_walktrap
    #0 calls Louvain - cluster_louvain
    # anything >= 1 is the iterative for cluster_fast_greedy
#   cln.alg - the algorithm used to define which clinic to consider the physician attached too
#   cln.alg.num - only needed when cln.alg='percent', so defaults to 0, but defines the % of time spent at clinic
create.network <- function(docPat, type, num, itr.max, cln.alg, cln.alg.num=0){
 
require(igraph)
   
if(!tolower(type) %in% c('absolute','relative','episode')){stop('Type must be absolute, relative or episode')}

docPairs <- data.frame(pr1=numeric(), pr2=numeric(), year=numeric(), prac_shr=numeric(), 
                       pt_shr=numeric(), cln.match=logical(), com1=numeric(), com2=numeric())
  
years <- sort(unique(docPat$YR))
clinic_data <- c()
network_data <- c()

for(m in years){
  net <- c()
  
  y <- docPat[which(docPat$YR == m),]
  
  idnum <- length(unique(y$PTNT_ID))
  pracnum <- length(unique(y$PHYSICIAN_ID))
  clncnum <- length(unique(y$CLINIC_ID))
  
  y$num <- factor(y$PTNT_ID)
  #releveling for easy 1 to n identifiiers
  y$num <- as.numeric(y$num) 
  
  y$pr <- factor(y$PHYSICIAN_ID)
  #releveling for easy 1 to n identifiers
  y$pr <- as.numeric(y$pr)
  
  newmat <- matrix(data=0, nrow=pracnum, ncol=idnum)
  
  
  for(i in 1:nrow(y)){
    newmat[y$pr[i],y$num[i]] <- 1
  }
  
  
  pracmat <- matrix(data=0, nrow=pracnum, ncol=pracnum)
  
  pl <- rowSums(newmat) 
  
  for(i in 1:idnum){
    for(j in 1:(pracnum-1)){
      if(newmat[j,i]==1){
        for(k in (j+1):pracnum){
          if(newmat[k,i] == 1){
            if(tolower(type)=='episode'){
              j.ep <- unique(y$episode[which(y$pr==j & y$num==i)])
              k.ep <- unique(y$episode[which(y$pr==k & y$num==i)])
              ep.share <- sum(j.ep %in% k.ep)
              if(ep.share > 0){
                pracmat[j,k] <- pracmat[j,k]+1
                pracmat[k,j] <- pracmat[j,k]
              }
            }
            else{
            pracmat[j,k] <- pracmat[j,k]+1
            pracmat[k,j] <- pracmat[j,k]
            }
          }
        }
      }
    }
  }
  
  old_pracmat <- pracmat
  
  if(tolower(type)=='absolute'){
    #for ease of speed, do not run if the absolute number is 1
  if(num > 1) {
    for(j in 1:(pracnum-1)){
      for(k in (j+1):pracnum){
        #skip step if it's already 0
        if(pracmat[j,k] < num & pracmat[j,k] != 0){
          #set sharing to 0 if less then the absolute number for a tie definition
            pracmat[j,k] <- 0
            pracmat[k,j] <- 0
          }
        }
      }
    }
  }
  
  
  if(tolower(type)=='relative'){
    
  if(num > 1){stop('For relative, leave percentage in decimal form')}
  
    #only keeps shared ties with atleast 20% shared pts
    for(j in 1:(pracnum-1)){
      for(k in (j+1):pracnum){
        if(pracmat[j,k] > 0){
          #conditions is the perctile of 1-num (ex- top20% would be 80th percentile)
          percn <- 1-num
          nonzeroj <- pracmat[j,which(pracmat[j,] != 0)]
          nonzerok <- pracmat[k,which(pracmat[k,] != 0)]
          cond1 <- unname(quantile(nonzeroj, probs=percn))
          cond2 <- unname(quantile(nonzerok, probs=percn))
          if(pracmat[j,k] < cond1 & pracmat[j,k] < cond2){
            pracmat[j,k] <- 0
            pracmat[k,j] <- 0
          }
        }
      }
    }
  }
  
  
  g <- graph_from_adjacency_matrix(pracmat,mode="undirected", weighted = TRUE)
#  plot(g,layout=layout_with_fr, vertex.size=3,vertex.label=NA)
  net$year <- m
  net$edges <- ecount(g)
  net$loneVert <- sum(degree(g)==0)
  net$vertice <- vcount(g)
  net$trans <- transitivity(g)
  net$maxDeg <- max(degree(g))
  net$avgDeg <- mean(degree(g))

  
  #calculate strength before adjusting the weights
  str <- strength(g, vids = V(g), weights = E(g)$weights)
  
  
  #returns all communities (albiet in character, but these can be directly turned back to numeric)
  if(itr.max >= 1){
    c <- cluster_fast_greedy(g)
    all.coms <- comm.recurs(j1=0,c1=c,comslist=list(),itr.max,pracmat)
    com_alg <- 'Girvan-Newman'
  }
  else if(itr.max == 0){
    c <- cluster_louvain(g)
    all.coms <- communities(c)
    com_alg <- 'Louvain'
  }
  else if(itr.max == -1){
    c <- cluster_walktrap(g)
    all.coms <- communities(c)
    com_alg <- 'Walktrap'
  }
  else if(itr.max == -2){
    c <- cluster_leiden(g, objective_function = 'CPM')
    all.coms <- communities(c)
    com_alg <- 'Leiden'
  }
  else if(itr.max == -3){
    c <- cluster_infomap(g)
    all.coms <- communities(c)
    com_alg <- 'Infomap'
  }
  else if(itr.max == -4){
    c <- cluster_label_prop(g)
    all.coms <- communities(c)
    com_alg <- 'LabelPropogation'
  }
  else if(itr.max == -5){
    c <- cluster_leading_eigen(g)
    all.coms <- communities(c)
    com_alg <- 'LeadingEigenvector'
  }
  else if(itr.max == -6){
    c <- cluster_spinglass(g, spins=pracnum)
    all.coms <- communities(c)
    com_alg <- 'SpinGlass'
  }
  else{stop('Iteration maximum unrecognized for options')}
  
  cmore1 <- all.coms[which(lengths(all.coms) > 1)]
  
  cnone <- as.numeric(unlist(all.coms[which(lengths(all.coms) == 1)]))
  
  net$loneCom <- length(cnone)
  net$totalCom <- length(cmore1) + length(cnone)
  net$comalg <- com_alg
  net$clnalg <- cln.alg
  net$clnnum <- cln.alg.num
  net$edgetype <- type
  net$edgenum <- num
  net$itr <- itr.max
  
  clust_coef <- c() 
  
  if(length(cmore1) > 0){
  
  for(i in 1:length(cmore1)){
    testc <- as.numeric(cmore1[[i]])
    matrem <- pracmat[testc,testc]
    
    g.new <- graph_from_adjacency_matrix(matrem,mode="undirected", weighted = TRUE)
    
    clust_coef[i] <- transitivity(g.new, type = "global", isolates = "zero")
  }
  }
  
  y$str <- rep(0, nrow(y))
  y$comsize <- rep(0, nrow(y))
  y$com <- rep(0, nrow(y))
  y$clust_coef <- rep(0, nrow(y))
  y$ties <- rep(0, nrow(y))
  
  #assigning the ties, strength and community info to each physician patient pairing
  for(i in 1:pracnum){
    y$str[which(y$pr==i)] <- str[i]
    y$ties[which(y$pr==i)] <- length(which(old_pracmat[i,] > 0))
    if(length(cmore1) > 0){
    for(j in 1:length(cmore1)){
      if(i %in% cmore1[[j]]){
        y$com[which(y$pr==i)] <- j
        y$clust_coef[which(y$pr==i)] <- clust_coef[j]
        y$comsize[which(y$pr==i)] <- length(cmore1[[j]])
        
      }
    }
    }
    if(sum(y$com[which(y$pr==i)])==0){
      for(l in 1:length(cnone)){
        if(i == cnone[l]){
          y$com[which(y$pr==i)] <- j+l
          y$comsize[which(y$pr==i)] <- 1
        }
      }
    }
  }
  
  
  y$patld <- rowSums(newmat)[y$pr]
  docPat$patload[which(docPat$YR == m)] <- y$patld
  docPat$comSize[which(docPat$YR == m)] <- y$comsize
  docPat$com[which(docPat$YR == m)] <- y$com
  docPat$clust_coef[which(docPat$YR == m)] <- y$clust_coef
  docPat$strength[which(docPat$YR == m)] <- y$str
  docPat$ties[which(docPat$YR == m)] <- y$ties
  
  
#for each practitioner, we want to compare with the other practitioners and track whether they are at the same clinic
  #and if they were grouped into the same community network
for(j in 1:(pracnum-1)){
    #obtain all the clinics prescriber j works at
    cln.j.presc <-  y[which(y$pr==j),] %>% group_by(CLINIC_ID) %>% summarize(n_clients=length(unique(PTNT_ID)))
    total_cln.j <- sum(cln.j.presc$n_clients)
    cln.j.presc$perc_client <- cln.j.presc$n_clients/total_cln.j
    com.j <- unique(y$com[y$pr==j])
    #find the connected pracs - even though which returns just the positions, the positions are numbered their relabelled id, so it works
    connect.pracs.j <- which(old_pracmat[j,] != 0)
    
    for(k in (j+1):pracnum){
      #obtained all the clinic prescriber k works at
      cln.k.presc <-  y[which(y$pr==k),] %>% group_by(CLINIC_ID) %>% summarize(n_clients=length(unique(PTNT_ID)))
      total_cln.k <- sum(cln.k.presc$n_clients)
      cln.k.presc$perc_client <- cln.k.presc$n_clients/total_cln.k
      if(tolower(cln.alg)=='max'){
        #if max, take the clinic ID with the most prescriptions
        cln.k <- cln.k.presc$CLINIC_ID[which.max(cln.k.presc$n_clients)]
        cln.j <- cln.j.presc$CLINIC_ID[which.max(cln.j.presc$n_clients)]
        cln.match <- cln.j==cln.k
      }
      if(tolower(cln.alg)=='any'){
        #just check if there are any overlapping ones
        #max will take the value TRUE/1 if there are any that match
        any.in <- max(cln.k.presc$CLINIC_ID %in% cln.j.presc$CLINIC_ID)
        cln.match <- any.in==1
      }
      if(tolower(cln.alg)=='percent'){
        if(cln.alg.num > 1){stop('For top percentage of clinics, leave percentage in decimal form')}
        cln.j <- cln.j.presc$CLINIC_ID[which(cln.j.presc$perc_client >= cln.alg.num)]
        cln.k <- cln.k.presc$CLINIC_ID[which(cln.k.presc$perc_client >= cln.alg.num)]
        #if there are none over the given percentage, then just take the max
        if(length(cln.j)==0){cln.j <- cln.j.presc$CLINIC_ID[which.max(cln.j.presc$n_clients)]}
        if(length(cln.k)==0){cln.k <- cln.k.presc$CLINIC_ID[which.max(cln.k.presc$n_clients)]}
        #similiarly to any - if there is any overlap, max will take the value of TRUE/1
        any.in <- max(cln.k %in% cln.j)
        cln.match <- any.in==1
      }
      #save the community for k and their true patient sharing
      com.k <- unique(y$com[y$pr==k])
      true_pt_shr <- old_pracmat[j,k]
      connect.pracs.k <- which(old_pracmat[,k] != 0)
      #measure the number of commenly connected practitioners
      common_pracs <- length(intersect(connect.pracs.j, connect.pracs.k))
      newvec <- c(j,k,m,common_pracs,true_pt_shr,cln.match,com.j,com.k)
      docPairs[nrow(docPairs)+1,] <- newvec
      }
    }
  #true clincs and main clinics will be different based on the true clinic definition.
  #Remove rows in Y that do not fit with the defintion.
  
  if(tolower(cln.alg) != 'any'){
    
    for(j in 1:pracnum){
      #obtain all the clinics prescriber j works at
      cln.j.presc <-  y[which(y$pr==j),] %>% group_by(CLINIC_ID) %>% summarize(n_clients=length(unique(PTNT_ID)))
      total_cln.j <- sum(cln.j.presc$n_clients)
      cln.j.presc$perc_client <- cln.j.presc$n_clients/total_cln.j
      
      if(nrow(y) < pracnum){stop('Not enough rows to exist at practitioner level:',j)}
      
      if(tolower(cln.alg)=='percent'){
        if(cln.alg.num > 1){stop('For top percentage of clinics, leave percentage in decimal form')}
        cln.j.length <- length(cln.j.presc$CLINIC_ID[which(cln.j.presc$perc_client >= cln.alg.num)])
        
        #if there are none over the given percentage, then just take the max
        if(cln.j.length==0){
          rows.delete <- which(y$pr==j & y$CLINIC_ID != cln.j.presc$CLINIC_ID[which.max(cln.j.presc$perc_client)])
        }
        else{
          rows.delete <- which(y$pr==j & y$CLINIC_ID %in% cln.j.presc$CLINIC_ID[which(cln.j.presc$perc_client < cln.alg.num)])
        }
        if(!is_empty(rows.delete)){
          y <- y[-rows.delete,]
        }
      }
      
      if(tolower(cln.alg)=='max'){
        rows.delete <- which(y$pr==j & y$CLINIC_ID != cln.j.presc$CLINIC_ID[which.max(cln.j.presc$perc_client)])
        if(!is_empty(rows.delete)){
          y <- y[-rows.delete,]
        }
      }
    }
  }
  
  com.y <- y %>% group_by(YR, com) %>% summarise(n_phys=length(unique(PHYSICIAN_ID)), n_client=length(unique(PTNT_ID)),
                                                  main_clinic=as.numeric(names(sort(table(CLINIC_ID),decreasing=TRUE))[1]),
                                                  n_clinic=length(unique(CLINIC_ID)), 
                                                  n_phys_main_clinic=length(unique(PHYSICIAN_ID[which(CLINIC_ID==main_clinic)])))
  
  cur.yr <- docPairs[which(docPairs$year==m),]
  net$true_match <- sum(cur.yr$cln.match==1)
  net$alg_match <- sum(ifelse(cur.yr$com1==cur.yr$com2,1,0))
  net$concord <- sum(ifelse((cur.yr$cln.match==1 & cur.yr$com1==cur.yr$com2) | 
                             (cur.yr$cln.match==0 & cur.yr$com1!=cur.yr$com2),1,0))
  net$pos_con <- sum(ifelse((cur.yr$cln.match==1 & cur.yr$com1==cur.yr$com2),1,0))
  net$neg_con <- sum(ifelse((cur.yr$cln.match==0 & cur.yr$com1!=cur.yr$com2),1,0))
  clinic_data <- rbind(clinic_data, com.y)
  network_data <- rbind(network_data, net)
}

#clinics truly match
docPairs$true_match <- ifelse(docPairs$cln.match==TRUE,1,0)
docPairs$alg_match <- ifelse(docPairs$com1==docPairs$com2,1,0)
#either clinic and matched, or not clininc and matched
docPairs$concord <- ifelse((docPairs$cln.match==TRUE & docPairs$com1==docPairs$com2) | 
                             (docPairs$cln.match==FALSE & docPairs$com1!=docPairs$com2),1,0)
docPairs$pos_con <- ifelse((docPairs$cln.match==TRUE & docPairs$com1==docPairs$com2),1,0)
docPairs$neg_con <- ifelse((docPairs$cln.match==FALSE & docPairs$com1!=docPairs$com2),1,0)

docPairs$type <- type
docPairs$num <- num
docPairs$alg_type <- com_alg

docPairs$cln_alg <- cln.alg
docPairs$cln_num <- cln.alg.num

docPat$adj_str <- docPat$strength/docPat$patload

require(dplyr)


prac_info <- docPat %>% group_by(YR,PHYSICIAN_ID) %>% summarise(n_client=length(unique(PTNT_ID)), n_clinic=length(unique(CLINIC_ID)),
                                                                main_clinic=as.numeric(names(sort(table(CLINIC_ID),decreasing=TRUE))[1]),
                                                                com=max(com), patload=max(patload), clust_coef=max(clust_coef),
                                                                adj_str=max(adj_str), comSize=max(comSize), ties=max(ties))
                                                                  #rural=max(rural) - add this once rural indicator is added
                                                                        #assuming a physician is always either rural or urban and not both
#note all the inputs used in max are unique to the physician so only has one occurrence

prac_info$type <- type
prac_info$num <- num
prac_info$alg_type <- com_alg
prac_info$cln_alg <- cln.alg
prac_info$cln_num <- cln.alg.num

pat_info <- unique(docPat[,c('YR','com','PTNT_ID','PHYSICIAN_ID')])
#patnt_com should provide one row per individuals for year and network community - but an individual can see doctors in different communities
#clinic_pat <- patnt_com %>% group_by(YR, com) %>% summarise(n_cl=n())


clinic_phys <- prac_info %>% group_by(YR,com) %>% summarise(adj_str=median(adj_str),clust_coef=max(clust_coef), ties=median(ties), patload=median(patload))
#clust_coef is specific to the community, so the physicians in the community will all have the same - max will grab the correct number


clinic_info <- merge(x=clinic_data, y=clinic_phys, by=c('YR','com'))

clinic_info$type <- type
clinic_info$num <- num
clinic_info$alg_type <- com_alg
clinic_info$cln_alg <- cln.alg
clinic_info$cln_num <- cln.alg.num

pat_info$type <- type
pat_info$num <- num
pat_info$alg_type <- com_alg
pat_info$cln_alg <- cln.alg
pat_info$cln_num <- cln.alg.num

if(itr.max==-1){
  docPairs$itr <- 0
  clinic_info$itr <- 0
  pat_info$itr <- 0
}
else{
  docPairs$itr <- itr.max
  clinic_info$itr <- itr.max
  pat_info$itr <- itr.max
}



return(list(physicianPairs=docPairs, comLevel=clinic_info, patLevel=pat_info, netInfo=network_data))
}



#function to create the network variables for the true clinics
#inputs:
#   docPat - data to pass to the method, the list of doctor patient pairings
#   year - the year that we want to plot.
#   cln.alg - the algorithm used to define which clinic to consider the physician attached too
#   cln.alg.num - only needed when cln.alg='percent', so defaults to 0, but defines the % of time spent at clinic
extract.clinic <- function(docPat, cln.alg, cln.alg.num=0){

  require(dplyr)
  require(rlang)
  require(igraph)
  
  years <- sort(unique(docPat$YR))
  clin_out <- c()

  
  for(m in years){
  
  y <- docPat[which(docPat$YR == m),]
  
  idnum <- length(unique(y$PTNT_ID))
  pracnum <- length(unique(y$PHYSICIAN_ID))
  clncnum <- length(unique(y$CLINIC_ID))
  
  y$num <- factor(y$PTNT_ID)
  #releveling for easy 1 to n identifiiers
  y$num <- as.numeric(y$num) 
  
  y$pr <- factor(y$PHYSICIAN_ID)
  #releveling for easy 1 to n identifiers
  y$pr <- as.numeric(y$pr)
  
  
  if(tolower(cln.alg) != 'any'){
  
    for(j in 1:pracnum){
      #obtain all the clinics prescriber j works at
      cln.j.presc <-  y[which(y$pr==j),] %>% group_by(CLINIC_ID) %>% summarize(n_clients=length(unique(PTNT_ID)))
      total_cln.j <- sum(cln.j.presc$n_clients)
      cln.j.presc$perc_client <- cln.j.presc$n_clients/total_cln.j
    
      if(nrow(y) < pracnum){stop('Not enough rows to exist at practitioner level:',j)}
      
      if(tolower(cln.alg)=='percent'){
        if(cln.alg.num > 1){stop('For top percentage of clinics, leave percentage in decimal form')}
        
        rows.delete <- which(y$pr==j & y$CLINIC_ID %in% cln.j.presc$CLINIC_ID[which(cln.j.presc$perc_client < cln.alg.num)])
        if(!is_empty(rows.delete)){
          y <- y[-rows.delete,]
        }
      }
    
      if(tolower(cln.alg)=='max'){
        rows.delete <- which(y$pr==j & y$CLINIC_ID != cln.j.presc$CLINIC_ID[which.max(cln.j.presc$perc_client)])
        if(!is_empty(rows.delete)){
          y <- y[-rows.delete,]
        }
      }
    }
    
      
  }
  # clin.y <- y %>% group_by(YR, CLINIC_ID) %>%
  #   summarize(NUM_CLIENTS=length(unique(PTNT_ID)), NUM_PHYSICIANS=length(unique(PHYSICIAN_ID)))
  # 
  # clinics <- rbind(clinics, clin.y)
  # }
  # clinics$cln_alg <- cln.alg
  # clinics$cln_num <- cln.alg.num
  # 
  # return(clinics)
  
  newmat <- matrix(data = 0, nrow = pracnum, ncol = idnum)
  
  
  for (i in 1:nrow(y)) {
    newmat[y$pr[i], y$num[i]] <- 1
  }
  
  
  pracmat <- matrix(data = 0, nrow = pracnum, ncol = pracnum)
  
  for (i in 1:idnum) {
    for (j in 1:(pracnum - 1)) {
      if (newmat[j, i] == 1) {
        for (k in (j + 1):pracnum) {
          if (newmat[k, i] == 1) {
            pracmat[j, k] <- pracmat[j, k] + 1
            pracmat[k, j] <- pracmat[j, k]
          }
        }
      }
    }
  }
  
  g <- graph_from_adjacency_matrix(pracmat, mode = "undirected", weighted = TRUE)
  
  #calculate strength before adjusting the weights
  str <- strength(g, vids = V(g), weights = E(g)$weights)
  
  clinics <- unique(y$CLINIC_ID)
  clust_coef <- c()
  
  for (i in 1:length(clinics)) {
    cl_pr <- unique(y$pr[which(y$CLINIC_ID == clinics[i])])
    
    matrem <- pracmat[cl_pr, cl_pr]
    
    if (length(cl_pr) > 1) {
      g.new <- graph_from_adjacency_matrix(matrem, mode = "undirected", weighted = TRUE)
      
      clust_coef[i] <- transitivity(g.new, type = "global", isolates = "zero")
    }
    else{
      clust_coef[i] <- 0
    }
  }
  
  y$str <- rep(0, nrow(y))
  y$clust_coef <- rep(0, nrow(y))
  y$ties <- rep(0, nrow(y))
  
  
  for (i in 1:pracnum) {
    y$str[which(y$pr == i)] <- str[i]
    y$ties[which(y$pr == i)] <- length(which(pracmat[i,] > 0))
  }
  
  for (i in 1:length(clinics)) {
    y$clust_coef[which(y$CLINIC_ID == clinics[i])] <- clust_coef[i]
  }
  
  
  prac_info <-
    y %>% group_by(YR, PHYSICIAN_ID, CLINIC_ID) %>% 
    summarise(patload = length(unique(PTNT_ID)),
              str = max(str),
              ties = max(ties),
              adj_str=str/patload)
  
  
  #add in any other information about the clients here to summarize the clinic characteristics
  
  clin.y <- y %>% group_by(YR, CLINIC_ID) %>%
    summarize(NUM_CLIENTS = length(unique(PTNT_ID)),
              NUM_PHYSICIANS = length(unique(PHYSICIAN_ID)),
              clust_coef = max(clust_coef))
  
  clinic_phys <- prac_info %>% group_by(YR, CLINIC_ID) %>% 
    summarise(adj_str = median(adj_str),
              ties = median(ties),
              patload = median(patload))
  
  
  
  clinic_info <-
    merge(x = clin.y,
          y = clinic_phys,
          by = c('YR', 'CLINIC_ID'))
  
  clin_out <- rbind(clin_out, clinic_info)
  }

clin_out$cln_alg <- cln.alg
clin_out$cln_num <- cln.alg.num



return(clin_out)
  
}
  



#function to create the network variables for the true clinics
#inputs:
#   docPat - data to pass to the method, the list of doctor patient pairings with clinics and communities derived
#         This data also has the algorithm name and percentages used to define the true clinic pairing
#   clinic - clinicID to plot
#   comID - com to plot.

plot.sidebyside <- function(pat_data,clinic,comID){
  
  require(igraph)
  require(dplyr)
  require(rlang)
  
  type <- unique(pat_data$type)
  num <- unique(pat_data$num)
  itr.max <- unique(pat_data$itr)
  cln.alg <- unique(pat_data$cln_alg)
  cln.alg.num <- unique(pat_data$cln_num)
  
  y <- pat_data
  idnum <- length(unique(y$PTNT_ID))
  pracnum <- length(unique(y$PHYSICIAN_ID))
  clncnum <- length(unique(y$CLINIC_ID))
  
  pat_load <- y %>% group_by(PHYSICIAN_ID) %>% 
    summarize(n_pats=length(unique(PTNT_ID)), n_clins=length(unique(CLINIC_ID)))
  
  y$num <- factor(y$PTNT_ID)
  #releveling for easy 1 to n identifiiers
  y$num <- as.numeric(y$num) 
  
  y$pr <- factor(y$PHYSICIAN_ID)
  #releveling for easy 1 to n identifiers
  y$pr <- as.numeric(y$pr)
  
  newmat <- matrix(data=0, nrow=pracnum, ncol=idnum)
  
  for(i in 1:nrow(y)){
    newmat[y$pr[i],y$num[i]] <- 1
  }
  
  pracmat <- matrix(data=0, nrow=pracnum, ncol=pracnum)  
  
  for(i in 1:idnum){
    for(j in 1:(pracnum-1)){
      if(newmat[j,i]==1){
        for(k in (j+1):pracnum){
          if(newmat[k,i] == 1){
            if(tolower(type)=='episode'){
              j.ep <- unique(y$episode[which(y$pr==j & y$num==i)])
              k.ep <- unique(y$episode[which(y$pr==k & y$num==i)])
              ep.share <- sum(j.ep %in% k.ep)
              if(ep.share > 0){
                pracmat[j,k] <- pracmat[j,k]+1
                pracmat[k,j] <- pracmat[j,k]
              }
            }
            else{
              pracmat[j,k] <- pracmat[j,k]+1
              pracmat[k,j] <- pracmat[j,k]
            }
          }
        }
      }
    }
  }
  
  if(tolower(type)=='absolute'){
    #for ease of speed, do not run if the absolute number is 1
    if(num > 1) {
      for(j in 1:(pracnum-1)){
        for(k in (j+1):pracnum){
          #skip step if it's already 0
          if(pracmat[j,k] < num & pracmat[j,k] != 0){
            #set sharing to 0 if less then the absolute number for a tie definition
            pracmat[j,k] <- 0
            pracmat[k,j] <- 0
          }
        }
      }
    }
  }
  
  
  if(tolower(type)=='relative'){
    
    if(num > 1){stop('For relative, leave percentage in decimal form')}
    
    #only keeps shared ties with atleast 20% shared pts
    for(j in 1:(pracnum-1)){
      for(k in (j+1):pracnum){
        if(pracmat[j,k] > 0){
          #conditions is the perctile of 1-num (ex- top20% would be 80th percentile)
          percn <- 1-num
          nonzeroj <- pracmat[j,which(pracmat[j,] != 0)]
          nonzerok <- pracmat[k,which(pracmat[k,] != 0)]
          cond1 <- unname(quantile(nonzeroj, probs=percn))
          cond2 <- unname(quantile(nonzerok, probs=percn))
          if(pracmat[j,k] < cond1 & pracmat[j,k] < cond2){
            pracmat[j,k] <- 0
            pracmat[k,j] <- 0
          }
        }
      }
    }
  }
  
  y1 <- y[which(y$com==comID),]
  #create subset matrix just for community and plot prior to cleaning data for the clinic
  pracmat.n <- pracmat[unique(y1$pr),unique(y1$pr)]
  g.net <- graph_from_adjacency_matrix(pracmat.n,mode="undirected", weighted = TRUE)
  V(g.net)$name <- unique(y1$PHYSICIAN_ID)
  
  pat1 <- pat_load$n_pats[unique(y1$pr)]
  V(g.net)$size <- pat1
  
  if(tolower(cln.alg) != 'any'){
    
    for(j in 1:pracnum){
      #obtain all the clinics prescriber j works at
      cln.j.presc <-  y[which(y$pr==j),] %>% group_by(CLINIC_ID) %>% summarize(n_clients=length(unique(PTNT_ID)))
      total_cln.j <- sum(cln.j.presc$n_clients)
      cln.j.presc$perc_client <- cln.j.presc$n_clients/total_cln.j
      
      if(nrow(y) < pracnum){stop('Not enough rows to exist at practitioner level:',j)}
      
      if(tolower(cln.alg)=='percent'){
        if(cln.alg.num > 1){stop('For top percentage of clinics, leave percentage in decimal form')}
        cln.j.length <- length(cln.j.presc$CLINIC_ID[which(cln.j.presc$perc_client >= cln.alg.num)])
        
        #if there are none over the given percentage, then just take the max
        if(cln.j.length==0){
          rows.delete <- which(y$pr==j & y$CLINIC_ID != cln.j.presc$CLINIC_ID[which.max(cln.j.presc$perc_client)])
        }
        else{
          rows.delete <- which(y$pr==j & y$CLINIC_ID %in% cln.j.presc$CLINIC_ID[which(cln.j.presc$perc_client < cln.alg.num)])
        }
        if(!is_empty(rows.delete)){
          y <- y[-rows.delete,]
        }
      }
      
      if(tolower(cln.alg)=='max'){
        rows.delete <- which(y$pr==j & y$CLINIC_ID != cln.j.presc$CLINIC_ID[which.max(cln.j.presc$perc_client)])
        if(!is_empty(rows.delete)){
          y <- y[-rows.delete,]
        }
      }
    }
  }
  
  
  
  
  y2 <- y[which(y$CLINIC_ID==clinic),]
  #create subset matrix just for communit7
  pracmat.c <- pracmat[unique(y2$pr),unique(y2$pr)]
  
  g.clin <- graph_from_adjacency_matrix(pracmat.c,mode="undirected", weighted = TRUE)
  V(g.clin)$name <- unique(y2$PHYSICIAN_ID)
  pat2 <- pat_load$n_pats[unique(y2$pr)]
  V(g.clin)$size <- pat2
  
  #add to the list physicians to add the physician tag onto.
  return(list(comGraph=g.net, clinicGraph=g.clin))
}

