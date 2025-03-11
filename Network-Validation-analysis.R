#make sure to install all the packages first
#install.packages(c('openxlsx','dplyr','epiR','irrCAC','igraph','lubridate','tidyr','purrr'))
library(igraph)
library(openxlsx)
library(dplyr)
library(epiR)
library(irrCAC)
library(lubridate)
library(tidyr)
library(purrr)


options(dplyr.summarise.inform = FALSE)

#for ease - just update this to where it is for you - then should just need to updated the file name to read in as well
setwd('H:/Code/prescriber networks/Ontario validation')
#read in the functions to create the network
source("Network-Validation-functions.R")

medWin <- openxlsx::read.xlsx('data/medication_windowsv20x.xlsx')

medWin$start_window <- as.Date(medWin$start_window, origin = "1899-12-30")
medWin$endwindow <- as.Date(medWin$endwindow, origin = "1899-12-30")

medWin$f_yr <- year(medWin$start_window)
medWin$l_yr <- year(medWin$endwindow)


docPat2 <- medWin %>%
  mutate(YR = map2(f_yr, l_yr, seq)) %>%
  unnest(YR) %>%
  select(YR, ptnt_id, cc_id, cp_id, windowcount)


colnames(docPat2) <- c('YR','PTNT_ID','CLINIC_ID','PHYSICIAN_ID','episode')

#confirmed previously 2021 was only up to Feb.
docPat <- docPat2[which(docPat2$YR %in% c(2013:2020)),]
docPat$com <- rep(0, nrow(docPat))
docPat$patload <- rep(0, nrow(docPat))
docPat$clust_coef <- rep(0, nrow(docPat))
docPat$strength <- rep(0, nrow(docPat))
docPat$comSize <- rep(0, nrow(docPat))
docPat$ties <- rep(0, nrow(docPat))

test_table <-
  docPat2 %>% group_by(YR) %>%
  summarize(n_pats=length(unique(PTNT_ID)), n_phys=length(unique(PHYSICIAN_ID)),
            n_clinic=length(unique(CLINIC_ID)))

full_data <- c()
clinic_data <- c()
pat_data <- c()
network_data <- c()
#number of iterations for the re-iteration method to run
#-5 is to call leading eigenvector - cluster_leading_eigen
#-4 is to call label propogation - cluster_label_prop
#-3 us to call infomap - cluster_infomap
#-2 is to call Leiden - cluster_leiden
#-1 is to call the walktrap algorithm instead - cluster_walktrap
#0 calls Louvain - cluster_louvain
# anything >= 1 is the iterative for cluster_fast_greedy
itr.max <- c(-4,-1, 0, 1, 2, 3, 4)
#how we define the ties in the network, absolute and relative
tie.type <- c('absolute', 'relative','episode')
#top% of ties to keep for the relative tie definition
rel.num <- c(.20, .40, .60, .80)
#number of shared patients required to define a tie in the absolute defition
abs.num <- c(1, 2, 5, 10, 20)
#how to define the main clinic for the physicians in the pairing
cln.alg <- c('max', 'any', 'percent')
#when we consider percent, what percent of time should be spent at the clinics?
#NOTE: anything > 50% would be the same as max.
cln.alg.perc <- c(0.45, 0.30, 0.15)

docPat_noep <- unique(docPat[,-5])

#loop through all the different settings to create the networks
#if we want to add other options add to the 

for (i in 1:length(itr.max)) {
  for (j in 1:length(cln.alg)) {
    for (k in 1:length(tie.type)) {
      #absolute and relative have different vector for numbers to loop through
      if (tie.type[k] == 'absolute') {
        for (l in 1:length(abs.num)) {
          if (cln.alg[j] == 'percent') {
            for (m in 1:length(cln.alg.perc)) {
              network <- create.network(docPat_noep,
                                        type = tie.type[k],
                                        num = abs.num[l],
                                        itr.max = itr.max[i],
                                        cln.alg = cln.alg[j],
                                        cln.alg.num = cln.alg.perc[m])
              full_data <- rbind(full_data, network$physicianPairs)
              clinic_data <- rbind(clinic_data, network$comLevel)
              pat_data <- rbind(pat_data, network$patLevel)
              network_data <- rbind(network_data, network$netInfo)
            }
          }
          else{
            network <- create.network(docPat_noep,
                                      type = tie.type[k],
                                      num = abs.num[l],
                                      itr.max = itr.max[i],
                                      cln.alg = cln.alg[j])
            full_data <- rbind(full_data, network$physicianPairs)
            clinic_data <- rbind(clinic_data, network$comLevel)
            pat_data <- rbind(pat_data, network$patLevel)
            network_data <- rbind(network_data, network$netInfo)
          }
        }
      }
      if (tie.type[k] == 'relative') {
        for (l in 1:length(rel.num)) {
          if (cln.alg[j] == 'percent') {
            for (m in 1:length(cln.alg.perc)) {
              network <- create.network(docPat_noep,
                                        type = tie.type[k],
                                        num = rel.num[l],
                                        itr.max = itr.max[i], 
                                        cln.alg = cln.alg[j],
                                        cln.alg.num = cln.alg.perc[m])
              full_data <- rbind(full_data, network$physicianPairs)
              clinic_data <- rbind(clinic_data, network$comLevel)
              pat_data <- rbind(pat_data, network$patLevel)
              network_data <- rbind(network_data, network$netInfo)
            }
          }
          else{
            network <- create.network(docPat_noep,
                                      type = tie.type[k],
                                      num = rel.num[l],
                                      itr.max = itr.max[i],
                                      cln.alg = cln.alg[j])
            full_data <- rbind(full_data, network$physicianPairs)
            clinic_data <- rbind(clinic_data, network$comLevel)
            pat_data <- rbind(pat_data, network$patLevel)
            network_data <- rbind(network_data, network$netInfo)
          }
        }
      }
      if (tie.type[k] == 'episode') {
          if (cln.alg[j] == 'percent') {
            for (m in 1:length(cln.alg.perc)) {
              network <- create.network(docPat,
                                        type = tie.type[k],
                                        num = NA,
                                        itr.max = itr.max[i], 
                                        cln.alg = cln.alg[j],
                                        cln.alg.num = cln.alg.perc[m])
              full_data <- rbind(full_data, network$physicianPairs)
              clinic_data <- rbind(clinic_data, network$comLevel)
              pat_data <- rbind(pat_data, network$patLevel)
              network_data <- rbind(network_data, network$netInfo)
            }
          }
          else{
            network <- create.network(docPat,
                                      type = tie.type[k],
                                      num = NA,
                                      itr.max = itr.max[i],
                                      cln.alg = cln.alg[j])
            full_data <- rbind(full_data, network$physicianPairs)
            clinic_data <- rbind(clinic_data, network$comLevel)
            pat_data <- rbind(pat_data, network$patLevel)
            network_data <- rbind(network_data, network$netInfo)
            
          }
      }
    }
  }
}


write.csv(clinic_data, file = "output/clinic_data_v6.csv", row.names = FALSE)
write.csv(network_data, file = "output/network_data_v6.csv", row.names = FALSE)


pat_data2 <- pat_data[which(pat_data$itr==3 & pat_data$cln_num==0.15 & pat_data$num==0.60),]

#obtain the number of physciains at the true clinic
docPat_clin <- unique(docPat_noep[,c(1:4)])
pat_data3 <- merge(x=pat_data2,y=docPat_clin,by.x=c('PTNT_ID','PHYSICIAN_ID','YR'),by.y=c('PTNT_ID','PHYSICIAN_ID','YR'),all.x=TRUE)

write.csv(pat_data3, file = "output/patient_data_v5.csv", row.names = FALSE)

#create the table for output by algorithm
table_dat <-
  full_data %>% group_by(type, num, alg_type,itr, cln_alg, cln_num) %>%
  summarize(n = n(),n_pos_real = sum(true_match), n_neg_real = sum(1 - true_match),
    n_agree = sum(concord), n_pos_agr = sum(pos_con), n_neg_agr = sum(neg_con),
    false_pos = sum((1 - true_match) * alg_match), false_neg = sum(true_match * (1 - alg_match)),
    overall_agr = round(n_agree / n, digits = 3),
    pos_agr = round((2 * n_pos_agr) / (n_pos_real + n_pos_agr + false_pos), digits = 3),
    neg_agr = round((2 * n_neg_agr) / (n_neg_real + n_neg_agr + false_neg), digits = 3),
    kappa_est = NA, kappa_ci = NA, gwet_est = NA,  gwet_ci = NA,
    sens_est = NA, sens_lb = NA, sens_ub = NA, spec_est = NA, spec_lb = NA, spec_ub = NA,
    ppv_est = NA, ppv_lb = NA, ppv_ub = NA, npv_est = NA, npv_lb = NA, npv_ub = NA
  )

#calculate the kappa, gwet, sensitivity, specificity, npv and ppv
for (i in 1:nrow(table_dat)) {
  table_row <- table_dat[i, ]
  xtab <- with(table_row, as.table(rbind(c(n_pos_agr, false_pos), c(false_neg, n_neg_agr))))
  
  kap <- kappa2.table(xtab)
  table_dat$kappa_est[i] <- round(kap[2], digits = 3)
  table_dat$kappa_ci[i] <- kap[4]
  
  gwet <- gwet.ac1.table(xtab)
  table_dat$gwet_est[i] <- round(gwet[2], digits = 3)
  table_dat$gwet_ci[i] <- gwet[4]

  epit <- epi.tests(xtab)
  table_dat$sens_est[i] <- round(epit$detail[3,]$est, digits = 3)
  table_dat$sens_lb[i] <- round(epit$detail[3,]$lower, digits = 3)
  table_dat$sens_ub[i] <- round(epit$detail[3,]$upper, digits = 3)
  
  table_dat$spec_est[i] <- round(epit$detail[4,]$est, digits = 3)
  table_dat$spec_lb[i] <- round(epit$detail[4,]$lower, digits = 3)
  table_dat$spec_ub[i] <- round(epit$detail[4,]$upper, digits = 3)
  
  table_dat$ppv_est[i] <- round(epit$detail[9,]$est, digits = 3)
  table_dat$ppv_lb[i] <- round(epit$detail[9,]$lower, digits = 3)
  table_dat$ppv_ub[i] <- round(epit$detail[9,]$upper, digits = 3)
  
  table_dat$npv_est[i] <- round(epit$detail[10,]$est, digits = 3)
  table_dat$npv_lb[i] <- round(epit$detail[10,]$lower, digits = 3)
  table_dat$npv_ub[i] <- round(epit$detail[10,]$upper, digits = 3)
}


#save the output
table_dat2 <- data.frame(lapply(table_dat, as.character), stringsAsFactors = FALSE)
write.csv(table_dat2, file = "output/concordance_stats_v7.csv", row.names = FALSE)


clinics <- c()


#loop through all the different settings to create the networks
for (j in 1:length(cln.alg)) {
  if (cln.alg[j] == 'percent') {
     for (m in 1:length(cln.alg.perc)) {
       clinic_raw <- extract.clinic(docPat_noep,cln.alg = cln.alg[j], cln.alg.num = cln.alg.perc[m])
       clinics <- rbind(clinics, clinic_raw)
      }
  }
  else{
    clinic_raw <- extract.clinic(docPat_noep,cln.alg = cln.alg[j])
    clinics <- rbind(clinics, clinic_raw)
  }
}

clinics_merge <- clinics[,c(1:4,9,10)]
#obtain the number of physciains at the true clinic
coms.wclin <- merge(x=clinic_data,y=clinics_merge,
                    by.x=c('main_clinic','YR','cln_alg','cln_num'),
                    by.y=c('CLINIC_ID','YR','cln_alg','cln_num'),
                    all.x=TRUE)

#calculated the recall and 1-purity
coms.wclin$recall <- coms.wclin$n_phys_main_clinic/coms.wclin$NUM_PHYSICIANS
coms.wclin$purity <- (coms.wclin$n_phys_main_clinic/coms.wclin$n_phys)
coms.wclin$minus_purity <- 1-coms.wclin$purity

#calculated F measure
coms.wclin$F_meas <- (2*coms.wclin$recall*coms.wclin$purity)/(coms.wclin$recall + coms.wclin$purity)

coms.wclin2 <- coms.wclin[-which(coms.wclin$cln_num %in% c(0.30,0.45)),]

#calculate the medians, 1st and 3rd quartiles for all measures
coms.sum <- coms.wclin2 %>% group_by(YR,type, num, alg_type,itr, cln_alg, cln_num) %>% 
  summarise(n_clins=length(unique(com)), n_phys_med=median(n_phys), n_phys_q1 = quantile(n_phys, probs=c(0.25)), n_phys_q3 = quantile(n_phys, probs=c(0.75)),
            n_pats_med=median(n_client), n_pats_q1 = quantile(n_client, probs=c(0.25)), n_pats_q3 = quantile(n_client, probs=c(0.75)),
            patload_med=median(patload), patload_q1 = quantile(patload, probs=c(0.25)), patload_q3 = quantile(patload, probs=c(0.75)),
            ties_med=median(ties), ties_q1 = quantile(ties, probs=c(0.25)), ties_q3 = quantile(ties, probs=c(0.75)),
            adj_str_med=median(adj_str), adj_str_q1 = quantile(adj_str, probs=c(0.25)), adj_str_q3 = quantile(adj_str, probs=c(0.75)),
            clust_coef_med=median(clust_coef), clust_coef_q1 = quantile(clust_coef, probs=c(0.25)), clust_coef_q3 = quantile(clust_coef, probs=c(0.75)),
            recall_mean=mean(recall), recall_sd=sd(recall),
            recall_med=median(recall), recall_q1 = quantile(recall, probs=c(0.25)), recall_q3 = quantile(recall, probs=c(0.75)),
            minus_p_mean=mean(minus_purity), minus_p_sd=sd(minus_purity),
            minus_purity_med=median(minus_purity), minus_purity_q1 = quantile(minus_purity, probs=c(0.25)), minus_purity_q3 = quantile(minus_purity, probs=c(0.75)),
            F_meas_mean=mean(F_meas), F_meas_sd=sd(F_meas),
            F_meas_med=median(F_meas), F_meas_q1 = quantile(F_meas, probs=c(0.25)), F_meas_q3 = quantile(F_meas, probs=c(0.75)))
            #apologize for the messy code here - but add any other variables output from coms.wclin such as num of rural physicians
            # - grab the median, 1st and 3rd quartile for all variables, then this will output data to use in a table for the paper


#output the table. I've kept the results for every network algorithm run - but we can select the few to present in Table2
coms.sum2 <- data.frame(lapply(coms.sum, as.character), stringsAsFactors = FALSE)
write.csv(coms.sum2, file = "output/network_community_stats_v3.csv", row.names = FALSE)

#create the network measures for the real clinic

#calculated the medians, 1st and 3rd quartile for the true clinics
clinic.sum <- clinics %>% group_by(YR, cln_alg, cln_num) %>%
  summarise(n_clins=length(unique(CLINIC_ID)), n_phys_med=median(NUM_PHYSICIANS), n_phys_q1 = quantile(NUM_PHYSICIANS, probs=c(0.25)), n_phys_q3 = quantile(NUM_PHYSICIANS, probs=c(0.75)),
            n_pats_med=median(NUM_CLIENTS), n_pats_q1 = quantile(NUM_CLIENTS, probs=c(0.25)), n_pats_q3 = quantile(NUM_CLIENTS, probs=c(0.75)),
            patload_med=median(patload), patload_q1 = quantile(patload, probs=c(0.25)), patload_q3 = quantile(patload, probs=c(0.75)),
            ties_med=median(ties), ties_q1 = quantile(ties, probs=c(0.25)), ties_q3 = quantile(ties, probs=c(0.75)),
            adj_str_med=median(adj_str), adj_str_q1 = quantile(adj_str, probs=c(0.25)), adj_str_q3 = quantile(adj_str, probs=c(0.75)),
            clust_coef_med=median(clust_coef), clust_coef_q1 = quantile(clust_coef, probs=c(0.25)), clust_coef_q3 = quantile(clust_coef, probs=c(0.75)))

#output the real clinic results
clinic.sum2 <- data.frame(lapply(clinic.sum, as.character), stringsAsFactors = FALSE)
write.csv(clinic.sum2, file = "output/network_clinic_stats_v3.csv", row.names = FALSE)


plot_com_prep <- pat_data3[which(pat_data3$YR==2020),]


#plot.sidebyside requires the known clinic ID and the equivalent community identifier
gs <- plot.sidebyside(plot_com_prep,149,19)
g.com<- gs$comGraph
com.weights <- log(E(g.com)$weight)+1
com.size <- log(V(g.com)$size)+1

g.clin <- gs$clinicGraph
clin.weights <- log(E(g.clin)$weight)+1
clin.size <- log(V(g.clin)$size)+1

phys <- unique(c(V(g.com)$name, V(g.clin)$name))
labels <- LETTERS[1:length(phys)]

#assign the labels to the unique physician ID for de-identification in visualization
com.labs <- c()
for(i in 1:length(V(g.com)$name)){
  let.lbl <- labels[which(phys==V(g.com)$name[i])]
  com.labs <- c(com.labs, let.lbl)
}

#repeat for clinics
clin.labs <- c()
for(i in 1:length(V(g.clin)$name)){
  let.lbl <- labels[which(phys==V(g.clin)$name[i])]
  clin.labs <- c(clin.labs, let.lbl)
}

par(mfrow=c(1,2))

#set the seed to ensure reproducibility for the network layout.
set.seed(2345)
plot(g.com,layout=layout_with_dh, main="Network Identified Clinic",
     edge.width=com.weights, vertex.size=com.size,
     vertex.label=com.labs, vertex.label.cex=2, vertex.label.dist=1.4,
     vertex.color='grey26', edge.color='lightgrey')
set.seed(2345)
plot(g.clin,layout=layout_with_dh, main = "Comparable True Clinic",
     edge.width=clin.weights, vertex.size=clin.size,
     vertex.label=clin.labs, vertex.label.cex=2, vertex.label.dist=1.4,
     vertex.color='grey26', edge.color='lightgrey')
