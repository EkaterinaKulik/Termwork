library(sna)




bills_sum <- read.table('/Users/ekaterinakulik/Downloads/co_sponsorship_matrix-4.csv', header = TRUE, sep = ',', row.names = 1)
bills_net <- read.table('/Users/ekaterinakulik/Downloads/Matrix_for_termwork/binary_co_sponsorship_matrix-3.csv', header = TRUE, sep = ',', row.names = 1)
sex_net <- read.table('/Users/ekaterinakulik/Downloads/Matrix_for_termwork/co_sex_matrix-2.csv', header = TRUE, sep = ',', row.names = 1)
age_net <- read.table('/Users/ekaterinakulik/Downloads/Matrix_for_termwork/co_age_matrix-2.csv', header = TRUE, sep = ',', row.names = 1)
Edu_net <- read.table('/Users/ekaterinakulik/Downloads/co_education_matrix-3.csv', header = TRUE, sep = ',', row.names = 1)
Spec_net <- read.table('/Users/ekaterinakulik/Downloads/co_specialisation_matrix-3.csv', header = TRUE, sep = ',', row.names = 1)
regions_net <- read.table('/Users/ekaterinakulik/Downloads/Matrix_for_termwork/co_occurrence_matrix-2.csv', header = TRUE, sep = ',', row.names = 1)
faction_UR_net <- read.table('/Users/ekaterinakulik/Downloads/intersection_matrix_Единая Россия-3.csv', header = TRUE, sep = ',', row.names = 1)
faction_KPRF_net <- read.table('/Users/ekaterinakulik/Downloads/intersection_matrix_КПРФ-3.csv', header = TRUE, sep = ',', row.names = 1)
faction_LDPR_net <- read.table('/Users/ekaterinakulik/Downloads/intersection_matrix_ЛДПР-3.csv', header = TRUE, sep = ',', row.names = 1)
faction_SR_net <- read.table('/Users/ekaterinakulik/Downloads/intersection_matrix_Справедливая Россия-3.csv', header = TRUE, sep = ',', row.names = 1)
Committee_net <- read.table('/Users/ekaterinakulik/Downloads/co_committees_matrix-3.csv', header = TRUE, sep = ',', row.names = 1)
work_net <- read.table('/Users/ekaterinakulik/Downloads/Matrix_for_termwork/co_work_matrix-2.csv', header = TRUE, sep = ',', row.names = 1)
coalition_net <- read.table('/Users/ekaterinakulik/Downloads/Matrix_for_termwork/co_coalitions_matrix.csv', header = TRUE, sep = ',', row.names = 1)
coalition_net_4 <- read.table('/Users/ekaterinakulik/Downloads/co_coalitions_matrix_for_4_clusters.csv', header = TRUE, sep = ',', row.names = 1)
if ("X" %in% names(coalition_net)) {
  rownames(coalition_net) <- coalition_net$X
  coalition_net <- coalition_net[ , -which(names(coalition_net) == "X")]
}
exists("QAP_df") 

nl<-netlogit(bills_net,           # Dependent variable/network
          list(sex_net, age_net, Edu_net, Spec_net, regions_net, 
               faction_UR_net, faction_KPRF_net, faction_LDPR_net, 
               faction_SR_net, Committee_net, work_net, coalition_net), # List the independent variables/networks
          reps=100)     
print("QAP Analysis complete")
summary(nl)

nl_4<-netlogit(bills_net,           # Dependent variable/network
             list(sex_net, age_net, Edu_net, Spec_net, regions_net, 
                  faction_UR_net, faction_KPRF_net, faction_LDPR_net, 
                  faction_SR_net, Committee_net, work_net, coalition_net_4), # List the independent variables/networks
             reps=100)     
print("QAP Analysis complete")
summary(nl_4)
df_4 <- data.frame(Betas=nl_4$coefficients, P.values=nl_4$pgreqabs, tstat=nl_4$tstat)
df_4$index <- seq(from = nrow(df_4), to = 1, by = -1)
df_4$Names <- c(`1` = "Константа (edges/density)",  `2` ="Общий пол", `3` ="Разница в возрасте", `4` ="Общее образование", `5` ="Общая специальность",`6` ="Общий регион", `7` ="Фракция Единая Россия",`8` ="Фракция КПРФ", `9` ="Фракция ЛДПР", `10` ="Фракция СР", `11` ="Общий комитет", `12` = "Общее место предыдущей работы", `13` = "Общая коалиция") 
print("QAP Analysis complete")

df_total_4 <- df_4
df_total_4$P.values<- round(df_total_4$P.values, digits = 4)
df_total_4$SEs<- 0

write.csv(df_total_4, file = "model_data_for_4_clusters.csv", row.names = FALSE)


df <- data.frame(Betas=nl$coefficients, P.values=nl$pgreqabs, tstat=nl$tstat)
df$index <- seq(from = nrow(df), to = 1, by = -1)
df$Names <- c(`1` = "Константа (edges/density)",  `2` ="Общий пол", `3` ="Разница в возрасте", `4` ="Общее образование", `5` ="Общая специальность",`6` ="Общий регион", `7` ="Фракция Единая Россия",`8` ="Фракция КПРФ", `9` ="Фракция ЛДПР", `10` ="Фракция СР", `11` ="Общий комитет", `12` = "Общее место предыдущей работы", `13` = "Общая коалиция") 


print("QAP Analysis complete")

df_total <- df
df_total$P.values<- round(df_total$P.values, digits = 4)
df_total$SEs<- 0

write.csv(df_total, file = "model_data.csv", row.names = FALSE)


g2 <- print(ggforestplot::forestplot(
  df = df_total,
  name = Names,
  estimate = Betas,
  se = SEs,
  pvalue = P.values,
  psignif = 0.1,
  title = paste("Результаты QAP модели", sep=" "),
  colour = model
)+labs(
  title = "Результаты QAP модели",
  subtitle = "VII созыв",
)+ theme(
  title=element_text(size=12),
  axis.text=element_text(size=10),
  axis.title=element_text(size=10,face="bold"),
  plot.title.position = "plot", 
  plot.caption.position =  "plot"))
ggsave(paste("Результаты QAP модели", sep=""), width = 20, height = 16, units = c("cm"), dpi = 300 )



install.packages("ergm")
install.packages("robustbase")

library(ergm)
library(network)
row.names(bills_net)<-seq(1,length(bills_net))
colnames(bills_net)<-seq(1,length(bills_net))



clean_adjacency_matrix <- function(adj_matrix) {
  diag(adj_matrix) <- 0  
  adj_matrix[adj_matrix > 1] <- 1  
  return(adj_matrix)
}


create_network <- function(matrix) {
  network(matrix, directed = TRUE, loops = FALSE, multiple = FALSE)
}


bills_net_n_clean <- create_network(clean_adjacency_matrix(as.matrix(bills_net)))
sex_net_n_clean <- create_network(clean_adjacency_matrix(as.matrix(sex_net)))
age_net_n_clean <- create_network(clean_adjacency_matrix(as.matrix(age_net)))
Edu_net_n_clean <- create_network(clean_adjacency_matrix(as.matrix(Edu_net)))
Spec_net_n_clean <- create_network(clean_adjacency_matrix(as.matrix(Spec_net)))
regions_net_n_clean <- create_network(clean_adjacency_matrix(as.matrix(regions_net)))
faction_UR_net_n_clean <- create_network(clean_adjacency_matrix(as.matrix(faction_UR_net)))
faction_KPRF_net_n_clean <- create_network(clean_adjacency_matrix(as.matrix(faction_KPRF_net)))
faction_LDPR_net_n_clean <- create_network(clean_adjacency_matrix(as.matrix(faction_LDPR_net)))
faction_SR_net_n_clean <- create_network(clean_adjacency_matrix(as.matrix(faction_SR_net)))
Committee_net_n_clean <- create_network(clean_adjacency_matrix(as.matrix(Committee_net)))
work_net_n_clean <- create_network(clean_adjacency_matrix(as.matrix(work_net)))
coalition_net_n_clean <- create_network(clean_adjacency_matrix(as.matrix(coalition_net)))
coalition_net_n_clean_4 <- create_network(clean_adjacency_matrix(as.matrix(coalition_net_4)))



model_2.1 <- ergm(bills_net_n_clean ~ edges+
                    dgwesp(0, fixed=T, type='ISP')+
                    edgecov(sex_net_n_clean)+ 
                    edgecov(age_net_n_clean)+ 
                    edgecov(Edu_net_n_clean)+ 
                    edgecov(Spec_net_n_clean)+ 
                    edgecov(regions_net_n_clean)+ 
                    edgecov(faction_UR_net_n_clean)+ 
                    edgecov(faction_KPRF_net_n_clean)+ 
                    edgecov(faction_LDPR_net_n_clean)+ 
                    edgecov(faction_SR_net_n_clean)+ 
                    edgecov(Committee_net_n_clean)+
                    edgecov(work_net_n_clean)+
                    edgecov(coalition_net_n_clean),
                  control=control.ergm(parallel=8, parallel.type="PSOCK")
)

model4_2.1 <- ergm(bills_net_n_clean ~ edges+
                    dgwesp(0, fixed=T, type='ISP')+
                    edgecov(sex_net_n_clean)+ 
                    edgecov(age_net_n_clean)+ 
                    edgecov(Edu_net_n_clean)+ 
                    edgecov(Spec_net_n_clean)+ 
                    edgecov(regions_net_n_clean)+ 
                    edgecov(faction_UR_net_n_clean)+ 
                    edgecov(faction_KPRF_net_n_clean)+ 
                    edgecov(faction_LDPR_net_n_clean)+ 
                    edgecov(faction_SR_net_n_clean)+ 
                    edgecov(Committee_net_n_clean)+
                    edgecov(work_net_n_clean)+
                    edgecov(coalition_net_n_clean_4),
                  control=control.ergm(parallel=8, parallel.type="PSOCK")
)

library(network)
df_M_1.1 <- summary(model_2.1)
df_M_1.1<-as.data.frame(df_M_1.1$coefficients)
df_ergm_1 <- data.frame(Betas=df_M_1.1$Estimate, SEs=df_M_1.1$`Std. Error`, P.values=df_M_1.1$`Pr(>|z|)`)
df_ergm_1$P.values[is.na(df_ergm_1$P.values)]<-1
df_ergm_1$P.values<- round(df_ergm_1$P.values, digits = 4)
df_ergm_1$index <- seq(from = nrow(df_ergm_1), to = 1, by = -1)
df_ergm_1$CI_low <- df_ergm_1$Betas - (qnorm(0.975)*df_ergm_1$SEs)
df_ergm_1$CI_high <- df_ergm_1$Betas + (qnorm(0.975)*df_ergm_1$SEs)
df_ergm_1[is.na(df_ergm_1)]<-0

p.cols <- vector()
for (l in 1:nrow(df_ergm_1)) {
  if (df_ergm_1$CI_low[l]  < 0 & df_ergm_1$CI_high[l] > 0) {
    p.cols[l] <- "#d7191c" # Color for insignificant estimates 
  } else {
    p.cols[l] <- "#2c7bb6" # Color for significant estimates 
    if (df_ergm_1$CI_low[l] == 0) {
      p.cols[l] <- "#d7191c" # Color for insignificant estimates 
    }
  }
}; rm(l)

df_ergm_1$model <- paste(i, " Созыв", sep= "")

#### добавить названия переменных
df_ergm_1$Names <- c(`1` = "Константа (edges/density)",`2` ="gwesp_0", `3` ="Общий пол", `4` ="Разница в возрасте", `5` ="Общее образование", `6` ="Общая специальность",`7` ="Общий регион", `8` ="Фракция Единая Россия",`9` ="Фракция КПРФ", `10` ="Фракция ЛДПР", `11` ="Фракция СР", `12` ="Общий комитет", `13` = "Общее место предыдущей работы", `14` = "Общая коалиция") 
rownames(df_ergm_1)<-df_ergm_1$index
write.csv(df_ergm_1, file = "df_ergm_1.csv", row.names = FALSE)

df_ergm_1_array<-append(df_ergm_1_array, df_ergm_1)


g3 <- print(ggforestplot::forestplot(
  df = df_total1,
  name = Names,
  estimate = Betas,
  se = SEs,
  pvalue = P.values,
  psignif = 0.05,
  title = paste("Результаты ERGM модели", sep=" "),
  colour = model
)+labs(
  title = "Результаты ERGM модели",
  subtitle = "VII созыв Государственной Лумы РФ",
)+ theme(
  title=element_text(size=12),
  axis.text=element_text(size=10),
  axis.title=element_text(size=10,face="bold"),
  plot.title.position = "plot", 
  plot.caption.position =  "plot"))
ggsave(paste("Результаты ERGM модели", sep=""), width = 20, height = 16, units = c("cm"), dpi = 300 )

