# Set root and main directories
root <- "/mmfs1/home/david.adeleke/Thesis"
main <- "SNV"
proj.scope <- c("PAAD-US", "PACA-CA", "PACA-AU")
# Create directory if it doesn't exist and set working directory
dir.create(file.path(root, main), showWarnings = FALSE, recursive = TRUE)



setwd(file.path(root, main))

# Load necessary libraries
required_packages <- c(
  "circlize", "tidyverse", "ComplexHeatmap", "RTCGA.clinical",
  "SummarizedExperiment", "TCGAbiolinks", "ggplot2", "DT", "tidyr",
  "reshape2", "plyr", "plotly", "tibble", "hrbrthemes", "viridis",
  "GGally", "kableExtra", "survival", "survminer", "grid", "gridExtra",
  "rmdformats","R.utils","readr","data.table")


# Load libraries, install if not available
lapply(required_packages, function(pkg) {
  if (!require(pkg, character.only = TRUE)) {
    install.packages(pkg, dependencies = TRUE)
    library(pkg, character.only = TRUE)
  }
})

message(paste0("Working directory is set successfully at ", getwd()))

# Initialize junk list
junck.list <- c()


# Function to convert SNV data to matrix

source("alphamisScript.R")



maxl=1
for(i in 1:maxl){
  coho=i
# Initialize the cohort


# Load SNV and clinical data
snv_files <- c("dataset/PAAD-US_ICGC_snv.rds", "dataset/PACA-CA_ICGC_snv.rds", "dataset/PACA-AU_ICGC_snv.rds")

clinical_files <- c("dataset/PAAD-US_ICGC_clin.rds", "dataset/PACA-CA_ICGC_clin.rds", "dataset/PACA-AU_ICGC_clin.rds")

snv <- readRDS(file.path(root, snv_files[coho]))
snv <- snv %>% dplyr::filter(Hugo_Symbol != "UnknownGene")
#snv_k= snv %>% dplyr::filter(Hugo_Symbol == "KRAS", cds_mutation %in% c("34G>C", "35G>A","35G>T" ))
clinical <- readRDS(file.path(root, clinical_files[coho]))

# Process SNV data
snv_mat <- snv2mat(snv = snv, i = coho)

snv_path <- snv2path(snv = snv, i = coho)

# Create tumor ID table
clin=make_tumorid(output_file=".rds", i = coho, cutoff=0.3)

# Relevel factors
ref <- "KpTp"
clin$grp <- as.factor(clin$grp)
clin <- within(clin, grp <- relevel(grp, ref = ref))


clin$grpbin <- as.factor(clin$grpbin)
clin <- within(clin, grpbin <- relevel(grpbin, ref = ref))

#"mutClinPAAD-US.rds"
# coxph(Surv(time, status) ~ grpbin, data = clin[clin$time >=10,]) %>% 
#   gtsummary::tbl_regression(exp = TRUE) 
# 
# coxph(Surv(time, status) ~ grp, data = clin[clin$time >=10,]) %>% 
#   gtsummary::tbl_regression(exp = TRUE)

# proj.scope <- c("PAAD-US", "PACA-CA", "PACA-AU")
# ppt= create_upset_plot(snv_path)
# png(file=paste0("images/", proj.scope[coho],"upset.png"))
# ppt

# htmp=createHeatmap(snv_path, clin, 25)
# png(file=paste0("images/", proj.scope[coho],"heatmp.png"))
# htmp

}

misalpha$id = paste0(misalpha$CHROM,misalpha$POS, misalpha$REF, misalpha$ALT)


subalpha=misalpha[,c(9:11)]

subalpha$asin_am <- asin(subalpha$am_pathogenicity)



snv$id= paste0("chr",snv$Chromosome,snv$Start_Position,snv$Reference_Allele, snv$Tumor_Seq_Allele2)

c.snp=intersect(snv$id, subalpha$id)

snv1 = snv[snv$id %in% c.snp,]
subalpha = subalpha[subalpha$id %in% c.snp,]
mis.snv=merge(snv1, subalpha, by ="id")

mis.snv <- mis.snv %>%
  select(icgc_donor_id, am_pathogenicity) %>%
  group_by(icgc_donor_id) %>%
  summarise(
    tmb = n(),
    am_pathogenicity = sum(am_pathogenicity, na.rm = TRUE),
    .groups = "drop"
  )

clinical <- clinical %>% filter(project_code %in% proj.scope[i])



clin.mis.snv= merge(clinical, mis.snv, by="icgc_donor_id")

clin.mis.snv=clin.mis.snv[,c(1,22:25)]



# Define the list of variables
class_grp1 <- c("tmb","am_pathogenicity")

# Iterate over each variable and apply transformation



clin.all=readRDS(( "/mmfs1/home/david.adeleke/dissertation/data/processed_data/tcga_us/clin.all.rds"))


clin.all <- clin.all %>%
  mutate(across(
    all_of(class_grp1),
    ~ case_when(
      . >= quantile(., 0.50, na.rm = TRUE) ~ "HighPath",
      . >= quantile(., 0.15, na.rm = TRUE) & . < quantile(., 0.50, na.rm = TRUE) ~ "MildPath",
      TRUE ~ "LowPath"
    ),
    .names = "class_{.col}"
  ))



clin.all$Patho=clin.all$class_am_pathogenicity

clin.all$Patho <- factor(clin.all$Patho, levels = c("LowPath", setdiff(unique(clin.all$Patho), "LowPath")))

sfit <- survfit(Surv(time, Overall_Survival_Status) ~ Patho, data = clin.all)

# Plot survival curves
ggsurvplot(
  sfit,
  risk.table = TRUE,
  pval = TRUE,
  conf.int = FALSE,
  #palette = c("#47BC82", "#FAA928", "#B95D4B", "#A18CBA"),
  xlim = c(0, 2500),
  break.time.by = 200,
  #ggtheme = theme_RTCGA(),
  risk.table.y.text.col = TRUE,
  risk.table.y.text = FALSE,
  risk.table.height = 0.4,
  legend.title = "Percentiles")


table(clin.all$Group, clin.all$amp)

saveRDS(clin.all, "/mmfs1/home/david.adeleke/dissertation/data/processed_data/tcga_us/clin.all.rds")


formula <- as.formula(paste("Surv(time, Overall_Survival_Status) ~", "Group"))

cox_model <- coxph(formula, data = clin.all)

cox_model %>%  gtsummary::tbl_regression(exponentiate = TRUE)




clin.all$am_pathogenicity[clin.all$am_pathogenicity > 50] <- 50
clin.all  %>%
  group_by(Group) %>%  # Optional if you want summary stats or prep before plotting
  ggplot(aes(x = Group, y = am_pathogenicity, fill = Group)) +
  
  geom_violin(trim = FALSE) +
  geom_boxplot(width = 0.1, fill = "white", outlier.shape = NA) +
  labs(
    title = "Violin Plot of Pathscore by Group",
    x = "Group",
    y = "Pathscore"
  ) +
  theme_minimal() +
  theme(legend.position = "none")















# Dynamically apply the factor conversion and releveling to each specified column
clin <- clin %>%
  mutate(across(all_of(columns_to_relevel), 
                ~ factor(.x) ))


# Loop over the columns and relevel
for(i in columns_to_relevel) {
  clin[[i]] <- relevel(clin[[i]], ref = "Top 10% score")
}

clinbk= clin

load("/mmfs1/home/david.adeleke/dissertation/data/processed_data/tcga_us/dataset.RData")

clinbk2=clinical_surv[,c(1,25)]

clin.all = merge(clin , clinbk2, by="icgc_donor_id" )

#saveRDS(clin.all, "/mmfs1/home/david.adeleke/dissertation/data/processed_data/tcga_us/clin.all.rds")
saveRDS(clin.all, "/mmfs1/home/david.adeleke/dissertation/data/processed_data/tcga_us/clin.all.rds")

# Fit survival model
sfit <- survfit(Surv(time, Overall_Survival_Status) ~ class_am_pathogenicity, data = clin.all)

# Plot survival curves
ggsurvplot(
  sfit,
  risk.table = TRUE,
  pval = TRUE,
  conf.int = FALSE,
  #palette = c("#47BC82", "#FAA928", "#B95D4B", "#A18CBA"),
  xlim = c(0, 2500),
  break.time.by = 200,
  #ggtheme = theme_RTCGA(),
  risk.table.y.text.col = TRUE,
  risk.table.y.text = FALSE,
  risk.table.height = 0.4,
  legend.title = "Percentiles")

clin.all$am_pathogenicity[clin.all$am_pathogenicity > 50] <- 50


clin.all  %>%
  group_by(class_am_pathogenicity) %>%  # Optional if you want summary stats or prep before plotting
  ggplot(aes(x = class_am_pathogenicity, y = am_pathogenicity, fill = class_am_pathogenicity)) +
  
  geom_violin(trim = FALSE) +
  geom_boxplot(width = 0.1, fill = "white", outlier.shape = NA) +
  labs(
    title = "Violin Plot of Pathscore by Group",
    x = "Group",
    y = "Pathscore"
  ) +
  theme_minimal() +
  theme(legend.position = "none")




for(i in columns_to_relevel) {
  
  formula <- as.formula(paste("Surv(time, Overall_Survival_Status) ~", i))
  
  cox_model <- coxph(formula, data = clin)
  
  cox_model %>% 
    gtsummary::tbl_regression(exp = TRUE) %>% 
    print()
}




columns_to_relevel_sub = columns_to_relevel[4:6]


for(i in  columns_to_relevel_sub ) {
  
  formula <- as.formula(paste("Surv(OS_MONTHS, status) ~", i))
  
  cox_model <- coxph(formula, data = clin)
  
  t1=cox_model %>% 
    gtsummary::tbl_regression(exponentiate = TRUE) 
  
  assign(i,t1)
}

gtsummary::tbl_merge(
  tbls = list(get( columns_to_relevel_sub [1]),get( columns_to_relevel_sub [2]) ,get( columns_to_relevel_sub [3])),
  tab_spanner = c("**TMB**", "**Alpha**","**AlphaTMB**" ))



# Clinicopathological characteristics of the cohort

clin %>% tbl_summary(include = c(CANCER_TYPE, SEX, status_label, AGE_GROUP))

#saveRDS(clin, file="clin_snv_alpha.RData")

clin=readRDS("clin_snv_alpha.RData")

clin$sTMB_NONSYNONYMOUS=scale(clin$TMB_NONSYNONYMOUS)


clin %>% 
  tbl_summary(
    include = c(CANCER_TYPE, SEX, status_label, AGE_GROUP),
    statistic = list(
      all_continuous() ~ "{mean} ({sd})",
      all_categorical() ~ "{n} / {N} ({p}%)"
    ),
    digits = all_continuous() ~ 2,
    
    missing_text = "(Missing)")



###########Correlation plot
clin %>% 
  ggplot(aes(x = alpham, y = TMB_NONSYNONYMOUS)) +    
  geom_point() +  
  geom_smooth(method="lm", formula=y ~ x) + 
  annotate("text", 
           y = max(clin$TMB_NONSYNONYMOUS) * 0.95,  # Adjust y to avoid overlap
           x = median(clin$alpham), 
           hjust = "inward", vjust = "inward", 
           label = sprintf("spearman rho = %.3f, p-value = %.3f", 
                           cor(clin$alpham, clin$TMB_NONSYNONYMOUS, method="spearman"), 
                           cor.test(clin$alpham, clin$TMB_NONSYNONYMOUS, exact=FALSE, method="spearman")$p.value)) + 
  labs(title = "   ", 
       y = "Tumor Mutation Burden", 
       x = "Alphamissense Score")+
  theme(
    plot.background = element_rect(fill = "white", color = "white"),  # White background
    panel.background = element_rect(fill = "white", color = "white"), # White panel background
    panel.grid.major = element_line(color = "gray", size = 0.1), # Optional: grid lines
    panel.grid.minor = element_blank()  # Optional: remove minor grid lines
  )

###################Pie chart - proportion


# Data
total <- 17539
no_alpha <- 5367
alpha_variants <- 12172

# Create a data frame for ggplot
data <- data.frame(
  category = c("No Alpha", "Alpha Variants"),
  value = c(no_alpha, alpha_variants)
)

# Calculate the percentages
data <- data %>%
  mutate(percentage = value / total * 100,  # Calculate percentage
         label = paste(value, "(", round(percentage, 1), "%)", sep = ""))  # Create label in 'n (percent)' format

ggplot(data, aes(x = "", y = value, fill = category)) +
  geom_bar(stat = "identity", width = 1) +
  coord_polar(theta = "y") +
  labs(title = "Proportion of All TMB variants with Alphamissense scores") +
  theme_void() +  # Remove axis lines and labels
  scale_fill_brewer(palette="Set1") +
  geom_text(aes(label = label), position = position_stack(vjust = 0.5), 
            color = "white", size = 5) +
  theme(legend.position = "top") 

#########kernel density################



library(ggplot2)

p1 =ggplot(data = clin, aes(x = sAlpha, group = status_label, fill = status_label)) +
  geom_density(adjust = 1.5, alpha = 0.4) +
  ggtitle("Density Distribution of Alphamissense Score by Survival Status") +
  xlab("Alphamissense Score") +
  ylab("Density") +
  theme_minimal() +
  theme(
    plot.title = element_text(hjust = 0.5, size = 14, face = "bold"),  # Centered and bold title
    panel.background = element_rect(fill = "white"),
    plot.background = element_rect(fill = "white"),
    legend.title = element_text(face = "bold"),  # Bold legend title for clarity
    legend.position = "top"  # Moves legend to the top
  ) + scale_fill_discrete(name = "Survival Status")

p2 = ggplot(data = clin, aes(x = TMB, group = status_label, fill = status_label)) +
  geom_density(adjust = 1.5, alpha = 0.4) +
  ggtitle("Density Distribution of Alphamissense Score by Survival Status") +
  xlab("Tumor Mutation Burden") +
  ylab("Density") +
  theme_minimal() +
  theme(
    plot.title = element_text(hjust = 0.5, size = 14, face = "bold"),  # Centered and bold title
    panel.background = element_rect(fill = "white"),
    plot.background = element_rect(fill = "white"),
    legend.title = element_text(face = "bold"),  # Bold legend title for clarity
    legend.position = "top"  # Moves legend to the top
  ) + scale_fill_discrete(name = "Survival Status")




























clinical <- clinical %>% filter(project_code %in% proj.scope[i])
snv <- snv %>% filter(project_code %in% proj.scope[i]) %>%
  filter(!(Variant_Classification %in% excl_var))

snv <- as.data.frame(snv)

snv$gr <- paste0(snv$Chromosome, snv$Start_Position, snv$End_Position) 
snv_ <- snv[, c(1, 22, 50)]

misalpha$CHROM_ <- trimws(substring(misalpha$CHROM, 4, 5))
misalpha$gr <- paste0(misalpha$CHROM_, misalpha$POS, misalpha$POS) 
misalpha_ <- misalpha[, c(12, 9)]

joined_df <- merge(snv_, misalpha_, by.x = "gr", by.y = "gr")

snv_alpha <- snv_alpha %>% dplyr::select(c(2, 3, 4)) 
names(snv_alpha) <- c("gene", "id", "pathscore")
snv_alpha=snv_alpha[!is.na(snv_alpha$pathscore),]
snv_alpha_ <- snv_alpha %>% select(id, pathscore)  %>%
  group_by(id) %>%
  mutate(score = sum(pathscore, na.rm = TRUE))











ref <- "KnTn"
clin$cds_mutation=as.character(clin$cds_mutation)
clin$cds_mutation<- as.factor(clin$cds_mutation)


clin$grp=as.character(clin$grp)
clin$grp <- as.factor(clin$grp)
clin <- within(clin, grp<- relevel(grp, ref = ref))


clin <- within(clin, cds_mutation<- relevel(cds_mutation, ref = ref))


"/mmfs1/home/david.adeleke/Thesis/dataset/clin_worked.RData"
save(clin, file="clin_worked.RData")


coxph(Surv(time, status) ~ cds_mutation, data = clin[clin$time >=1 & !(clin$cds_mutation %in% c("131A>G","34G>A", "35G>C","37G>T" , "34G>T")),]) %>% 
  gtsummary::tbl_regression(exp = TRUE)


coxph(Surv(time, status) ~ grp, data = clin[clin$time >=1 & !(clin$cds_mutation %in% c("131A>G","34G>A", "35G>C","37G>T" , "34G>T")),]) %>% 
  gtsummary::tbl_regression(exp = TRUE)




sfit <- survfit(Surv(time, status)~ cds_mutation, data = clin[clin$time >=1 & !(clin$cds_mutation %in% c("131A>G","34G>A", "35G>C","37G>T" , "34G>T")),])


ggsurvplot(
  sfit,                     # survfit object with calculated statistics.
  pval = TRUE,             # show p-value of log-rank test.
  conf.int = F,         # show confidence intervals for 
  # point estimaes of survival curves.
  #conf.int.style = "step",  # customize style of confidence intervals
  xlab = "Time in days",   # customize X axis label.
  break.time.by = 200,     # break X axis in time intervals by 200.
  ggtheme = theme_light(), # customize plot and risk table with a theme.
  risk.table = "abs_pct",  # absolute number and percentage at risk.
  risk.table.y.text.col = T,# colour risk table text annotations.
  risk.table.y.text = FALSE,# show bars instead of names in text annotations
  # in legend of risk table.
  ncensor.plot = F,      # plot the number of censored subjects at time t
  surv.median.line = "hv",  # add the median survival pointer.
  #legend.labs = 
  #c("Low Risk", "High Risk"),    # change legend labels.
  palette =  c("blue","red",  "black", "orange", "green" ) )
