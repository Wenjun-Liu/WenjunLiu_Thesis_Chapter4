tumorSample <- read.delim("/Users/wenjunliu/20131906_HickeyT_JC_NormalBreast/config/tumor_sample.csv", header = TRUE)
normalSample <- read.delim("/Users/wenjunliu/20131906_HickeyT_JC_NormalBreast/config/normal_samples.tsv", header = TRUE)
combineSample <- tumorSample %>%
    mutate(patient = vapply(.$SampleName, function(x){
        str_split(x, "-")[[1]][1]
        }, character(1)),
    patient = str_replace(patient, "TH", "TH-"), 
    sample = str_remove(Filename, ".r_1.fq.gz"), 
    sample = str_remove(sample, ".tophat.homo_sapiens.bam"), 
    flowcell_id = str_sub(SourceID, 20,28), 
    treat = vapply(.$SampleName, function(x){
        str_split(x, "-")[[1]][2]
    }, character(1)), 
    treat = case_when(treat == "D" ~ "DHT", 
                      treat == "E" ~ "E2", 
                      treat == "ED" ~ "E2+DHT", 
                      treat == "V" ~ "Veh"), 
    Tumor = TRUE, 
    desc = paste(patient, treat,sep = " ")) %>%
    dplyr::rename(name = SourceID) %>%
    dplyr::select(-c("SampleName", "Filename")) %>%
    rbind(normalSample %>%
              mutate(Tumor = FALSE))
write_tsv(combineSample, file = "/Users/wenjunliu/20131906_HickeyT_JC_NormalBreast/config/samples.tsv")    
