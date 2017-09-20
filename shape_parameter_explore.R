library(readr)
library(stringr)
library(dplyr)
library(ggplot2)
library(reshape2)

shape_parameter <- read_csv("C:/Users/srwisots/Google Drive/Muse/Mitochondria/Mito_Proj/shape_parameter.csv")

shape_parameter <- mutate(shape_parameter,
                 Order = toupper(str_extract_all(shape_parameter$FILE, "\\w+(?=-)", simplify = T)[,1]),
                 gene = toupper(str_extract_all(shape_parameter$FILE, "\\w+(?=-)", simplify = T)[,2]),
                 Model = str_extract_all(shape_parameter$FILE, "\\w+(?=\\.)", simplify = T)[,4],
                 cat = str_extract_all(shape_parameter$FILE, "\\d+x\\d+", simplify = T)[,1]
)


shape_parameter$Syn.Shape <- as.numeric(str_extract_all(shape_parameter$Syn.Shape, "[-+]?\\d*\\.\\d+|[-+]?\\d+", simplify = T)
                                        [,1])
shape_parameter$Nonsyn.Shape <- as.numeric(str_extract_all(shape_parameter$Nonsyn.Shape, "[-+]?\\d*\\.\\d+|[-+]?\\d+", simplify = T)
                                           [,1])

shape_parameter <- rename(shape_parameter,File=FILE)

shape_parameter$Order[which(shape_parameter$Order == "CARNIVORES")] <-"CARNIVORA"
shape_parameter$Order[which(shape_parameter$Order == "GASTEROSTEIFORMES")] <-"GASTEROSTEALES"
shape_parameter$Order[which(shape_parameter$Order == "CHIMAERIFORMS")]<-"CHIMAERIFORMES"


#Stats and Data doesn't have the 3x3 4x4 orders

Stats_and_Data_56_orders_mtDNA <- read_csv("C:/Users/srwisots/Google Drive/Muse/Mitochondria/Mito_Proj/Stats_and_Data_56_orders_mtDNA.csv", 
                                           col_types = cols(X1 = col_skip(), avg.seq.length = col_double()))
Stats_and_Data_56_orders_mtDNA$cat <- str_replace_all(Stats_and_Data_56_orders_mtDNA$cat, " ","")

same_rate_cats_all <- read_csv("C:/Users/srwisots/Google Drive/Muse/Mitochondria/Mito_Proj/same_rate_cats_all.csv", 
                               col_types = cols(X1 = col_skip()))

same_rate_cats_all$Order <- toupper(same_rate_cats_all$Order)
same_rate_cats_all$gene <- toupper(same_rate_cats_all$gene)


same_rate_cats_all$Order[which(same_rate_cats_all$Order == "CARNIVORES")] <-"CARNIVORA"
same_rate_cats_all$Order[which(same_rate_cats_all$Order == "GASTEROSTEIFORMES")] <-"GASTEROSTEALES"
same_rate_cats_all$Order[which(same_rate_cats_all$Order == "CHIMAERIFORMS")]<-"CHIMAERIFORMES"

same_min <- same_rate_cats_all %>% filter(Model == "Dual Variable Rates") %>% select(Order,gene,cat,Synonymous.CV,NS.CV)
Stats_and_Data_min <-Stats_and_Data_56_orders_mtDNA %>% filter(Model == "Dual Variable Rates") %>% select(Order,gene,cat,Synonymous.CV,NS.CV)

x <- bind_rows(same_min, Stats_and_Data_min)

x <- full_join(x, shape_parameter, by = c("Order","gene","cat"))

x <- mutate(x, syn.rate.cats = as.numeric(str_extract_all(x$cat, "\\d+(?=x)", simplify = T)[,1]),
            nonsyn.rate.cats = as.numeric(str_extract_all(x$cat, "(?<=x)\\d+", simplify = T)[,1]))

x %>% filter((cat %in% c("3x3","4x4","5x5","7x7","10x10"))) %>% ggplot(aes(x = Syn.Shape,y = Nonsyn.Shape)) + 
  geom_point()  + facet_wrap(~cat)

x %>% filter((cat %in% c("3x3","4x4","5x5","7x7","10x10"))) %>% ggplot(aes(x = Syn.Shape,y = Synonymous.CV)) + 
  geom_point()  + facet_wrap(~cat)

x %>% filter((cat %in% c("3x3","4x4","5x5","7x7","10x10"))) %>% ggplot(aes(x = Nonsyn.Shape,y = NS.CV)) + 
  geom_point()  + facet_wrap(~cat)

x %>% ggplot(aes(x =syn.rate.cats, y = Syn.Shape)) + geom_boxplot(aes(group = syn.rate.cats))

x %>% ggplot(aes(x =nonsyn.rate.cats, y = Nonsyn.Shape)) + geom_boxplot(aes(group = nonsyn.rate.cats))

x %>% ggplot(aes(x =syn.rate.cats, y = Synonymous.CV)) + geom_boxplot(aes(group = syn.rate.cats))

x %>% ggplot(aes(x =nonsyn.rate.cats, y = NS.CV)) + geom_boxplot(aes(group = nonsyn.rate.cats))

x <- mutate(x, calculated.syn.cv = sqrt(Syn.Shape)/Syn.Shape,
                calculated.nonsyn.cv = sqrt(Nonsyn.Shape)/Nonsyn.Shape)

x %>% ggplot(aes(x=calculated.syn.cv,y=Synonymous.CV))+geom_point()+geom_abline()
x %>% ggplot(aes(x=calculated.nonsyn.cv,y=NS.CV))+geom_point()+geom_abline()+facet_wrap(~nonsyn.rate.cats)
