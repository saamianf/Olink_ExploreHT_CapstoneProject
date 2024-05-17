library(readxl)
library(ggVennDiagram)
library(ggplot2)
library(gplots)
library(plyr)
library(data.table)

data <- read_excel('Olink_panels_comparisons.xlsx')
ExploreHT <- na.omit(data[,1])
OncologyIII_T96 <- na.omit(data[,2])
ImmuneResponse_T96 <- na.omit(data[,3])
ImmunoOncology_T96 <- na.omit(data[,4])
Flex48_A <- na.omit(data[,5])
Flex48_B <- na.omit(data[,6])

venn_fn <- function(x, col1, col2, atitle){
  venn_dia <- ggVennDiagram(x, label_alpha=0) + scale_fill_gradient(low=col1, high=col2)  + scale_x_continuous(expand = expansion(mult = .2)) + labs(title=atitle, caption=Sys.Date())
  dev.new()
  print(venn_dia)
  venn_data <- venn(x, show.plot=FALSE)
  venn_int <- attr(venn_data, 'intersection')
  venn_table <- transpose(ldply(venn_int, rbind))
  full_atitle <- paste(atitle, '_Proteins.csv', sep='')
  write.csv(venn_table, full_atitle, na='')
}

all_venn <- venn_fn(c(ExploreHT, OncologyIII_T96, ImmuneResponse_T96, ImmunoOncology_T96, Flex48_A, Flex48_B), 'white', 'green', 'Olink_Panels')

library(eulerr)

fit <- euler(c('ExploreHT&IR' = 74, 'ExploreHT&OncIII' = 80, 'ExploreHT&Flex48A' = 6, 'ExploreHT&Flex48B' = 3,
               'ExploreHT&ImmunoOnc' = 41, 'IR' = 5, 'ExploreHT&IR&ImmunoOnc' = 1, 'ExploreHT&OncIII&Flex48B' = 1,
               'ExploreHT&Flex48A&Flex48B' = 1, 
               'OncIII' = 11, 'ImmunoOnc' = 26, 'ExploreHT' = 5112,  
               'OncIII&ImmunoOnc' = 0, 'IR&ImmunoOnc' = 0, 
               'IR&OncIII' = 0, 
               'ExploreHT&IR&Flex48A' = 1, 
               'ExploreHT&ImmunoOnc&Flex48A' = 5, 'ExploreHT&ImmunoOnc&Flex48B' = 8, 'ExploreHT&IR&ImmunoOnc&Flex48A' = 4,
               'IR&ImmunoOnc&Flex48A&Flex48B' = 1, 'ExploreHT&IR&ImmunoOnc&Flex48B' = 3, 
               'ExploreHT&IR&ImmunoOnc&Flex48A&Flex48B' = 3, 'Flex48B' = 1))
# values in above are from the saved csv
pl <- plot(fit, quantities=list(cex=2), legend = TRUE, 
           fills = c('lightgrey', 'lightgreen', 'lightblue', 'purple', 'red', 'orange'),
           adjust_labels = TRUE)
dev.new()
print(pl)