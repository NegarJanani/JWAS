library(ggplot2)
install_github("juliedwhite/miamiplot", build_vignettes = TRUE)
library(miamiplot)
vignette("miamiplot")

# Miamai Plot for Airway AA

MIMAAA<- read.csv("~/Desktop/UQ/AA/MIMAAA.csv")

MIMAAA$X <- NULL

my_upper_colors <- RColorBrewer::brewer.pal(4, "Paired")[1:2]
my_lower_colors <- RColorBrewer::brewer.pal(4, "Paired")[3:4]


ggmiami(data = MIMAAA , split_by = "study", split_at = "DCI", 
        upper_ylab = "DCI Results", lower_ylab = "GWAS Results", chr_colors = NULL,
        upper_chr_colors = my_upper_colors, lower_chr_colors = my_lower_colors,
        genome_line_color = "red", suggestive_line_color = "#A9A9A9",  genome_line = 5e-8)



# Miamai Plot for Emphysema AA


MIMEAA<- read.csv("~/Desktop/UQ/AA/MIMEAA.csv")

MIMEAA$X <- NULL

my_upper_colors <- RColorBrewer::brewer.pal(4, "Paired")[1:2]
my_lower_colors <- RColorBrewer::brewer.pal(4, "Paired")[3:4]


ggmiami(data = MIMEAA , split_by = "study", split_at = "DCI", 
        upper_ylab = "DCI Results", lower_ylab = "GWAS Results", chr_colors = NULL,
        upper_chr_colors = my_upper_colors, lower_chr_colors = my_lower_colors,
        genome_line_color = "red", suggestive_line_color = "#A9A9A9", genome_line = 5e-8)



# Miamai Plot for Airway NHW

MIMANHW<- read.csv("~/Desktop/UQ/NHW/MIMANHW.csv")

MIMANHW$X <- NULL

my_upper_colors <- RColorBrewer::brewer.pal(4, "Paired")[1:2]
my_lower_colors <- RColorBrewer::brewer.pal(4, "Paired")[3:4]


ggmiami(data = MIMANHW , split_by = "study", split_at = "DCI", 
        upper_ylab = "DCI Results", lower_ylab = "GWAS Results", chr_colors = NULL,
        upper_chr_colors = my_upper_colors, lower_chr_colors = my_lower_colors,genome_line_color = "red", 
        hits_label = c("rs2941504") ,
        hits_label_col = "rsid",suggestive_line_color = "#A9A9A9", genome_line = 5e-8)


# Miamai Plot for Emphysema NHW

MIMENHW<- read.csv("~/Desktop/UQ/NHW/MIMENHW.csv")

MIMENHW$X <- NULL

my_upper_colors <- RColorBrewer::brewer.pal(4, "Paired")[1:2]
my_lower_colors <- RColorBrewer::brewer.pal(4, "Paired")[3:4]


ggmiami(data = MIMENHW , split_by = "study", split_at = "DCI", 
        upper_ylab = "DCI Results", lower_ylab = "GWAS Results", chr_colors = NULL,
        upper_chr_colors = my_upper_colors, lower_chr_colors = my_lower_colors,
        genome_line_color = "red", suggestive_line_color = "#A9A9A9", genome_line = 5e-8)

