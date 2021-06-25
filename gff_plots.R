# Plots sequence coordinates in GFF format in multiple tracks Each track
# corresponds to a GFF file

#Load libraries
library(ggplot2)
library(gridExtra)

# Set GFF files as variables
gff_ref_file = "chr21.all.fil.gff"
gff_cmp_file_3 = "chr21_RepeatScout_RMasker.gff"
gff_cmp_file_2 = "chr21_RepeatModeler_RMasker.gff"
gff_cmp_file_1 = "chr21_REPET_RMasker.gff"
gff_cmp_file_4 = "chr21.phRaider.gff"
gff_cmp_file_5 = "chr21.Pclouds.gff"
gff_cmp_file_6 = "chr21.Red.gff"

# Load GFF files into multiple data frames
gff_ref = as.data.frame(read.csv(gff_ref_file, sep = "\t", header = 0))
gff_cmp_1 = as.data.frame(read.csv(gff_cmp_file_1, sep = "\t", header = 0))
gff_cmp_2 = as.data.frame(read.csv(gff_cmp_file_2, sep = "\t", header = 0))
gff_cmp_3 = as.data.frame(read.csv(gff_cmp_file_3, sep = "\t", header = 0))
gff_cmp_4 = as.data.frame(read.csv(gff_cmp_file_4, sep = "\t", header = 0))
gff_cmp_5 = as.data.frame(read.csv(gff_cmp_file_5, sep = "\t", header = 0))
gff_cmp_6 = as.data.frame(read.csv(gff_cmp_file_6, sep = "\t", header = 0))

# Config scientific notation
options(scipen = 10000)

# Set up color scheme
green = "#4DAF4A"
blue = "#377EB8"
purple = "#984EA3"
red = "#E41A1C"
yellow = "#F79616"
cyan = "#08B2C4"
brown = "#C73302"
sw_labels = c("reference", "REPET", "RModeler", "RScout", "phRaider", "Pclouds", "Red")
sw_colors = c("REPET" = blue, "RepeatModeler" = purple, "RepeatScout" = red, "RepeatMasker" = green, 
	      "Red" = brown, "phRaider" = yellow, "Pclouds" = cyan)

#Function to draw plots in multiple tracks
create_plot <- function(xmin, xmax) {
    df <- data.frame()
    base_plot <- ggplot(df) + geom_point() + coord_cartesian(ylim=c(0.5,7.5), xlim=c(xmin, xmax), expand = FALSE)  +
        scale_y_continuous(limits=c(1,8), breaks=seq(1,7),labels=sw_labels) 
    main_plot <- base_plot + 
        geom_segment(aes(x=xmin, xend=xmax,y=1,yend=1),data=df, size=0.7, alpha=0.5) + 
        geom_segment(aes(x=V4,xend=V5,y=1,yend=1,colour=V2),data=subset(gff_ref,V4>xmin & V5<xmax),size=3) + 
        geom_segment(aes(x=V4,xend=V5,y=2,yend=2,colour=V2),data=subset(gff_cmp_1, V4>xmin & V5<xmax),size=3, alpha=0.7) +
        geom_segment(aes(x=V4,xend=V5,y=3,yend=3,colour=V2),data=subset(gff_cmp_2, V4>xmin & V5<xmax),size=3, alpha=0.7) +
        geom_segment(aes(x=V4,xend=V5,y=4,yend=4,colour=V2),data=subset(gff_cmp_3, V4>xmin & V5<xmax),size=3, alpha=0.7) +
        geom_segment(aes(x=V4,xend=V5,y=5,yend=5,colour=V2),data=subset(gff_cmp_4, V4>xmin & V5<xmax),size=3, alpha=0.7) +
        geom_segment(aes(x=V4,xend=V5,y=6,yend=6,colour=V2),data=subset(gff_cmp_5, V4>xmin & V5<xmax),size=3, alpha=0.7) +
        geom_segment(aes(x=V4,xend=V5,y=7,yend=7,colour=V2),data=subset(gff_cmp_6, V4>xmin & V5<xmax),size=3, alpha=0.7) +
       theme(
          axis.ticks = element_blank(),axis.title.y = element_blank(),
          axis.title.x = element_blank(), legend.position = "none") +
       theme(plot.margin = unit(c(0.3, 1, 0.3, 1), "cm")) + scale_colour_manual(values = sw_colors)
  return(main_plot)
}

#Set starting point and step of coordinates to plot
start <- 9000000
step <- 50000

# Create plots
p1 <- create_plot(start, start + step)
p2 <- create_plot(start + step, start + step * 2)
p3 <- create_plot(start + step * 2, start + step * 3)

#Display plots in grid
grid.arrange(ncol = 1, p1, p2, p3)

