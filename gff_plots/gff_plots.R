# Plots sequence coordinates from GFF format files in multiple tracks 

# Load libraries
library(ggplot2)
library(gridExtra)
library(yaml)

# Read configuration file
config <- yaml.load_file("config.yml")

# Load GFF files
list_gff <- lapply(config$gff, function(x) {as.data.frame(read.csv(x, sep = "\t", header = 0))})

# Config scientific notation
options(scipen = 10000)

# Set up color scheme and labels
sw_labels <- config$labels
sw_colors <- config$colors 
plot_colors <- setNames(sw_colors, sw_labels)

# Function to draw plots in multiple tracks
create_plot <- function(xmin, xmax, list_gff) {
    daf <- data.frame()
    #Base plot to set coordinates
    base_plot <- ggplot(daf) + geom_point() + 
                 coord_cartesian(ylim=c(0.5, 7.5), xlim=c(xmin, xmax), expand=FALSE) +
                 scale_y_continuous(limits=c(1, 8), breaks=seq(1, 7), labels=sw_labels) 
    #Iterate over GFFs to draw each one in a track
    y_pos <- 0
    gplot <- lapply(seq_along(list_gff), function(i){
                    x <- list_gff[[i]]
                    geom_segment(aes(x=V4, xend=V5, y=i, yend=i, colour=V2),
                    data=subset(x,V4>xmin & V5<xmax), size=3, alpha=0.7)
	      })
    main_plot <- base_plot + gplot
    theme(axis.ticks=element_blank(), axis.title.y=element_blank(),
          axis.title.x=element_blank(), legend.position="none") +
       theme(plot.margin=unit(c(0.3, 1, 0.3, 1), "cm")) + scale_colour_manual(values=plot_colors)
    return(main_plot)
}

# Set coordinates for starting point, step and number of plots
start_coor <- config$coor[['start']]
step_size <- config$coor[['step']]
repeats <- config$coor[['repeat']]

# Iterate over function to create multiple plots
glist <- lapply(seq(1:repeats), function(x){create_plot(start_coor, start_coor + step_size*x, list_gff)})

# Save plot to file and show it
png("plot.png")
grid.arrange(ncol=1, grobs=glist)
dev.off()
browseURL("plot.png")
