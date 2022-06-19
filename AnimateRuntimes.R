# AnimateRuntimes.R
# Creates an animated GIF of the progress of Condor jobs

# Parameters ------------------------------------------------------------------------

duration_of_anim = 30 # seconds
fps = 10

# Parse the log file to get jobs' submit time, finishing time, and batch ID ---------
raw_log = data.table(raw = readLines("Wavefront.log"))

data = raw_log[grep("Job submitted from host", raw), .(submission = raw)]
data[, DAG_node := raw_log[grep("DAG Node", raw), raw]]
data[, stage := tstrsplit(DAG_node, "_")[2]]
data[, job := tstrsplit(DAG_node, "_")[4]]
data[, `:=`(stage = as.integer(stage), job = as.integer(job))]
data[, submitted := str_extract(submission, "2021-[0-9]{2}-[0-9]{2} [0-9]{2}:[0-9]{2}:[0-9]{2}")]
data[, jobID := str_extract(submission, "(?<=000 \\()[0-9]{8}\\.[0-9]{3}\\.[0-9]{3}(?=\\))")]

job_terminations = raw_log[grep("Job terminated\\.", raw), .(termination = raw)]
job_terminations[, jobID := str_extract(termination, "(?<=005 \\()[0-9]{8}\\.[0-9]{3}\\.[0-9]{3}(?=\\))")]
job_terminations[, finished := str_extract(termination, "2021-[0-9]{2}-[0-9]{2} [0-9]{2}:[0-9]{2}:[0-9]{2}")]

data[job_terminations, on = .(jobID), finished := i.finished]

# Format things so they'll work in ggplot --------------------------------------------

data[, size := 1
   ][, id := .I
   ][, submitted := as.POSIXct(submitted)
   ][, finished := as.POSIXct(finished)]

# store every grid point that will be required for the animation and create a long
# data table that includes 
time_grid = data.table(time_considered = seq(min(data$submitted) - 1, max(data$finished),
                                             length.out = duration_of_anim * fps))

cartesian_prod = CJ.dt(time_grid, data)
setorder(cartesian_prod, time_considered, stage, -submitted, -finished)

cartesian_prod[, status_at_time := as.character((submitted <= time_considered) + (finished <= time_considered))]

cartesian_prod[, time_grid_ind := .GRP, time_considered
             ][, id := .I
             ][, seconds_elapsed := as.integer(time_considered - min(time_considered))
             ][, hours_minutes_elapsed := as.character(round_hms(as_hms(time_considered - min(time_considered)), 60))
             ][, hm_elapsed := substr(hours_minutes_elapsed, 1, nchar(hours_minutes_elapsed) - 3)]

# Create fill colors that correspond to status ---------------------------------------

myColors = c("grey", "cornflowerblue", "#107010FF")
names(myColors) = c("0", "1", "2")
fillScale = scale_fill_manual(name = "status_at_time", values = myColors)

# Create black borders that fade out as the job density increases --------------------

semi_transparent_borders = paste0("#000000", format(as.hexmode(round(seq(255, 40, length.out = 11))),
                                                    width = 2, upper.case = TRUE))
myColors = c(semi_transparent_borders,
             rep("#00000000", max(data$stage) - 2 * length(semi_transparent_borders)), #fully transparent in the middle
             rev(semi_transparent_borders))
names(myColors) = as.character(1:max(data$stage))
colScale = scale_color_manual(name = "stage", values = myColors)


# Plot every frame -------------------------------------------------------------------

clusters = makeCluster(7)
registerDoSNOW(clusters)

plotlist = foreach(this_time_grid_ind = unique(cartesian_prod$time_grid_ind), .combine = c) %dopar% {
  data_here = cartesian_prod[time_grid_ind == this_time_grid_ind]
  ggplot(data_here) +
    geom_col(aes(x = stage, y = size, fill = status_at_time, group = id, color = as.character(stage)),
             position = "fill", width = 0.79) +
    colScale +
    fillScale +
    theme_cowplot() +
    theme(axis.line = element_blank(),
          axis.text.x = element_blank(),
          axis.text.y = element_blank(),
          axis.ticks = element_blank(),
          axis.title.x = element_blank(),
          axis.title.y = element_blank(),
          # plot.background = element_rect(fill = "white"),
          legend.position = "none",
          plot.title = element_text(hjust = 0)) +
    coord_cartesian(expand = FALSE) +
    ggtitle(paste0("Time Elapsed   ", data_here[, first(hm_elapsed)]))
    filename = paste0("Images/WavefrontPlot-", this_time_grid_ind - 1, ".png")
    ggsave(filename)
    filename
}

stopCluster(clusters)


# Knit the frames together into an animation using magick -----------------------------

all_images = lapply(plotlist, image_read)
animated = image_animate(image_join(all_images), fps = fps, loop = 1)
image_write(image = animated, path = "AnimateRuntimes.gif")
