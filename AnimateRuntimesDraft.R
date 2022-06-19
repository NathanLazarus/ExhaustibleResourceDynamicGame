library(ggplot2)
library(gganimate)
library(gifski)
library(cowplot)
library(RColorBrewer)
library(stringr)
library(av)
library(magick)
library(foreach)

asdf = data.table(raw = readLines("Wavefront.log"))
hm = asdf[grep("Job submitted from host", raw), .(submission = raw)]
hm[, DAG_node := asdf[grep("DAG Node", raw), raw]]
hm[, stage := tstrsplit(DAG_node, "_")[2]]
hm[, job := tstrsplit(DAG_node, "_")[4]]
hm[, `:=`(stage = as.integer(stage), job = as.integer(job))]
hm[, submitted := str_extract(submission, "2021-[0-9]{2}-[0-9]{2} [0-9]{2}:[0-9]{2}:[0-9]{2}")]
hm[, jobID := str_extract(submission, "(?<=000 \\()[0-9]{8}\\.[0-9]{3}\\.[0-9]{3}(?=\\))")]

hm2 = asdf[grep("Job terminated\\.", raw), .(termination = raw)]
hm2[, jobID := str_extract(termination, "(?<=005 \\()[0-9]{8}\\.[0-9]{3}\\.[0-9]{3}(?=\\))")]
hm2[, finished := str_extract(termination, "2021-[0-9]{2}-[0-9]{2} [0-9]{2}:[0-9]{2}:[0-9]{2}")]

hm[hm2, on = .(jobID), finished := i.finished]
data = hm


# data = data.table(
#   stage = rep(1:5, times = 1:5),
#   submitted = c(
#     "2021-08-08 04:05:43",
#     "2021-08-08 04:05:44",
#     "2021-08-08 04:05:49",
#     "2021-08-08 04:05:57",
#     "2021-08-08 04:06:24",
#     "2021-08-08 04:05:40",
#     "2021-08-08 04:07:24",
#     "2021-08-08 04:07:34",
#     "2021-08-08 04:07:44",
#     "2021-08-08 04:08:24",
#     "2021-08-08 04:08:34",
#     "2021-08-08 04:08:44",
#     "2021-08-08 04:09:24",
#     "2021-08-08 04:09:34",
#     "2021-08-08 04:09:44"
#   ),
#   finished = c(
#     "2021-08-08 04:07:43",
#     "2021-08-08 04:07:44",
#     "2021-08-08 04:14:49",
#     "2021-08-08 04:07:57",
#     "2021-08-08 04:06:54",
#     "2021-08-08 04:06:40",
#     "2021-08-08 04:08:24",
#     "2021-08-08 04:08:34",
#     "2021-08-08 04:12:44",
#     "2021-08-08 04:08:24",
#     "2021-08-08 04:09:34",
#     "2021-08-08 04:12:44",
#     "2021-08-08 04:11:24",
#     "2021-08-08 04:11:34",
#     "2021-08-08 04:11:44"
#   )
# )

data[, size := 1
   ][, id := .I
   ][, submitted := as.POSIXct(submitted)
   ][, finished := as.POSIXct(finished)]

# all_times = setkey(data.table(time_considered = unique(c(data$submitted, data$finished))), time_considered)
# 
# all_times[, include := (.I %% 50 == 0 | .I == 1)]

dur = 20
fps = 10
all_times = data.table(time_considered = seq(min(data$submitted) - 1, max(data$finished), length.out = dur * fps))

CJ.dt = function(X,Y) {
  stopifnot(is.data.table(X),is.data.table(Y))
  k = NULL
  X = X[, c(k = 1, .SD)]
  setkey(X, k)
  Y = Y[, c(k = 1, .SD)]
  setkey(Y, NULL)
  X[Y, allow.cartesian = TRUE][, k := NULL][]
}

cartesian_prod = CJ.dt(all_times, data)
setorder(cartesian_prod, time_considered, stage, -submitted, -finished)

cartesian_prod[, status_at_time := as.character((submitted <= time_considered) + (finished <= time_considered))]

myColors = c("grey", "cornflowerblue", "forestgreen")
names(myColors) = c("0", "1", "2")
fillScale = scale_fill_manual(name = "status_at_time", values = myColors)

semi_transparent_borders = paste0("#000000", format(as.hexmode(round(seq(255, 40, length.out = 11))), width = 2, upper.case = TRUE))
myColors = c(semi_transparent_borders, rep("#00000000", max(data$stage) - 2 * length(semi_transparent_borders)), rev(semi_transparent_borders))
names(myColors) = as.character(1:max(data$stage))
colScale = scale_color_manual(name = "stage", values = myColors)


cartesian_prod[, help := .GRP, time_considered]
cartesian_prod[, id := .I]
cartesian_prod[, seconds_elapsed := as.integer(time_considered - min(time_considered))]
cartesian_prod[, hours_minutes_elapsed := as.character(round_hms(as_hms(time_considered - min(time_considered)), 60))]
cartesian_prod[, hm_elapsed := substr(hours_minutes_elapsed, 1, nchar(hours_minutes_elapsed) - 3)]

timediffs = as.integer(all_times[, time_considered - shift(time_considered)][-1])

plotlist = foreach(this_help = unique(cartesian_prod$help), .combine = c) %do% {
ggplot(cartesian_prod[help == this_help]) +
  geom_col(aes(x = stage, y = size, fill = status_at_time, group = id, color = as.character(stage)), position = "fill", width = 0.79) + # , color = "black"
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
  ggtitle(paste0("Time Elapsed   ", cartesian_prod[help == this_help, first(hm_elapsed)]))
  filename = paste0("Images/WavefrontPlot-", this_help - 1, ".png")
  ggsave(filename)
  filename
}

img_list = lapply(plotlist, image_read)
img_animated = image_animate(image_join(img_list), fps = fps, loop = 1)
# img_animated
image_write(image = img_animated, path = "testmagick.gif")








test = ggplot(cartesian_prod) +
  geom_col(aes(x = stage, y = size, fill = status_at_time, group = id), position = "fill") + # , color = "black"
  colScale +
  theme_cowplot() +
  theme(axis.line = element_blank(),
        axis.text.x = element_blank(),
        axis.text.y = element_blank(),
        axis.ticks = element_blank(),
        axis.title.x = element_blank(),
        axis.title.y = element_blank(),
        legend.position = "none",
        plot.title = element_text(hjust = 0)) +
  transition_states(time_considered, transition_length = 0,
                    state_length = timediffs
                    ) +
  labs(title = "Time Elapsed: {closest_state}") # +
  # transition_events(time_considered) +
  # transition_events(start = as.integer(submitted), enter_length = 0L, exit_length = 0L) # +
  # ease_aes("linear")

b = animate(test, duration = dur, fps = fps, renderer = av_renderer())
anim_save("testing.mp4", b)

# anim_save("test2.gif", b, renderer = gifski_renderer(loop = FALSE))

# dat <- data.table(
#   x = 1:3,
#   y = runif(3),
#   begin = runif(3, 1, 100),
#   length = runif(3, 5, 20),
#   enter = runif(3, 5, 10),
#   exit = runif(3, 5, 10)
# )
# 
# anim <- ggplot(dat, aes(x, y)) +
#   geom_col() +
#   transition_events(start = begin,
#                     end = begin + length,
#                     enter_length = enter,
#                     exit_length = exit)

# library(gapminder)
# 
# ggplot(gapminder, aes(gdpPercap, lifeExp, size = pop, colour = country)) +
#   geom_point(alpha = 0.7, show.legend = FALSE) +
#   scale_colour_manual(values = country_colors) +
#   scale_size(range = c(2, 12)) +
#   scale_x_log10() +
#   facet_wrap(~continent) +
#   # Here comes the gganimate specific bits
#   labs(title = "Year: {frame_time}", x = "GDP per capita", y = "life expectancy") +
#   transition_time(year) +
#   ease_aes("linear")
# 
# anim_save("test.gif")
# 
# sam<-data.frame(population=c(rep("PRO",8),rep("SOM",4)),
#                 allele=c("alele1","alele2","alele3","alele4",rep("alele5",2),
#                             rep("alele3",2),"alele2","alele3","alele3","alele2"), 
#                 frequency=rep(c(10,5,4,6,7,16),2) #,rep(1,6)))
# )
# 
# sam <- setDT(sam)[, .(frequencySum=sum(frequency)), by=.(population,allele)]
# 
# sam <- sam[order(sam$frequency, decreasing = TRUE),] # (1)
# 
# # (2)
# sam$frequency<-factor(sam$frequency, levels = unique(sam$frequency) )
# 
# library(ggplot2)
# ggplot(sam)+
#   geom_bar(aes(x=population, y=frequencySum, group=frequency, fill=allele), # (3)
#            stat="identity", color="white") 
# 
# 
# 
# 
# 
# d <- read.table(text="Day Length Amount
#                 1 3 1
#                 1 4 2
#                 3 3 2
#                 3 5 1",header=T)
# 
# 
# d$Amount<-as.factor(d$Amount) 
# 
# d <- d[order(d$Length, decreasing = TRUE),] # (1)
# 
# d$LengthFactor<-factor(d$Length, levels= unique(d$Length) ) # (2)
# 
# ggplot(d)+
#   geom_bar(aes(x=Day, y=Length, group=LengthFactor, fill=Amount), # (3)
#            stat="identity", color="white") 