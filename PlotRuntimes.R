library(foreach)
library(ggplot2)

setwd('C:/Users/Nathan/Downloads/Exhaustible Resource Dynamic Game/FromCHTC')

output_files = list.files(pattern = 'iter.*')

iter_times = foreach(file = output_files, .combine = rbind)%do%{
  fread(file)
}
iter_times = rbind(iter_times, data.table('2021-09-21 19:17:32.5193', 4))

setnames(iter_times, c('time', 'iter'))
setorder(iter_times, -iter)
iter_times[, posixTime := as.POSIXct(time)]
iter_times[, running_time := as.numeric(shift(posixTime, n = -1) - posixTime, units = 'mins')]
Njobs = c(1, 5, 15, 35, 70, 126, 210, 330, 491, 695, 941, 1225, 1540, 1876, 2220, 2556, 2871, 3155, 3401, 3605, 3766, 3886, 3970, 4026, 4061, 4081, 4091, 4095, 4096) -
  c(0, 1, 5, 15, 35, 70, 126, 210, 330, 491, 695, 941, 1225, 1540, 1876, 2220, 2556, 2871, 3155, 3401, 3605, 3766, 3886, 3970, 4026, 4061, 4081, 4091, 4095)
iter_times[, Njobs := Njobs][, job_num := .I]
iter_times[17, running_time := running_time - 38]
iter_times[29, running_time := 5.227]

ggplot(data = iter_times) +
  geom_bar(aes(x = job_num, y = Njobs), stat = 'identity') +
  geom_line(aes(x = job_num, y = running_time / 0.033), color = 'red', size = 2) + 
  scale_y_continuous(
    "Number of Jobs", 
    sec.axis = sec_axis(~ . * 0.033, name = "Running Time (minutes)")
  )+
  theme(axis.title.x=element_blank(),
        axis.text.x=element_blank(),
        axis.ticks.x=element_blank())

ggsave('DAGmanQueuingActualRun.png')
