Njobs = c(1, 5, 15, 35, 70, 126, 210, 330, 491, 695, 941, 1225, 1540, 1876, 2220, 2556, 2871, 3155, 3401, 3605, 3766, 3886, 3970, 4026, 4061, 4081, 4091, 4095, 4096) -
  c(0, 1, 5, 15, 35, 70, 126, 210, 330, 491, 695, 941, 1225, 1540, 1876, 2220, 2556, 2871, 3155, 3401, 3605, 3766, 3886, 3970, 4026, 4061, 4081, 4091, 4095)
runtimes = data.table(Njobs)[, job_num := .I][, running_time := 3 + 0.0075 * Njobs + 0.5*(Njobs > 100) + 0.0075*Njobs*(Njobs > 200)][]
library(ggplot2)
ggplot(data = runtimes) +
  geom_bar(aes(x = job_num, y = Njobs), stat = 'identity') +
  geom_line(aes(x = job_num, y = running_time / 0.033), color = 'red', size = 2) + 
  scale_y_continuous(
    "Number of Jobs", 
    sec.axis = sec_axis(~ . * 0.033, name = "Guessed Running Time (minutes)")
  )+
  theme(axis.title.x=element_blank(),
        axis.text.x=element_blank(),
        axis.ticks.x=element_blank())

ggsave('DAGmanQueuingSketch.png')