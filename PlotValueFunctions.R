library(ncdf4)
library(foreach)
library(iterators)
library(ggplot2)
library(cowplot)
library(RColorBrewer)


files = data.table(filename = list.files('FromCHTC/OilGameSmallerRho/'),
                   path = list.files('FromCHTC/OilGameSmallerRho/', full.names = TRUE))
files[, firstbit := tstrsplit(filename, '_')[1]]
files = files[firstbit == 'solarray']
files[, no_extension := gsub('\\.nc', '', filename)]
files[, Nplayers := lapply(tstrsplit(no_extension, '_')[2], as.integer)]
files[, paste0('DS_', 1:max(files$Nplayers)) := lapply(tstrsplit(no_extension, '_')[3:(2 + max(files$Nplayers))], as.integer)]

firstpass = files[(DS_1 == 10 & DS_3 == 10 & DS_4 == 10) |
                    (DS_1 == 10 & DS_3 == 10 & DS_4 == 5) |
                    (DS_1 == 10 & DS_3 == 5 & DS_4 == 10) |
                    (DS_1 == 5 & DS_3 == 10 & DS_4 == 10) |
                    (DS_1 == 10 & DS_3 == 5 & DS_4 == 5) |
                    (DS_1 == 5 & DS_3 == 5 & DS_4 == 5) |
                    (DS_1 == 1 & DS_3 == 1 & DS_4 == 1)]

data = 
  foreach(row = iter(firstpass, by = 'row'), .combine = rbind) %do% {
    solvals = data.table(t(ncvar_get(nc_open(row$path))))
    setnames(solvals, c(t(outer(c('DS_', 'control_', 'value_'), 1:row$Nplayers, FUN = 'paste0'))))
  }

to_plot = data[(DS_1 == 100 & DS_3 == 100 & DS_4 == 100) |
                    (DS_1 == 100 & DS_3 == 100 & DS_4 == 50) |
                    (DS_1 == 100 & DS_3 == 50 & DS_4 == 100) |
                    (DS_1 == 50 & DS_3 == 100 & DS_4 == 100) |
                    (DS_1 == 100 & DS_3 == 50 & DS_4 == 50) |
                    (DS_1 == 50 & DS_3 == 50 & DS_4 == 50) |
                    (DS_1 == 1 & DS_3 == 1 & DS_4 == 1)]

to_plot[, paste0('oil_remaining_', 1:max(firstpass$Nplayers)) := (100 - .SD) / 10, .SDcols = c('DS_1', 'DS_2', 'DS_3', 'DS_4')]
to_plot[, `Other Firms'\nOil Remaining` := paste0('(', paste(oil_remaining_1, oil_remaining_3, oil_remaining_4, sep = ', '), ')')]
to_plot[, `Oil Remaining` := oil_remaining_2]
to_plot[, `Value` := value_2]

ggplot(to_plot) +
  geom_line(aes(x = `Oil Remaining`, `Value`, color = `Other Firms'\nOil Remaining`), size = 1.3) +
  ggtitle('Value') + # (ρ = 0.02)') +
  theme_cowplot() +
  scale_color_brewer(palette = "Dark2") +
  theme(plot.title = element_text(hjust = 0)) +
  scale_y_continuous(expand = c(0, 0), limits = c(0, 6.35)) +
  scale_x_continuous(expand = c(0, 0), limits = c(0, 10)) +
  theme(axis.title.y = element_blank()) # +
  # theme(legend.position = c(0.8, 0.8))

ggsave('ValueFunctionsOilGameRho02.png')


# to_plot[, `Control` := control_2]
to_plot[oil_remaining_1 > (control_1 + control_2 + control_3 + control_4), `Planner` := 0]
to_plot[oil_remaining_1 < (control_1 + control_2 + control_3 + control_4) & oil_remaining_1 > 0, `Planner` := (control_1 + control_2 + control_3 + control_4) - oil_remaining_1]
to_plot[oil_remaining_1 == 0, `Planner` := (control_1 + control_2 + control_3 + control_4)]

to_plot2 = rbind(cbind(to_plot, data.table(Production = to_plot$control_2, Solution = "Equilibrium")),
                 cbind(to_plot, data.table(Production = to_plot$Planner, Solution = "Planner"))) # , to_plot[][, `:=`(Production = Planner, Solution = "Planner")])

to_plot2[, make_transparent := ifelse(DS_1 == 100 & DS_3 == 100 & DS_4 == 100 & Solution == "Equilibrium", "yes", "no")]

ggplot(to_plot2) +
  geom_line(aes(x = `Oil Remaining`, `Production`, color = `Other Firms'\nOil Remaining`, linetype = Solution, alpha = make_transparent), size = 1.3) +
  ggtitle('Output') +
  theme_cowplot() +
  scale_color_brewer(palette = "Dark2") +
  scale_alpha_manual(values = c(1, 0.5), guide = "none") +
  theme(plot.title = element_text(hjust = 0)) +
  # scale_y_continuous(expand = c(0, 0), limits = c(0, 6.35)) +
  scale_x_continuous(expand = c(0, 0), limits = c(0, 10)) +
  theme(axis.title.y = element_blank()) # +
  # theme(legend.position = c(0.8, 0.8))

ggsave('PolicyFunctionsOilGameRho02.png')

