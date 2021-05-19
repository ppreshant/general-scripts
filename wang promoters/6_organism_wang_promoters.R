# Ad hoc script : Get and plot harris wang promoter-rbs efficacy across organisms from the paper

# Author: Prashant, # Date: 19/May/2021


# preliminaries ----

# Loading libraries, functions and user inputs
source('./general_functions.R') # Source the general_functions file

fylname.dat <- '6 organisms, Tx and Tl data.xlsx'  # file name for data input
fylpth.dat <- 'C:/Users/new/Zotero/storage/HVDSQTYF' # filepath - it's a zotero folder 
fylsheets.dat <- c('Log2 Transcription', 'Log10 Translation') # sheets to read

general.path <- 'wang promoters/'
query_for_subset <- 'wang_promoters_query.csv'

# Data import ----

enter.dat <- map_dfr(fylsheets.dat, 
                     ~ readxl::read_xlsx(path = str_c(fylpth.dat, fylname.dat, sep = '/'), 
                                         sheet = .x,
                                         range = 'A1:G242',
                                         col_types = 'numeric') %>% 
                       mutate(Measurement = .x)       
)

enter.query <- paste0(general.path, query_for_subset) %>% read_csv

# wrangling ----

exact.query <- enter.query %>% 
  mutate(across(promoter_name, 
                ~ str_remove(., 'HW') %>%  # remove letters from the prmototer name
                  as.numeric)  # convert to numbers to match to the raw dataset
  )
rm(enter.query) # removing original query - cleanup clutter

sbset.dat <- enter.dat %>% 
  filter(id %in% exact.query$promoter_name)
rm(enter.dat) # removing raw data from pipeline - too large, slowing down R 

slected.dat <- sbset.dat %>% 
  # select(id, tx_raw, `protein (log10)`, organism) %>%  # select only relevant columns to plot
  rename(promoter_name = 'id') %>% 
  left_join(exact.query) %>% # join query plasmid names here
  pivot_longer(cols = -c(promoter_name, `Plasmid ID`, Measurement), names_to = 'organism') %>% 
  group_by(organism) %>% # perform operations within each measurement
  arrange(Measurement, organism, value) %>% 
  mutate(across(c('Plasmid ID', 'promoter_name'), ~ as.character(.) %>%  fct_inorder))
# pivot_wider()


# mid.all.dat <- enter.dat %>% 
#   select(id, tx_raw, `protein (log10)`, organism) %>%  # select only relevant columns to plot
#   rename(promoter_name = 'OLIGO ID') %>% 
#   left_join(exact.query) %>% 
#   group_by(organism) %>% 
#   arrange(organism, `protein (log10)`) %>% 
#   mutate(across(c('promoter_name'), ~ as.character(.) %>%  fct_inorder))


# plotting ----

# plasmid id vs protein
plt_all <- slected.dat %>% 
  {ggplot(.,aes(value,`Plasmid ID`, colour = organism)) + 
      geom_point() +
      geom_line(aes(group = organism), orientation = 'y') + 
      facet_grid(~ Measurement, scales = 'free_x') +
      scale_color_brewer(palette = 'Dark2') +
      ggtitle('Broad host expression', subtitle = 'Selected harris wang promoters')} %>% 
  print()

ggsave(str_c(general.path, '6_organism_existing.png'), width = 8, height = 4.6 )


# HW number vs protein
plt_all_HWnumber <- slected.dat %>% 
  {ggplot(.,aes(value, promoter_name, colour = organism)) + 
      geom_point() +
      geom_line(aes(group = organism), orientation = 'y') + 
      facet_grid(~ Measurement, scales = 'free_x') +
      scale_color_brewer(palette = 'Dark2') +
      ggtitle('Broad host expression', subtitle = 'Selected harris wang promoters')} %>% 
  print()

ggsave(str_c(general.path, '6_organism_existing_HW.png'), width = 8, height = 4.6 )


# interactive plots
# ggplotly(plt_tx, dynamicTicks = T)

# plotting all data
# # HW number vs protein
# all.plt_prot_HW <- ggplot(mid.all.dat, aes(`protein (log10)`, promoter_name, colour = organism)) + 
#   geom_point(alpha = 0.2) +
#   # geom_line(aes(group = organism), orientation = 'y') + 
#   ggtitle('Translation', subtitle = 'Selected harris wang promoters')
# ggsave(str_c(general.path, 'protein plot_all HW.png'), plot = all.plt_prot_HW, width = 4, height = 4.6)
# 
# 
# # plasmid id vs transcription 
# all.plt_tx <- mid.all.dat %>% 
#   ggplot(aes(tx_raw, promoter_name, colour = organism)) + 
#   geom_point(alpha = 0.2) +
#   # geom_line(aes(group = organism), orientation = 'y') + 
#   ggtitle('Transcription', subtitle = 'Selected harris wang promoters')
# ggsave(str_c(general.path, 'Transcription plot_all HW.png'), plot = all.plt_tx, width = 4, height = 4.6)

