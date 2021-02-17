# Ad hoc script : Get and plot harris wang promoter-rbs efficacy across organisms from the paper

# Author: Prashant, # Date: 10/Jan/2021


# preliminaries ----

# Loading libraries, functions and user inputs
source('./general_functions.R') # Source the general_functions file

fylname.dat <- '11,319 set Bs, Ec, Pa (table S2).xlsx'  # file name for data input
fylpth.dat <- 'C:/Users/new/Zotero/storage/SS9XI6IS' # filepath - it's a zotero folder 
fylsheets.dat <- c('BS', "EC", 'PA') # sheets to read

general.path <- 'wang promoters/'
query_for_subset <- 'wang_promoters_query.csv'

# Data import ----

enter.dat <- map_dfr(fylsheets.dat, 
                     ~ readxl::read_xlsx(path = str_c(fylpth.dat, fylname.dat, sep = '/'), 
                                 sheet = .x, 
                                 col_types = 'numeric') %>% 
                       mutate(organism = .x)       
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
  filter(`OLIGO ID` %in% exact.query$promoter_name)
rm(enter.dat) # removing raw data from pipeline - too large, slowing down R 

slected.dat <- sbset.dat %>% 
  select(`OLIGO ID`, tx_raw, `protein (log10)`, organism) %>%  # select only relevant columns to plot
  rename(promoter_name = 'OLIGO ID') %>% 
  left_join(exact.query) %>% 
  group_by(organism) %>% 
  arrange(organism, `protein (log10)`) %>% 
  mutate(across(c('Plasmid ID', 'promoter_name'), ~ as.character(.) %>%  fct_inorder))
  # pivot_wider()
  
# plotting ----

# plasmid id vs protein
plt_prot <- slected.dat %>% 
  {ggplot(.,aes(`protein (log10)`,`Plasmid ID`, colour = organism)) + 
  geom_point() +
  geom_line(aes(group = organism), orientation = 'y') + 
  ggtitle('Translation', subtitle = 'Selected harris wang promoters')} %>% 
  print()

ggsave(str_c(general.path, 'protein plot.png'), width = 4, height = 4.6 )


# HW number vs protein
plt_prot_HW <- ggplot(slected.dat, aes(`protein (log10)`, promoter_name, colour = organism)) + 
  geom_point() +
  geom_line(aes(group = organism), orientation = 'y') + 
  ggtitle('Translation', subtitle = 'Selected harris wang promoters')
ggsave(str_c(general.path, 'protein plot with HW number.png'), plot = plt_prot_HW, width = 4, height = 4.6)


# plasmid id vs transcription 
plt_tx <- slected.dat %>% 
  ggplot(aes(tx_raw,`Plasmid ID`, colour = organism)) + 
  geom_point() +
  geom_line(aes(group = organism), orientation = 'y') + 
  ggtitle('Transcription', subtitle = 'Selected harris wang promoters')
ggsave(str_c(general.path, 'transcription plot.png'))

# combined plot
# combplt <- plt_prot + plt_tx

# interactive plots
# ggplotly(plt_tx, dynamicTicks = T)
