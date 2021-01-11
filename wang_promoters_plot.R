# Ad hoc script : Get and plot harris wang promoter-rbs efficacy across organisms from the paper

# Author: Prashant, # Date: 10/Jan/2021


# preliminaries ----

# Loading libraries, functions and user inputs
source('./general_functions.R') # Source the general_functions file

fylname.dat <- '11,319 set Bs, Ec, Pa (table S2).xlsx'  # file name for data input
fylpth.dat <- 'C:/Users/new/Zotero/storage/SS9XI6IS' # filepath - it's a zotero folder 
fylsheets.dat <- c('BS', "EC", 'PA') # sheets to read

query_for_subset <- 'wang_promoters_query.csv'

# Data import ----

enter.dat <- map_dfr(fylsheets.dat, 
                     ~ read_xlsx(path = str_c(fylpth.dat, fylname.dat, sep = '/'), 
                                 sheet = .x, 
                                 col_types = 'numeric') %>% 
                       mutate(organism = .x)       
)
   
enter.query <- read_csv(query_for_subset)

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
  left_join(exact.query %>% 
              rename('OLIGO ID' = promoter_name)) %>% 
  group_by(organism) %>% 
  arrange(`protein (log10)`)
  # pivot_wider()
  
# plotting ----

slected.dat %>% 
  ggplot(aes(`protein (log10)`,`Plasmid ID`, colour = organism)) + 
  geom_point() +
  geom_line(aes(group = organism))
