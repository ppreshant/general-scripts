# FUnctions to load qPCR data and manipulate it. The functions can be called from another R file

# read in excel file (.xls) of qPCR exported from Quantstudio 3 
  # Make sure to include raw data as well

# calling libraries ; make sure they are installed (install.packages)
library(readxl); library(magrittr); library(tidyverse); library(ggrepel); library(googlesheets4); library(rlang); library(lubridate); library(plotly) 

sheeturls <- list(templates = 'https://docs.google.com/spreadsheets/d/19oRiRcRVS23W3HqRKjhMutJKC2lFOpNK8aNUkC-No-s/edit#gid=478762118',
                  biobot_id = 'https://docs.google.com/spreadsheets/d/1ghb_GjTS4yMFbzb65NskAlm-2Gb5M4SNYi4FHE4YVyI/edit#gid=233791008',
                  sample_registry = 'https://docs.google.com/spreadsheets/d/1mJcCt1wMiOuBic6sRlBZJf8KSNu2y-B5PjzCUu7jPM8/edit#gid=521099478',
                  
                  data_dump = 'https://docs.google.com/spreadsheets/d/1ouk-kCJHERRhOMNP07lXfiC3aGB4wtWXpnYf5-b2CI4/edit#gid=0',
                  raw_ddpcr = 'https://docs.google.com/spreadsheets/d/1jdO_P9SZGezSTLiIARtSmA7qaUuX3wA-jCe7YiQ1sCI/edit#gid=0',
                  complete_data = 'https://docs.google.com/spreadsheets/d/1ltvW7xZf2aUPoBchD4NFGuV0gGWPkdsOW3a0Hxalm-Y/edit#gid=1363292517',
                  wwtp_only_data = 'https://docs.google.com/spreadsheets/d/1dBESjgWSFsOBodFFpYNhWIOAQ2a3wvoLkrn12V_rFck/edit#gid=0' 
)

# reading files and manipulating columns ----

# read in the excel file (from row 36 onwards)
readqpcr <- function(flnm)
{
  fl <- flnm %>%  
    excel_sheets() %>% 
    set_names(.,.) %>% 
    map(read_excel, path = flnm, skip = 38)
  
  # convert CT values into numeric 
  class(fl$Results$CT) <- 'numeric'
  fl
}

columnwise_index <- function(fl)
{ # this orders the samples columnwise in the PCR plate or strip : order will give indices for A1,B1,C1 ... A2,B2,C2... (wherever samples exist)
  fl$Results$`Well Position` %>%  str_sub(2) %>% as.integer() %>% order() # Take the well position A1 etc., extract the number, read as integer (instead of char), order it (1,1,1.. 2,2,... 12,12...) and regurn indices 
}

# lookup table of primer pairs and respective targets - not useful as same information is in `Target` already
primer_table <- c('q1-3' = 'Flipped', 'q4-5' = 'Flipped', 
                  'q5-11' = 'Unflipped', 'q1-2' = 'Unflipped', 'q9-10' = 'Unflipped', 'q4-2' = 'Unflipped', 'q2-4' = 'Unflipped', 'q5-11' = 'Unflipped', 'q6-7' = 'Unflipped',
                  'q12-13' = 'Backbone')

# function to back-calculate CT using standard curve parameters
absolute_backcalc <- function(df, std_par)
{
  target_current <- df$Target %>% unique()
  std_current <- std_par %>% filter(str_detect(target_current, Target))
  
  df %>% mutate(`Copy #` = 10^( (CT - std_current$y_intercept)/std_current$Slope) )
}

read_plate_to_column <- function(data_tibble, val_name)
{ # transforms a plate reader table into a column (named after the top left cell, unless mentioned)
  # eliminates plate row,column numbering ; Select 1 row above the plate (even if it doesn't contain a label)
  
  val_name <- enquo(val_name)
  # colnames(data_tibble) <- data_tibble[1,] # set column names as the first row
  data_tibble[,] %>% gather(key = 'col_num', value = !!val_name, -`<>`) %>% rename(row_num = `<>`) %>% unite('Well Position', c('row_num', 'col_num'), sep = '') %>% drop_na()}

# mutates a subset of data and returns a new array (works for multiple conditions)
mutate_when <- function(data, ...) 
{ # Source: Stackoverflow - https://stackoverflow.com/a/34170176/9049673
  
  dots <- eval(substitute(alist(...)))
  for (i in seq(1, length(dots), by = 2)) 
  {
    condition <- eval(dots[[i]], envir = data)
    mutations <- eval(dots[[i + 1]], envir = data[condition, , drop = FALSE])
    data[condition, names(mutations)] <- mutations
  }
  data
}


# Reads single column of data and converts into a 96 well plate format (Baylor Sample_names)
baylor_col_to_plate <- function(sheetnm)
{
  target <- 'BCoV' # appends to the Sample_name for direct pasting in qPCR template sheet
  flnm <- 'baylor_labels'
  registry_sheet <- 'https://docs.google.com/spreadsheets/d/1mJcCt1wMiOuBic6sRlBZJf8KSNu2y-B5PjzCUu7jPM8/edit#gid=1011717808'
  
  # bb_sheet <- c('Week 13 (7/6)')
  # 
  # # week_match <- flnm %>% str_extract('[:digit:]+(?=-|_)') 
  # 
  # 
  # biobot_url <- 'https://docs.google.com/spreadsheets/d/1ghb_GjTS4yMFbzb65NskAlm-2Gb5M4SNYi4FHE4YVyI/edit#gid=233791008' 
  # 
  # # Getting biobot names: named vector for backconversion
  # biobot_translator <- read_sheet(biobot_url, sheet = bb_sheet) %>% 
  #   rename('Biobot ID' = matches('Biobot|Comments', ignore.case = T), 
  #          'WWTP' = contains('SYMBOL', ignore.case = T), 
  #          'FACILITY NAME' = matches('FACILITY NAME', ignore.case = T)) %>%
  #   
  #   drop_na(WWTP) %>% 
  #   mutate('biobot_baylor' = str_replace(`Biobot ID`,'\\.', '_'), WWTP = as.character(WWTP)) %>%  
  #   mutate(bb_translator = set_names(biobot_baylor , WWTP)) %>% 
  #   pull(bb_translator)
  
  
  # Read sample sheet
  
  baylor_names <- str_c('excel files/Baylor/', flnm, '.xlsx') %>% 
    read_xlsx (sheet = sheetnm) %>% 
    rename('Sample' = matches('Site')) %>% 
    select('Well', 'Sample') %>% 
    
    mutate('row' = str_match(Well, '[:upper:]'), 'col' = str_match(Well, '[:digit:]+')) %>% 
    select(-Well) %>% 
    
    # correct names - put a dot before replicate number
    mutate_at('Sample', ~str_match(., '(^[:upper:]+|^[:digit:]+).*([:digit:])') %>% {str_c(.[,2], .[,3], sep = '.')}) %>%
    
    
    # substitute biobot ids in
    # mutate_at('Sample', str_replace_all , biobot_translator) %>% 
    
    # attach template name and week ID
    mutate_at('Sample', ~str_c(target, '-', str_extract(sheetnm, '[:digit:]*(?=-[:digit:])'), '_' , .))
  
  # Convert column to table - 96 well
  
  baylor_table <- baylor_names %>% 
    pivot_wider(names_from = 'col', values_from = 'Sample') %>% 
    rename('<>' = row) %>% 
    mutate('-' = '', .before = 1)
  
  View(baylor_table) # manually copy paste wherever
  
  # write to sample registry 96 well sheet
  sheet_append(registry_sheet, sheet = '96 well RNA', tibble(c('', '', ''))) # append 3 empty rows
  sheet_append(registry_sheet, sheet = '96 well RNA', baylor_table)
  
}

# mutates a subset of data and returns a new array (does multiple mutations on same condition)
mutate_cond <- function(.data, condition, ..., envir = parent.frame()) 
{ # Source: Stackoverflow -  https://stackoverflow.com/a/34096575/9049673
  
  condition <- eval(substitute(condition), .data, envir)
  .data[condition, ] <- .data[condition, ] %>% mutate(...)
  .data
}

# Gets the 96 well layout with template names matching the experiment ID from filename in a google sheet
get_template_for <- function(bait, sheet_url = sheeturls$templates)
{ # Looking for WWx or Stdx - example WW21 or Std7 within the filename; Assumes plate spans from row B to N (1 row below the matching ID)
 
   # Finding the plate to be read
  plate_names_row <- read_sheet(sheet_url, sheet = 'Plate layouts', range = 'C:C', col_types = 'c')
  m_row <- plate_names_row %>% unlist() %>% as.character() %>% 
    # find the row with standard beginings matching the filename
    str_detect(., str_c('^', bait %>% str_match('^(WW|Std|dd.WW)[:alnum:]*') %>% .[1]) ) %>% 
    which() + 1
  range_to_get <- str_c('B', m_row + 1, ':N', m_row + 9)
  
  # Eror message and terminate if plate ID is not unique
  if(length(m_row) > 1) stop( str_c('Plate ID of :', bait, 'repeats in', paste0(m_row, collapse = ' & '), 'row numbers. Please fix and re-run the script', sep = ' '))
  
  # read the template corresponding to the file name
  plate_template_raw <- read_sheet(sheet_url, sheet = 'Plate layouts', range = range_to_get)
  
  # Convert the 96 well into a single column, alongside the Well position
  plate_template <- read_plate_to_column(plate_template_raw, 'Sample_name') # convert plate template (Sample_names) into a single vector, columnwise
  
}

# Check for neighboring dates and merge them
harmonize_week <- function(week_cols) 
{
  
  # Pick numeric entries in column (the rest will be restored as is)
  num_week <- week_cols %>% str_extract('[:digit:]{3,4}') %>%  as.numeric() %>% unique() %>% .[!is.na(.)]
  
  # Check for consecutive dates
  repl_week <- num_week %>% 
    str_c() %>% 
    str_replace('([:digit:]+)([:digit:]{2})', '\\1/\\2/20') %>% 
    mdy() %>%  # convert to dates
    map_dbl(., function(x) num_week[. %in% c(x, x-1)] %>% min()) %>% # make the min entry of consecutive dates
    as.character() %>% 
    set_names(nm = num_week)
  
  if(repl_week %>% is_empty() %>% {!.}) new_week_cols <- str_replace_all(week_cols, repl_week)
  else week_cols
  
}

# Data writing output ----

# This function writes to the specified google sheet if the current sheet does
# not exist. If the sheet does exist it will ask the user before writing.
check_ok_and_write <- function(data, sheet_URL, title_name)
{
  write_ok <- TRUE
  sheet_dne <- FALSE
  
  # this block sets sheet_dne to true if the sheet does not exist
  sheet_dne <- tryCatch(
    expr = {
      read_sheet(sheet_URL, sheet = title_name)
      message("Sheet already exists")
      FALSE
    },
    error = function(e){
      message('Sheet does not exist')
      return(TRUE)
      print(e)
    }
  )
  
  # if the sheet exists (sheet_dne is false), then ask the user if
  # they want to overwrite. If the user selects cancel, then abort
  if (!sheet_dne) {
    # write_ok <- askYesNo(paste("A sheet with the name", title_name, "already exists. Do you want to overwrite?", sep=" "))
    write_ok <- menu(c('Yes', 'No'), title = paste("A sheet with the name", title_name, "already exists. Do you want to overwrite?", sep=" "))
    if (write_ok == 2){
      stop("Cancel selected, script aborted.")
    }
  }
  
  if (write_ok == 1) {
    write_sheet(data, sheet_URL, sheet=title_name)
  }

}

# standard curve and regressions ----

# Plot Standard curve
plotstdcurve <- function(results_qpcr, plttitle, xlabel)
{
  plt <- results_qpcr %>% 
ggplot(.) + aes(x = log10(Quantity), y = CT, color = `Target`) + geom_point() +
theme_classic() + scale_color_brewer(palette="Set1") + theme(plot.title = element_text(hjust = 0.5)) + 
ggtitle(plttitle) + xlab(xlabel) + ylab(expression(C[q])) +
stat_smooth(data = filter(results_qpcr, Task == 'STANDARD'), method ="lm", se = F) # plots linear regression line
}

# getting regression line and R2 values to put into the standard curve plot 
# Source: https://stackoverflow.com/a/7549819/9049673

# least square fitting - linear regression for standard curve
lm_std_curve <- function(df, trig = 0)
{
  x = df %>% pull(Quantity) %>% log10 ; y = df %>% pull(CT)
  m <- lm(y ~ x, df)
  lm_eqn(m, trig)
}


# linear regression equation
lm_eqn <- function(m, trig = 0){
  
  eq <- substitute(italic(y) == b %.% italic(x)+ a*","~~italic(r)^2~":"~r2, 
                   list(a = format(unname(coef(m)[1]), digits = 4), 
                        b = format(unname(coef(m)[2]), digits = 3), 
                        r2 = format(summary(m)$r.squared, digits = 3)))
  # if(trig == 'coeff') c(round(coef(m)[2], 2), round(summary(m)$r.squared, 2))
  if(trig == 'coeff') tibble(slope = round(coef(m)[2], 2), y_intercept = round(coef(m)[1], 2), r_square = round(summary(m)$r.squared, 3))
  else as.character(as.expression(eq)); 
}

  
  optional1 <- function()
    {# output the difference between consecutive CT values
    tsumrev <- trev %>% group_by(`Sample_name`) %>% summarise(CT = mean(CT), Quantity = mean(Quantity), CT_sd = sd(CT))
    diff(tsumrev$CT) %>% round(2)}

# Tm plots ----
# plotting functions for Melting temperature

# plots all the Tm's if samples have multiple peaks in the melting curve
plotalltms <- function(results_relevant)
{ 
  # Gather the Tm's into another data frame and merge into 1 column
  tmfl <- results_relevant %>% select(`Sample_name`, starts_with('Tm')) %>% gather('Peak number','Tm',-`Sample_name`)
  
  # plot the Tm ; Graph will now show
  plttm2 <- tmfl %>% ggplot(.) + aes(x = `Sample_name`, y = Tm) + geom_point(aes(color = `Peak number`), size = 2) +
    theme_classic() + scale_color_brewer(palette="Set1") + 
    theme(plot.title = element_text(hjust = 0.5),axis.text.x = element_text(angle = 90, hjust = 1, vjust = .3)) + 
    ggtitle(paste(title_name,': Melting curves')) + facet_grid(~`Primer pair`, scales = 'free_x', space = 'free_x')
}

# plot the first Tm only ; Graph will now show
plottm1 <- function(results_relevant)
{ 
  plttm <- results_relevant %>% ggplot(.) + aes(x = `Sample_name`, y = Tm1) + geom_point(color = 'red', size = 2) +
    theme_classic() + scale_color_brewer(palette="Set1") + 
    theme(plot.title = element_text(hjust = 0.5),axis.text.x = element_text(angle = 90, hjust = 1, vjust = .3)) + 
    ggtitle(paste(title_name,': Melting curves')) + facet_grid(~`Primer pair`, scales = 'free_x', space = 'free_x')
}


# Plotting functions ----


# Plotting mean, sd and individual replicates jitter
plot_mean_sd_jitter <- function(.data_list = long_processed_minimal, 
                                long_format = TRUE, 
                                
                                measure_var = 'Copy #', 
                                colour_var = Target, x_var = assay_variable, y_var = `Copy #`, 
                                facet_var = `Sample_name`, 
                                
                                sample_var = '.*', exclude_sample = F, 
                                WWTP_var = wwtp_manhole_names, exclude_WWTP = F, 
                                target_filter = '.*', 
                                
                                ascending_order = FALSE, 
                                title_text = title_name, ylabel = 'Genome copies/ul RNA', xlabel = plot_assay_variable, 
                                facet_style = 'grid')

{ # Convenient handle for repetitive plotting in the same format; Specify data format: long vs wide (specify in long_format = TRUE or FALSE)
  
  .dat_filtered <- .data_list %>% map( ~ filter(.x, 
                                                if('Sample_name' %in% colnames(.x)) str_detect(`Sample_name`, sample_var, negate = exclude_sample) else TRUE, 
                                                if('WWTP' %in% colnames(.x)) str_detect(WWTP, WWTP_var, negate = exclude_WWTP) else TRUE, 
                                                str_detect(Target, target_filter))
  )
  
  # filtering data to be plotted by user inputs
  if(long_format) # use long format if not plotting Copy #s - ex. Recovery, % recovery etc.
  {
    
    .data_to_plot <- .dat_filtered %>% map(filter,
                                           Measurement == measure_var)
    if(ascending_order) .data_to_plot$summ.dat %<>% mutate_at('WWTP', as.character) %>% 
      arrange(`mean`) %>% 
      mutate_at('WWTP', as_factor)
    
    y_var <- sym('value') # default y variable is value
    
    summ_actual_spike_in <- .dat_filtered$summ.dat %>% filter(str_detect(Measurement,'Actual'))
    
  } else
    
  {
    .data_to_plot <- .dat_filtered 
    
    if(ascending_order) .data_to_plot$summ.dat %<>% mutate_at('WWTP', as.character) %>% 
      arrange(`mean`) %>% 
      mutate_at('WWTP', as_factor)
    
  }
  
  # Exit with a useful message if data is empty
  if(.data_to_plot %>% map_lgl(plyr::empty) %>% any()) return('Data does not exist')  
  
  # plotting
  plt1 <- .data_to_plot$summ.dat %>% ggplot(aes(x = {{x_var}}, y = mean, colour = {{colour_var}})) +
    geom_point(size = 2) + geom_errorbar(aes(ymin = mean - sd, ymax = mean + sd), width = .1) +
    
    # Individual data points
    geom_jitter(data = .data_to_plot$raw.dat, aes(y = {{y_var}}, alpha = map_chr({{y_var}}, ~. == 0), size = map_chr({{y_var}}, ~. == 0)), width = .2, show.legend = F ) +
    scale_alpha_manual(values = c(.3, 1)) + scale_size_manual(values = c(1, 2)) + # manual scale for emphasizing unamplified samples
    
    # Plotting actual spike ins (only for Recovery plot'; only with long format data )
    { if(measure_var == 'Recovered') list(geom_point(data = summ_actual_spike_in, colour = 'black', shape = 21), 
             geom_line(data = summ_actual_spike_in, aes(group = {{colour_var}})))
    } +
    
    # Facetting
    facet_grid(cols = vars({{facet_var}}), scales = 'free_x', space = 'free_x') +
    
    # experimental - conditional facetting (doesn't work for unknown reasons) : Just facet_grid the output to remove facets or add new!
    # { if (facet_style == 'grid') list(facet_grid(cols = vars({{facet_var}}), scales = 'free_x', space = 'free_x'))
    #   if (facet_style == 'wrap free') list(facet_wrap(facets =  vars({{facet_var}}), scales = 'free')) 
    #   else NULL
    # } +
    
    # Labelling
    ggtitle(title_text) + ylab(ylabel) #+ xlab(xlabel)

  plt1.formatted <- plt1 %>% format_classic() # clean formatting
  
}

# plotting individual replicates
plot_biological_replicates <- function(results_abs, title_text = title_name, xlabel = plot_assay_variable)
{ # Convenient handle for repetitive plotting 'Copy #' vs biological replicate
  
  plt <- results_abs %>% ggplot(aes(x = `Tube ID`, y = `Copy #`, color = Target)) + ylab('Copies/ul RNA extract') +    # Specify the plotting variables 
    geom_point(size = 2) + facet_grid(~`Sample_name`, scales = 'free_x', space = 'free_x') + # plot points and facetting
    ggtitle(title_text) + xlab(xlabel)
  plt.formatted <- plt %>% format_classic(.) %>% format_logscale_y() # formatting plot, axes labels, title and logcale plotting
}

# Scatter plot with a linear regression fit and equation
# Almost generalized: Only need to make it work with grouping_var - accounting for NULL and presence of a variable...

plot_scatter <- function(.data = processed_quant_data, 
                         x_var = N1_multiplex, y_var = N2_multiplex, colour_var = NULL, shape_var = NULL,
                         grouping_var = NULL, # CURRENTLY ONLY WORKS FOR NULL GROUP
                         already_pivoted_data = 'no',
                         title_text = title_name,
                         
                         show_y.equals.x_line = 'yes',
                         text_for_equation = 'Rsquare', # choice: "Rsquare" or "full equation"
                         measure_var = 'Copy #',
                         text_cols = minimal_label_columns,
                         sample_var = str_c(extra_categories, '|NTC|Vaccine'), 
                         exclude_sample = T)
  
{ # Convenient handle for repetitive plotting in the same format; 
  
  
  # Preliminaries
  
  
  # convert column names (target names) into string
  checkx <- strx <- paste(substitute(x_var)) 
  checky <- stry <- paste(substitute(y_var))
  
  modx <- function(x) x # unimplemented feature: modifier in case this function gets transformed variables
  
  if(already_pivoted_data == 'no')
  {  # filtering data for plotting according to function inputs
    .data_for_plot <- .data %>% 
      select(all_of(text_cols), biological_replicates, all_of(measure_var)) %>% 
      filter(str_detect(`Sample_name`, sample_var, negate = exclude_sample)) %>% # account for missing var
      pivot_wider(names_from = 'Target', values_from = all_of(measure_var)) %>% # Target can be generalized?
      ungroup() # why did you ungroup - for the lm ..?
  } else .data_for_plot <- .data #%>%  # direct carrying of data to next steps
  # {if(!is.null(grouping_var)) group_by(., {{grouping_var}}) else . } # DISABLED FOR CHECKING IF PLOTLY RUNS (GROUPS NOT WORKING RIGHT NOW
  
  
  
  # error handling
  
  # If plotting transformations of variables : Not fully implemented yet
  if(enexpr(x_var) %>% is.call()) {checkx <-  enexpr(x_var)[2] %>% paste(); }# modx <- eval(enexpr(x_var)[1])}
  if(enexpr(y_var) %>% is.call()) checky <-  enexpr(y_var)[2] %>% paste()
  
  # If X/Y vars are not present in the data, stop running and print a message
  if(.data_for_plot %>% names() %>% 
     {checkx %in% . & checky %in% . } %>% 
     !.) return('X or Y axis values for this scatterplot are not present in the data frame provided. 
Check if x_var and y_var are present in .data')
  
  # If duplicates of data exist in the data, stop running and print a message
  if((.data_for_plot %>% select(all_of(c(checkx, checky))) %>% 
      map_lgl(~class(.) == 'numeric') %>% 
      sum()) < 2) 
  { duplicated_data_points <- .data_for_plot %>% 
    filter(map({{x_var}}, length) > 1)
  print(duplicated_data_points)
  
  return('Repeated data instances for the same WWTP found in this scatterplot')
  }
  
  
  
  # For linear regression data
  
  # Making linear regression formula (source: https://stackoverflow.com/a/50054285/9049673)
  fmla <- new_formula(enexpr(y_var), enexpr(x_var))
  
  # Max and ranges for plotting
  xyeq <- .data_for_plot %>% summarise(across(where(is.numeric), max, na.rm = T)) %>% select(all_of(c(checkx, checky))) %>% min() %>% {.*0.9} %>% modx()
  
  # linear regression equation
  lin_reg_eqn <- .data_for_plot %>% mutate(across(all_of(c(checkx, checky)), ~if_else(.x == 0, NaN, .x))) %>% 
    lm(fmla, data = ., na.action = na.exclude) %>% lm_eqn(., trig = text_for_equation)
  
  
  
  
  # plotting part
  
  plt1 <- .data_for_plot %>% 
    ggplot(aes(x = {{x_var}}, y =  {{y_var}} )) +
    geom_point(size = 2, mapping = aes(colour = {{colour_var}}, shape = {{shape_var}})) +
    
    # linear regression
    geom_smooth(method = 'lm') + # .. , mapping = aes(group = {{grouping_var}})  # DISABLED GROUPS 
    geom_text(data = . %>% summarise(across(where(is.numeric), max, na.rm = T) ),
              # mapping = aes(group = {{grouping_var}}), # DISABLED FOR CHECKING IF PLOTLY RUNS (GROUPS NOT WORKING RIGHT NOW)
              label = lin_reg_eqn, parse = TRUE, show.legend = F, hjust = 'inward', nudge_x = -5) +
    
    # Dummy y = x line
    {if(show_y.equals.x_line == 'yes') list(geom_abline(slope = 1, intercept = 0, alpha = .4),
                                            annotate(geom = 'text', x = xyeq, y = xyeq, label = 'y = x', alpha = .3))
      } +
    
    # Labeling
    ggtitle(title_text, subtitle = measure_var)
}


# plot formatting ---- 
 

# Set theme universally : format as classic, colours = Set1
theme_set(theme_classic()) # theme
scale_colour_discrete <- function(...) { # palette
  scale_colour_brewer(..., palette="Set1")
}

# plot formatting function : format as classic, colours = Set1
format_classic <- function(plt)
{ # formats plot as classic, with colour palette Set1, centred title, angled x axis labels
  plt <- plt +
    theme_classic() + scale_color_brewer(palette="Set1")
}
 
  # plot formatting function : format as classic, colours = Set1
  format_classic <- function(plt)
  { # formats plot as classic, with colour palette Set1, centred title, angled x axis labels
    plt <- plt +
      theme_classic() + scale_color_brewer(palette="Set1")
  }
  
  # plot formatting function : format as logscale
format_logscale_y <- function(plt)
  { # extra comments
    plt <- plt +
      scale_y_log10(  # logscale for y axis with tick marks
        labels = fancy_scientific
        #labels = scales::trans_format("log10", scales::math_format(10^.x) )
      )
  }

# plot formatting function : format as logscale x
format_logscale_x <- function(plt)
{ # extra comments
  plt <- plt +
    scale_x_log10(  # logscale for y axis with tick marks
      labels = fancy_scientific
      #labels = scales::trans_format("log10", scales::math_format(10^.x) )
    )
}

  
  # formatting labels in logscale cleanly : a x 10^b
  fancy_scientific <- function(l) {
    # turn in to character string in scientific notation
    l <- format(l, scientific = TRUE)
    # quote the part before the exponent to keep all the digits
    l <- gsub("^(.*)e", "'\\1'e", l)
    # remove + after exponent, if exists. E.g.: (3x10^+2 -> 3x10^2) 
    l <- gsub("e\\+","e",l)
    # turn the 'e+' into plotmath format
    l <- gsub("e", "%*%10^", l)
    # convert 1x10^ or 1.000x10^ -> 10^ 
    l <- gsub("\\'1[\\.0]*\\'\\%\\*\\%", "", l)
    # return this as an expression
    parse(text=l)
  }
  
  # use as ggplot(df,aes(x,y)) + geom_point() + scale_y_log10(labels = fancy_scientific)
