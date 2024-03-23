# SS_nested_facets.R

# also called nested categorical structures?

# prereqs ----

library(tidyverse)


# example data ----
a <- tibble(a1 = rep(1:3, 10), # repeat same data 10 times
            a2 = a1 ^ 2 + rnorm(30, sd = 3), # add some random noise
            
            # make two sets of vectors for facets 
            inner = rep(c(rep(letters[1:3], 3), 'x'), each = 3),
            outer = rep(c(rep(LETTERS[1:3], each = 3), 'D'), each = 3)
)


# installing 
# pak::pkg_install('ggh4x')

# plot ----

ggplot(a, aes(a1, a2)) + 
  geom_point() + geom_line() +
  
  # facet with nesting
  ggh4x::facet_nested_wrap(facets = vars(inner, outer), nrow = 1) + 
  
  # this theme removes grey background and makes boxes around facet labels
  theme_classic()
            