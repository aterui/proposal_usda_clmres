
# setup -------------------------------------------------------------------

source("code/library.R")


# co-author list ----------------------------------------------------------

## your google id - can be found in the link bar of your google scholar page
sid <- "H9OuCKsAAAAJ"

## get publication list
pubs <- get_publications(sid) %>% 
  as_tibble()

## get coauthor list
coauthor <- pubs %>% 
  rowwise() %>%
  summarise(year = year,
            authors = get_complete_authors(sid,
                                           pubid,
                                           initials = FALSE))
## one author per row
df_coauthor <- coauthor %>% 
  rowwise() %>% 
  reframe(year = year,
          author = str_split(authors,
                             pattern = ",\\s") %>% 
            flatten_chr()) %>% 
  filter(year >= 2021,
         author != "Akira Terui") %>% 
  distinct(author) %>% 
  separate(author,
           into = c("a", "b", "last"),
           sep = "\\s",
           fill = "left") %>% 
  mutate(first = ifelse(is.na(a), b, paste(a, b)),
         full = paste(last, first, sep = ", ")) %>% 
  select(-a, -b) %>% 
  arrange(last)

write_csv(df_coauthor, "output/coauthor_terui.csv")
