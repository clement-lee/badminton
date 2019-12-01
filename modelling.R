### 00) prelim
library(tidyverse)
library(magrittr)
library(Rcpp)
sourceCpp("algo.cpp")



### 01) read-in data
df0.matches <- read_csv(
    "data.csv",
    col_types = c(
        player1 = col_character(),
        player2 = col_character(),
        games1 = col_integer(),
        games2 = col_integer(),
        year = col_integer()
    )
)
n <- df0.matches %>% nrow()

df0.players <- bind_rows(
    df0.matches %>% select(player = player1),
    df0.matches %>% select(player = player2)
) %>%
    count(player) %>%
    mutate(id = seq_along(player) %>% subtract(1L))
m <- df0.players %>% nrow()

df1.matches <- df0.matches %>%
    left_join(df0.players %>% select(player1 = player, id1 = id), "player1") %>%
    left_join(df0.players %>% select(player2 = player, id2 = id), "player2")



### 02)
