### 00) prelim
print(Sys.time())
library(tidyverse)
library(magrittr)
library(coda)
library(glue)
library(Rcpp)
sourceCpp("algo.cpp")





### 01) read-in data
## matches
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

## players
df0.players <- bind_rows(
    df0.matches %>% select(player = player1),
    df0.matches %>% select(player = player2)
) %>%
    count(player)
min.matches <- 3L #### can change
df1.players <- bind_rows(
    df0.players %>%
    filter(n >= min.matches) %>%
    mutate(id = seq_along(player)),
    df0.players %>%
    filter(n < min.matches) %>%
    mutate(id = 0L)
)
m <- df1.players$id %>% unique() %>% length()

## matches again
df1.matches <- df0.matches %>%
    left_join(df1.players %>% select(player1 = player, id1 = id), "player1") %>%
    left_join(df1.players %>% select(player2 = player, id2 = id), "player2")
y1 <- df1.matches$games1
y2 <- df1.matches$games2
x1 <- df1.matches$id1
x2 <- df1.matches$id2





### 02) Bradley-Terry model
K <- 3L
sds <- rep(0.3, K) #### can change
M0 <- (sds %*% t(sds)) * (diag(0.3, K, K) + 0.7)

set.seed(1000L)
t0 <- system.time({
    obj0 <- mh_bt(
        y1, y2, x1, x2, m, M0,
        N = 2e+4L,
        thin = 1e+1L,
        burnin = 5e+4L,
        print_freq = 1e+3L
    )
})
obj0$time <- t0

gamma0 <- obj0$gamma_par
df0.strength <- tibble(
    id = seq(m) - 1L,
    mean = gamma0 %>% colMeans,
    sd = gamma0 %>% apply(2, sd),
    lower = gamma0 %>% apply(2, quantile, 0.025),
    upper = gamma0 %>% apply(2, quantile, 0.975)
) %>%
    left_join(df1.players, ., "id") %>%
    mutate(player = fct_reorder(player, mean))

gg0.strength <- df0.strength %>%
    top_n(30, mean) %>%
    ggplot() +
    geom_segment(aes(lower, player, xend = upper, yend = player)) +
    geom_point(aes(mean, player)) +
    theme_bw(15L) +
    labs(x = "(log-)strength")






### 03) whole spectrum of beta = 1.0 / eta
set.seed(2000L)
t1 <- system.time({
    obj1 <- mh_model(
        y1, y2, x1, x2, m, M0,
        N = 2e+4L,
        thin = 2e+3L,
        burnin = 5e+6L,
        print_freq = 2e+3L
    )
})
obj1$time <- t1

gamma1 <- obj1$gamma_par
df1.strength <- tibble(
    id = seq(m) - 1L,
    mean = gamma1 %>% colMeans,
    sd = gamma1 %>% apply(2, sd),
    lower = gamma1 %>% apply(2, quantile, 0.025),
    upper = gamma1 %>% apply(2, quantile, 0.975)
) %>%
    left_join(df1.players, ., "id") %>%
    mutate(player = fct_reorder(player, mean))

gg1.strength <- df1.strength %>%
    top_n(30, mean) %>%
    ggplot() +
    geom_segment(aes(lower, player, xend = upper, yend = player)) +
    geom_point(aes(mean, player)) +
    theme_bw(15L) +
    labs(x = "(log-)strength")





### 04) save & exit
## Bradley-Terry
obj0 %>% write_rds("obj0.rds")
df0.strength %>% write_csv("strengths_B_T.csv")
print(gg0.strength)
ggsave("Bradley_Terry_30.jpg")

## whole spectrum
obj1 %>% write_rds("obj1.rds")
df1.strength %>% write_csv("strengths_all.csv")
print(gg1.strength)
ggsave("Flexible_beta_30.jpg")

##
print(Sys.time())
