library(tidyverse)


p1 <- read_csv('Path/to/p1_clone_counts.csv')
p2 <- read_csv('Path/to/p2_clone_counts.csv')
p3 <- read_csv('Path/to/p3_clone_counts.csv')
p4 <- read_csv('Path/to/p4_clone_counts.csv')
p5 <- read_csv('Path/to/p5_clone_counts.csv')
p6 <- read_csv('Path/to/p6_clone_counts.csv')
p7 <- read_csv('Path/to/p7_clone_counts.csv')
t1 <- read_csv('Path/to/t1_clone_counts.csv')
t2 <- read_csv('Path/to/t2_clone_counts.csv')
t3 <- read_csv('Path/to/t3_clone_counts.csv')
t4 <- read_csv('Path/to/t4_clone_counts.csv')
t5 <- read_csv('Path/to/t5_clone_counts.csv')
t6 <- read_csv('Path/to/t6_clone_counts.csv')
t7 <- read_csv('Path/to/t7_clone_counts.csv')

make_bubble <- function(data) {
  plot_table <- data %>%
    mutate(freq = count/sum(count)) %>%
    rownames_to_column(var = 'id') %>%
    mutate(id = as.integer(as.numeric(id)))
  p <- circleProgressiveLayout(plot_table$count)
  d <- circleLayoutVertices(p) %>%
    left_join(plot_table, by = 'id')
  
  return(d)
}

p1_bubble <- make_bubble(p1)
p2_bubble <- make_bubble(p2)
p3_bubble <- make_bubble(p3)
p4_bubble <- make_bubble(p4)
p5_bubble <- make_bubble(p5)
p6_bubble <- make_bubble(p6)
p7_bubble <- make_bubble(p7)
t1_bubble <- make_bubble(t1)
t2_bubble <- make_bubble(t2)
t3_bubble <- make_bubble(t3)
t4_bubble <- make_bubble(t4)
t5_bubble <- make_bubble(t5)
t6_bubble <- make_bubble(t6)
t7_bubble <- make_bubble(t7)

library(RColorBrewer)
coloring <- setNames(c( "#499b80",'#2a6958', "#0c3431"), c('None', 'Partly', 'Pair'))

plot_p1  <- ggplot(p1_bubble, aes(x, y)) + geom_polygon(aes(group = id, fill = common), show.legend = F) + theme_void()  + coord_cartesian(ylim = c(-40,40), xlim = c(-40,40)) + scale_fill_manual(values = coloring)
plot_p2  <- ggplot(p2_bubble, aes(x, y)) + geom_polygon(aes(group = id, fill = common), show.legend = F) + theme_void()  + coord_cartesian(ylim = c(-40,40), xlim = c(-40,40)) + scale_fill_manual(values = coloring)
plot_p3  <- ggplot(p3_bubble, aes(x, y)) + geom_polygon(aes(group = id, fill = common), show.legend = F) + theme_void()  + coord_cartesian(ylim = c(-40,40), xlim = c(-40,40)) + scale_fill_manual(values = coloring)
plot_p4  <- ggplot(p4_bubble, aes(x, y)) + geom_polygon(aes(group = id, fill = common), show.legend = F) + theme_void()  + coord_cartesian(ylim = c(-40,40), xlim = c(-40,40)) + scale_fill_manual(values = coloring)
plot_p5  <- ggplot(p5_bubble, aes(x, y)) + geom_polygon(aes(group = id, fill = common), show.legend = F) + theme_void()  + coord_cartesian(ylim = c(-40,40), xlim = c(-40,40)) + scale_fill_manual(values = coloring)
plot_p6  <- ggplot(p6_bubble, aes(x, y)) + geom_polygon(aes(group = id, fill = common), show.legend = F) + theme_void()  + coord_cartesian(ylim = c(-40,40), xlim = c(-40,40)) + scale_fill_manual(values = coloring)
plot_p7  <- ggplot(p7_bubble, aes(x, y)) + geom_polygon(aes(group = id, fill = common), show.legend = F) + theme_void()  + coord_cartesian(ylim = c(-40,40), xlim = c(-40,40)) + scale_fill_manual(values = coloring)
plot_t1  <- ggplot(t1_bubble, aes(x, y)) + geom_polygon(aes(group = id, fill = common), show.legend = F) + theme_void()  + coord_cartesian(ylim = c(-40,40), xlim = c(-40,40)) + scale_fill_manual(values = coloring)
plot_t2  <- ggplot(t2_bubble, aes(x, y)) + geom_polygon(aes(group = id, fill = common), show.legend = F) + theme_void()  + coord_cartesian(ylim = c(-40,40), xlim = c(-40,40)) + scale_fill_manual(values = coloring)
plot_t3  <- ggplot(t3_bubble, aes(x, y)) + geom_polygon(aes(group = id, fill = common), show.legend = F) + theme_void()  + coord_cartesian(ylim = c(-40,40), xlim = c(-40,40)) + scale_fill_manual(values = coloring)
plot_t4  <- ggplot(t4_bubble, aes(x, y)) + geom_polygon(aes(group = id, fill = common), show.legend = F) + theme_void()  + coord_cartesian(ylim = c(-40,40), xlim = c(-40,40)) + scale_fill_manual(values = coloring)
plot_t5  <- ggplot(t5_bubble, aes(x, y)) + geom_polygon(aes(group = id, fill = common), show.legend = F) + theme_void()  + coord_cartesian(ylim = c(-40,40), xlim = c(-40,40)) + scale_fill_manual(values = coloring)
plot_t6  <- ggplot(t6_bubble, aes(x, y)) + geom_polygon(aes(group = id, fill = common), show.legend = F) + theme_void()  + coord_cartesian(ylim = c(-40,40), xlim = c(-40,40)) + scale_fill_manual(values = coloring)
plot_t7  <- ggplot(t7_bubble, aes(x, y)) + geom_polygon(aes(group = id, fill = common), show.legend = F) + theme_void()  + coord_cartesian(ylim = c(-40,40), xlim = c(-40,40)) + scale_fill_manual(values = coloring)


plot_grid(plot_t1, plot_t2, plot_t3, plot_t4, plot_t5, plot_t6, plot_t7, plot_p1, plot_p2, plot_p3, plot_p4, plot_p5, plot_p6, plot_p7, nrow = 2)
