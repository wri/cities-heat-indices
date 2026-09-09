z_norm <- function(x) {
  mu <- global(x, "mean", na.rm = TRUE)[[1]]
  sd <- global(x, "sd",   na.rm = TRUE)[[1]]
  (x - mu) / sd
}


cat_z_sd <- function(z) {
  
  breaks <- c(-Inf, -2, -1, 1, 2, Inf)
  
  classify(
    z,
    rcl = cbind(
      breaks[-length(breaks)],
      breaks[-1],
      1:5
    )
  )
}



lst_z  <- z_norm(lst)
pop_z  <- z_norm(world_pop)
alb_z  <- z_norm(alb)
tree_z <- z_norm(tree)
fr_z   <- z_norm(fr)
# Invert protective infrastructure layers so higher = more vulnerable
alb_z_v  <- -alb_z
tree_z_v <- -tree_z
fr_z_v   <- -fr_z
# Categorize each input layer into 6 classes
lst_z_cat  <- cat_z_sd(lst_z)
pop_z_cat  <- cat_z_sd(pop_z)
alb_z_cat  <- cat_z_sd(alb_z_v)
tree_z_cat <- cat_z_sd(tree_z_v)
fr_z_cat   <- cat_z_sd(fr_z_v)
# Infrastructure component (continuous + 6 classes)
infra_z_cat <- (alb_z_cat + tree_z_cat + fr_z_cat) / 3

# Final index (continuous + 6 classes)
HVI_z <- (lst_z_cat + pop_z_cat + infra_z_cat) / 3
plot(HVI_z)

p_lst <- ggplot() +
  tidyterra::geom_spatraster(data = lst_z_cat) +
  scale_fill_viridis_c(
    option = "magma",
    name = "z score"
  ) +
  ggtitle("LST")
p_infra <- ggplot() +
  tidyterra::geom_spatraster(data = infra_z_cat) +
  scale_fill_viridis_c(
    option = "magma",
    name = "z score"
  ) +
  ggtitle("Infrastructure")
p_pop <- ggplot() +
  tidyterra::geom_spatraster(data = pop_z_cat) +
  scale_fill_viridis_c(
    option = "magma",
    name = "z score"
  ) +
  ggtitle("Population")
p_hvi <- ggplot() +
  tidyterra::geom_spatraster(data = HVI_z) +
  scale_fill_viridis_c(
    option = "magma",
    name = "z score"
  ) +
  ggtitle("HVI")

p_lst + p_infra + p_pop + p_hvi


p_lst_s <- ggplot() +
  tidyterra::geom_spatraster(data = lst) +
  scale_fill_viridis_c(
    option = "magma"
  ) +
  ggtitle("LST")
p_alb <- ggplot() +
  tidyterra::geom_spatraster(data = alb) +
  scale_fill_viridis_c(
    option = "magma",
    direction = -1) +
  ggtitle("Albedo")
p_fr <- ggplot() +
  tidyterra::geom_spatraster(data = fr) +
  scale_fill_viridis_c(
    option = "magma",
    direction = -1) +
  ggtitle("Fr")
p_tree <- ggplot() +
  tidyterra::geom_spatraster(data = tree) +
  scale_fill_viridis_c(
    option = "magma",
    direction = -1) +
  ggtitle("Tree")
p_pop_s <- ggplot() +
  tidyterra::geom_spatraster(data = world_pop) +
  scale_fill_viridis_c(
    option = "magma") +
  ggtitle("Population")


p_lst_s + p_alb + p_fr + p_tree + p_pop_s 

