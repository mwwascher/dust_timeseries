library(patchwork)
library(ggplot2)


setwd("C:/Users/matth/OneDrive/Documents/JuliaWD/Library")

hist_beta <- data.frame(x = beta.post)

beta.post1 = read.csv("beta_samp1.csv")
beta.post2 = read.csv("beta_samp2.csv")
beta.post3 = read.csv("beta_samp3.csv")
beta.post4 = read.csv("beta_samp4.csv")
beta.post1 = as.numeric(beta.post1[[1]])
beta.post2 = as.numeric(beta.post2[[1]])
beta.post3 = as.numeric(beta.post3[[1]])
beta.post4 = as.numeric(beta.post4[[1]])

beta.post = c(beta.post1[5001:15000],beta.post2[5001:15000],beta.post3[5001:15000],beta.post4[5001:15000])

n.post.1 = read.csv("n_samp_1.csv")
n.post.2 = read.csv("n_samp_2.csv")
n.post.3 = read.csv("n_samp_3.csv")
n.post.4 = read.csv("n_samp_4.csv")

hist_df1 <- data.frame(x = n.post.1[5001:15000,1])
hist_df2 <- data.frame(x = n.post.2[5001:15000,1])
hist_df3 <- data.frame(x = n.post.3[5001:15000,1])
hist_df4 <- data.frame(x = n.post.4[5001:15000,1])

p0 <- ggplot(hist_beta, aes(x = x)) +
  geom_histogram(
    bins = 16,
    fill = "darkorchid",
    color = "black",
    alpha = 0.5
  ) +
  labs(
    x = "beta",
    y = "Frequency"
  ) +
  theme_classic() +
  theme(
    panel.grid = element_blank(),
    axis.title.x = element_text(
      size = 10,
      color = "black",
      margin = margin(t = 0.1, unit = "cm")
    ),
    axis.title.y = element_text(
      size = 10,
      color = "black",
      margin = margin(r = 0.1, unit = "cm")
    ),
    axis.text.x = element_text(size = 9, color = "black"),
    axis.text.y = element_text(size = 9, color = "black"),
    strip.background = element_rect(fill = "gray80")
  )

p0 <- p0 +
  scale_y_continuous(expand = c(0, 0))

p0


p1 <- ggplot(hist_df1, aes(x = x)) +
  geom_histogram(
    bins = 16,
    fill = "dodgerblue2",
    color = "black",
    alpha = 0.5
  ) +
  labs(
    x = "infected individuals",
    y = "Frequency"
  ) +
  theme_classic() +
  theme(
    panel.grid = element_blank(),
    axis.title.x = element_text(
      size = 10,
      color = "black",
      margin = margin(t = 0.1, unit = "cm")
    ),
    axis.title.y = element_text(
      size = 10,
      color = "black",
      margin = margin(r = 0.1, unit = "cm")
    ),
    axis.text.x = element_text(size = 9, color = "black"),
    axis.text.y = element_text(size = 9, color = "black"),
    strip.background = element_rect(fill = "gray80")
  )

p1 <- p1 +
  scale_y_continuous(expand = c(0, 0))


p1 <- p1 +
  geom_vline(
    xintercept = 6,
    color = "red",
    linewidth = 0.8,
    linetype = "dashed"
  )


p2 <- ggplot(hist_df2, aes(x = x)) +
  geom_histogram(
    bins = 16,
    fill = "dodgerblue2",
    color = "black",
    alpha = 0.5
  ) +
  labs(
    x = "infected individuals",
    y = "Frequency"
  ) +
  theme_classic() +
  theme(
    panel.grid = element_blank(),
    axis.title.x = element_text(
      size = 10,
      color = "black",
      margin = margin(t = 0.1, unit = "cm")
    ),
    axis.title.y = element_text(
      size = 10,
      color = "black",
      margin = margin(r = 0.1, unit = "cm")
    ),
    axis.text.x = element_text(size = 9, color = "black"),
    axis.text.y = element_text(size = 9, color = "black"),
    strip.background = element_rect(fill = "gray80")
  )

p2 <- p2 +
  scale_y_continuous(expand = c(0, 0))

p2 <- p2 +
  geom_vline(
    xintercept = 12,
    color = "red",
    linewidth = 0.8,
    linetype = "dashed"
  )

p3 <- ggplot(hist_df3, aes(x = x)) +
  geom_histogram(
    bins = 16,
    fill = "dodgerblue2",
    color = "black",
    alpha = 0.5
  ) +
  labs(
    x = "infected individuals",
    y = "Frequency"
  ) +
  theme_classic() +
  theme(
    panel.grid = element_blank(),
    axis.title.x = element_text(
      size = 10,
      color = "black",
      margin = margin(t = 0.1, unit = "cm")
    ),
    axis.title.y = element_text(
      size = 10,
      color = "black",
      margin = margin(r = 0.1, unit = "cm")
    ),
    axis.text.x = element_text(size = 9, color = "black"),
    axis.text.y = element_text(size = 9, color = "black"),
    strip.background = element_rect(fill = "gray80")
  )

p3 <- p3 +
  scale_y_continuous(expand = c(0, 0))

p3 <- p3 +
  geom_vline(
    xintercept = 20,
    color = "red",
    linewidth = 0.8,
    linetype = "dashed"
  )

p4 <- ggplot(hist_df4, aes(x = x)) +
  geom_histogram(
    bins = 12,
    fill = "dodgerblue2",
    color = "black",
    alpha = 0.5
  ) +
  labs(
    x = "infected individuals",
    y = "Frequency"
  ) +
  theme_classic() +
  theme(
    panel.grid = element_blank(),
    axis.title.x = element_text(
      size = 10,
      color = "black",
      margin = margin(t = 0.1, unit = "cm")
    ),
    axis.title.y = element_text(
      size = 10,
      color = "black",
      margin = margin(r = 0.1, unit = "cm")
    ),
    axis.text.x = element_text(size = 9, color = "black"),
    axis.text.y = element_text(size = 9, color = "black"),
    strip.background = element_rect(fill = "gray80")
  )

p4 <- p4 +
  scale_y_continuous(expand = c(0, 0))

p4 <- p4 +
  geom_vline(
    xintercept = 19,
    color = "red",
    linewidth = 0.8,
    linetype = "dashed"
  )

(p1 | p2) /
  (p3 | p4) +
  plot_annotation(
    title = "Library",
    theme = theme(
      plot.title = element_text(
        hjust = 0.5
      )
    )
  )
