# Load libraries
library(ggplot2)

# Example data: 5 SNPs with effect sizes and standard errors
dat <- data.frame(
 SNP = paste0("rs", 1:6),
 beta_exposure = c(0.01, 0.02, 0.015, 0.025, 0.012, 0.018),
 se_exposure = c(0.005, 0.006, 0.005, 0.007, 0.005, 0.006),
 beta_outcome = c(0.02, 0.03, 0.022, 0.035, 0.018, 0.028),
 se_outcome = c(0.01, 0.012, 0.011, 0.013, 0.01, 0.012)
)


# IVW slope (forced through origin)
weights <- 1 / (dat$se_outcome^2)
ivw_beta <- sum(weights * dat$beta_exposure * dat$beta_outcome) / sum(weights * dat$beta_exposure^2)

# Determine maximum values for x and y (with some padding)
x_max <- max(dat$beta_exposure + dat$se_exposure) * 1.1
y_max <- max(dat$beta_outcome + dat$se_outcome) * 1.1

# Create plot
ggplot(dat, aes(x = beta_exposure, y = beta_outcome)) +
 geom_point(size = 3) +
 geom_errorbar(aes(ymin = beta_outcome - se_outcome, ymax = beta_outcome + se_outcome), width = 0) +
 geom_errorbarh(aes(xmin = beta_exposure - se_exposure, xmax = beta_exposure + se_exposure), height = 0) +
 geom_abline(intercept = 0, slope = ivw_beta, color = "#89CFF0", size = 1.2) +  # pastel blue
 labs(
  x = "SNP effect on exposure",
  y = "SNP effect on outcome",
  title = ""
 ) +
 scale_x_continuous(limits = c(0, x_max), expand = c(0, 0)) +
 scale_y_continuous(limits = c(0, y_max), expand = c(0, 0)) +
 theme_minimal(base_size = 14) +
 theme(
  panel.grid = element_blank(),
  axis.line = element_line(color = "black"),
  axis.ticks = element_line(color = "black"),
  panel.border = element_blank()
 )
