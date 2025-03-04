# Use ggpmisc package for regression fitting, including regression equation, R², and p-value
library(ggpmisc)
library(ggplot2)
library(ggrepel)

# Load data
arg <- read.csv('hidden_file_1.csv', header = TRUE, row.names = 1, sep = ',')
mge <- read.csv('hidden_file_2.csv', header = TRUE, row.names = 1, sep = ',')
genus <- read.table('hidden_file_3.txt', header = TRUE, row.names = 1, sep = '\t')

group <- read.table('hidden_file_4.txt', sep = '\t', header = TRUE)

### Data processing, adjust as needed ###
# Select the first 35 columns of the 'arg' and 'genus' data frames
# arg <- arg[, 1:35]
# genus <- genus[, 1:35]
# Calculate the sum of each column for arg and mge
arg_sum <- colSums(arg)
mge_sum <- colSums(mge)
# Calculate the count of non-zero values for each column in arg (detection count)
arg_counts <- colSums(arg != 0)
# Calculate the count of non-zero values for each column in mge (detection count)
mge_counts <- colSums(mge != 0)
# Calculate the count of non-zero values for each column in genus (detection count)
genus_counts <- colSums(t(genus) != 0)

# Create a new data frame
new_data <- data.frame(colnames = names(arg_sum), arg_sum = arg_sum, mge_sum = mge_sum,
                       genus_counts = genus_counts, mge_counts = mge_counts, arg_counts = arg_counts)
# Merge the 'group' data frame with the new data frame using the sample identifier
# new_data <- merge(new_data, group[1:35, c(1, 3)], by.x = "colnames", by.y = "V3")
new_data <- merge(new_data, group, by.x = "colnames", by.y = "Sample")

# Apply logarithmic transformation to gene abundance data
new_data$log_arg_sum <- log(new_data$arg_sum)
new_data$log_mge_sum <- log(new_data$mge_sum)
# Define color mapping rules
# color_rules <- c("chicken_gut" = "#E7DAD2", "sludge" = "#4daf4a", "water" = "#82B0D2", "fish_gut" = "#BEB8DC")
color_rules <- c("ck" = "#A6CEE3", "F" = "#B2DF8A", "pollution" = "#E31A1C")

# Plot scatter plot
# p1: no grouping, all samples fitted together
p1 <- ggplot(data = new_data, aes(x = log_arg_sum, y = log_mge_sum)) +
  geom_point(aes(color = LandM)) +
  scale_color_manual(values = color_rules) +
  theme_bw() +
  # geom_text adds labels; geom_text_repel prevents label overlapping,
  # nudge_y adjusts label position, size sets label size, max.overlaps limits maximum overlaps
  geom_text_repel(aes(label = colnames), nudge_y = 0.05, size = 4, max.overlaps = 15)
p1

### Identify obvious outliers (e.g., A2, A3) via fitted curve and remove them ###
# Adjust this step as needed based on your data
# new_data <- subset(new_data, !colnames %in% c("f3", "f4", "f15"))
new_data <- subset(new_data, !colnames %in% c("A2", "A3"))

# Plot again. p1_1: arg abundance vs. mge abundance
p1_1 <- ggplot(data = new_data, aes(x = log_arg_sum, y = log_mge_sum)) +
  geom_point(aes(color = LandM), alpha = 0.5, size = 3) +
  scale_color_manual(values = color_rules) +
  theme_bw()
  # To add labels without overlapping:
  # geom_text_repel(aes(label = colnames), nudge_y = 0.05, size = 4)
p1_1

# Simple example to distinguish groups with different colors. p1: arg abundance vs. mge abundance
p1 <- ggplot(data = new_data, aes(x = log_arg_sum, y = log_mge_sum, 
                                  group = LandM, color = LandM)) +
  geom_point(alpha = 0.5, size = 3) +
  scale_color_manual(values = color_rules) +
  theme_bw()
p1
# Note: p1 and p1_1 appear identical at this point, but p1_1 will later be fitted separately for each group.
#
#
# Examine the relationship between total arg abundance and genus counts. p1_2: arg abundance vs. genus counts
p1_2 <- ggplot(data = new_data, aes(x = log_arg_sum, y = genus_counts)) +
  geom_point(aes(group = LandM, color = LandM), alpha = 0.5, size = 3) +
  scale_color_manual(values = color_rules) +
  theme_bw()
p1_2

# Examine the relationship between total arg abundance and mge counts. p1_3: arg abundance vs. mge count
p1_3 <- ggplot(data = new_data, aes(x = log_arg_sum, y = mge_counts)) +
  geom_point(aes(group = LandM, color = LandM), alpha = 0.5, size = 3) +
  scale_color_manual(values = color_rules) +
  theme_bw()
p1_3


###
#
#
# Examine relationships with arg count. p2_1: arg count vs. mge abundance
p2_1 <- ggplot(data = new_data, aes(x = arg_counts, y = log_mge_sum)) +
  geom_point(aes(group = LandM, color = LandM), alpha = 0.5, size = 3) +
  scale_color_manual(values = color_rules) +
  theme_bw()
p2_1

# p2_2: arg count vs. genus counts
p2_2 <- ggplot(data = new_data, aes(x = arg_counts, y = genus_counts)) +
  geom_point(aes(group = LandM, color = LandM), alpha = 0.5, size = 3) +
  scale_color_manual(values = color_rules) +
  theme_bw()
p2_2

# p2_3: arg count vs. mge count
p2_3 <- ggplot(data = new_data, aes(x = arg_counts, y = mge_counts)) +
  geom_point(aes(group = LandM, color = LandM), alpha = 0.5, size = 3) +
  scale_color_manual(values = color_rules) +
  theme_bw()
p2_3

#
#
#
# Perform linear fitting: arg abundance vs. mge abundance
# Compute the coefficients using a linear model (lm)
fit <- lm(log_mge_sum ~ log_arg_sum, data = new_data)
# Compute the confidence interval; 'confidence' specifies the confidence interval level (here 99%)
conf_int <- predict(fit, interval = "confidence", level = 0.99)
# Plot
p1_1 + 
  # Parameters for the fitted curve
  geom_line(aes(y = fit$fitted.values), alpha = 0.65, color = "red", linewidth = 1.4) +
  # Parameters for the confidence interval: alpha sets transparency, fill sets the interval color
  geom_ribbon(aes(ymin = conf_int[, "lwr"], ymax = conf_int[, "upr"]), alpha = 0.4, fill = "lightblue") +
  # Display the fitted equation on the plot
  stat_poly_eq(aes(label = paste(eq.label, ..rr.label.., after_stat(p.value.label), sep = '~`,`~')),
               formula = y ~ x, parse = TRUE, label.x.npc = 'left', label.y.npc = 'top', size = 2.7)

# Perform linear fitting: arg abundance vs. mge count
fit <- lm(mge_counts ~ log_arg_sum, data = new_data)
conf_int <- predict(fit, interval = "confidence", level = 0.99)
# Plot
p1_3 + 
  geom_line(aes(y = fit$fitted.values), alpha = 0.65, color = "red", linewidth = 1.4) +
  geom_ribbon(aes(ymin = conf_int[, "lwr"], ymax = conf_int[, "upr"]), alpha = 0.4, fill = "lightblue") +
  stat_poly_eq(aes(label = paste(eq.label, ..rr.label.., after_stat(p.value.label), sep = '~`,`~')),
               formula = y ~ x, parse = TRUE, label.x.npc = 'left', label.y.npc = 'top', size = 2.7)

# Perform linear fitting: arg abundance vs. genus counts
fit <- lm(genus_counts ~ log_arg_sum, data = new_data)
conf_int <- predict(fit, interval = "confidence", level = 0.99)
# Plot
p1_2 + 
  geom_line(aes(y = fit$fitted.values), alpha = 0.65, color = "red", linewidth = 1.4) +
  geom_ribbon(aes(ymin = conf_int[, "lwr"], ymax = conf_int[, "upr"]), alpha = 0.4, fill = "lightblue") +
  stat_poly_eq(aes(label = paste(eq.label, ..rr.label.., after_stat(p.value.label), sep = '~`,`~')),
               formula = y ~ x, parse = TRUE, label.x.npc = 'left', label.y.npc = 'top', size = 2.7)

#
# Analysis related to arg count
#
# Perform linear fitting: arg count vs. mge abundance
fit <- lm(log_mge_sum ~ arg_counts, data = new_data)
conf_int <- predict(fit, interval = "confidence", level = 0.99)
# Plot
p2_1 + 
  geom_line(aes(y = fit$fitted.values), alpha = 0.65, color = "red", linewidth = 1.4) +
  geom_ribbon(aes(ymin = conf_int[, "lwr"], ymax = conf_int[, "upr"]), alpha = 0.4, fill = "lightblue") +
  stat_poly_eq(aes(label = paste(eq.label, ..rr.label.., after_stat(p.value.label), sep = '~`,`~')),
               formula = y ~ x, parse = TRUE, label.x.npc = 'left', label.y.npc = 'top', size = 2.7)

# Perform linear fitting: arg count vs. mge count
fit <- lm(mge_counts ~ arg_counts, data = new_data)
conf_int <- predict(fit, interval = "confidence", level = 0.99)
# Plot
p2_3 + 
  geom_line(aes(y = fit$fitted.values), alpha = 0.65, color = "red", linewidth = 1.4) +
  geom_ribbon(aes(ymin = conf_int[, "lwr"], ymax = conf_int[, "upr"]), alpha = 0.4, fill = "lightblue") +
  stat_poly_eq(aes(label = paste(eq.label, ..rr.label.., after_stat(p.value.label), sep = '~`,`~')),
               formula = y ~ x, parse = TRUE, label.x.npc = 'left', label.y.npc = 'top', size = 2.7)

#
# Perform linear fitting: arg count vs. genus counts
fit <- lm(genus_counts ~ arg_counts, data = new_data)
conf_int <- predict(fit, interval = "confidence", level = 0.99)
# Plot
p2_2 + 
  geom_line(aes(y = fit$fitted.values), alpha = 0.65, color = "red", linewidth = 1.4) +
  geom_ribbon(aes(ymin = conf_int[, "lwr"], ymax = conf_int[, "upr"]), alpha = 0.4, fill = "lightblue") +
  stat_poly_eq(aes(label = paste(eq.label, ..rr.label.., after_stat(p.value.label), sep = '~`,`~')),
               formula = y ~ x, parse = TRUE, label.x.npc = 'left', label.y.npc = 'top', size = 2.7)

  
# Grouped plotting.
p3 <- ggplot(data = new_data, aes(x = log_arg_sum, y = log_mge_sum, 
                                  group = LandM, color = LandM)) +
  geom_point(alpha = 0.5, size = 3) +
  scale_color_manual(values = color_rules) +
  theme_bw()

p3 + stat_poly_eq(aes(label = paste(..eq.label.., ..rr.label.., stat(p.value.label), sep = '~`,`~')),
               formula = y ~ x, parse = TRUE, label.x.npc = 'left', label.y.npc = 'top', size = 2.7) + 
  geom_smooth(method = 'lm', formula = y ~ x, se = TRUE, show.legend = FALSE,
              # Fitted curve parameters: line size, confidence interval transparency, fill color, confidence level, and linetype to hide fitted line
              alpha = 0.4, aes(color = LandM), size = 1, fill = "#f7f3e8", level = 0.99, linetype = "blank") +
  # If se = FALSE, the confidence interval is not plotted
  geom_smooth(method = 'lm', formula = y ~ x, se = FALSE, show.legend = FALSE,
              aes(color = LandM), size = 1) +
  stat_poly_eq(aes(label = paste(..eq.label.., ..rr.label.., stat(p.value.label), sep = '~`,`~')),
               formula = y ~ x, parse = TRUE, label.x.npc = 'left', label.y.npc = 'top', size = 2.7)

# Save output as needed
