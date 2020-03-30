library(readr)
library(dplyr)
library(ggplot2)

busco_output <- list.files("output/busco",
                           recursive = TRUE,
                           pattern = "short_summary.txt",
                           full.names = TRUE)

busco_results_list <- lapply(busco_output,
                             readr::read_tsv,
                             skip=2)
##this woill need to change - file no longer named with filter, but a folder a couple steps above
names(busco_results_list) <- gsub("output/busco/(.+)/run_endopterygota_odb10/full_table.tsv", "\\1", busco_output)

busco_all <- dplyr::bind_rows(busco_results_list, .id = "filename")

plot_data <- busco_all %>%
  group_by(filename, Status) %>%
  summarise(n = n()) %>%
  group_by(filename) %>%
  mutate(percentage = n/sum(n)*100)

#order variables
status_order <- c("Complete", "Fragmented", "Duplicated", "Missing")
plot_data$Status <- factor(plot_data$Status, levels = rev(status_order))

ggplot(plot_data, aes(x = filename, y = percentage, fill = Status)) +
  xlab(NULL) + ylab(NULL) + 
  theme(axis.text.x = element_text(angle = 45, hjust = 1)) +
  geom_col()
