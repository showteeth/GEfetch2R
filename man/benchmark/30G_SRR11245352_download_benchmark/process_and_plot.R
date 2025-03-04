# process the output

# library
library(tidyverse)

# read the download time
sra_down = read.table('./download/SRA_time_down.txt', header = T)
ena_down = read.table('./download/ENA_time_down.txt', header = T)
# remove the failed one
ena_down = ena_down %>% dplyr::filter(!(type == "sra" & batch =="9"))
all_down = rbind(sra_down, ena_down) %>% as.data.frame()
all_down$tool = factor(all_down$tool, levels = c("prefetch", "ascp", "ascp_fastq"))
all_down_plot = ggplot(all_down, aes(x = tool, y = time, color = tool)) +
  geom_boxplot(show.legend = FALSE) +
  geom_jitter(show.legend = FALSE) +
  theme_bw() +
  scale_x_discrete(labels=c("prefetch" = "prefetch (sra)", "ascp" = "ascp (sra)",
                            "ascp_fastq" = "ascp (fastq)")) +
  scale_y_continuous(breaks = seq(0, 6600, by = 500))
ggsave(filename = "all_down_plot.pdf", plot = all_down_plot, height = 5, width = 5)
write.table(x = all_down, file = 'SRA_ENA_time_down.txt', sep = "\t", row.names = F, quote = F)

# read the split time
all_split_file = list.files(path = "./split", pattern = "txt$", recursive = T, full.names = T)
all_split_list = lapply(all_split_file, function(x){
  read.table(file = x, header = T)
})
all_split_df = do.call(rbind, all_split_list)
# fastq
all_split_fq_df = all_split_df %>%
  dplyr::filter(type == "fq")
all_split_fq_df$tr = factor(all_split_fq_df$tr, levels = c(1,2,4,8,12))
all_split_fq_df$tool = factor(all_split_fq_df$tool, levels = c("fastq-dump", "fasterq-dump", "parallel-fastq-dump"))
# calculate mean tiem
all_split_fq_mean <- all_split_fq_df %>%
  dplyr::group_by(tool, tr) %>%
  dplyr::summarize(average = mean(time)) %>%
  as.data.frame()
all_split_fq_df_plot =
  ggplot() +
  geom_boxplot(data = all_split_fq_df, mapping = aes(x = tr, y = time, fill = tool), width=0.75) +
  geom_point(data = all_split_fq_mean,
             mapping = aes(x = tr, y = average, fill = tool), position = position_dodge(width=0.75),
             pch = 21, color = "black",size = 1.5) +
  geom_line(data = all_split_fq_mean,
            mapping = aes(x = tr, y = average, group = tool, color = tool), position = position_dodge(width=0.75)) +
  theme_bw() +
  scale_y_continuous(breaks = seq(0, 2200, by = 200)) +
  scale_fill_manual(values = c("fastq-dump" = "#F8766D", "fasterq-dump" = "#7CAE00", "parallel-fastq-dump" = "#C77CFF")) +
  scale_color_manual(values = c("fastq-dump" = "#F8766D", "fasterq-dump" = "#7CAE00", "parallel-fastq-dump" = "#C77CFF"))
ggsave(filename = "all_split_fq_df_plot.pdf", plot = all_split_fq_df_plot, height = 6, width = 7.5)

# fastq.gz
all_split_gz_df = all_split_df %>%
  dplyr::filter(type != "fq")
all_split_gz_df$tr = factor(all_split_gz_df$tr, levels = c(1,2,4,8,12))
all_split_gz_df = all_split_gz_df %>%
  dplyr::mutate(group = ifelse(tool == "fasterq-dump", paste(tool, type, sep = "_"), tool))
all_split_gz_df$group = factor(all_split_gz_df$group, levels = c("fastq-dump", "fasterq-dump_gzip", "fasterq-dump_pigz", "parallel-fastq-dump"))
all_split_gz_mean <- all_split_gz_df %>%
  dplyr::group_by(group, tr) %>%
  dplyr::summarize(average = mean(time)) %>%
  as.data.frame()
all_split_gz_df_plot =
  ggplot() +
  geom_boxplot(data = all_split_gz_df, mapping = aes(x = tr, y = time, fill = group), width=0.75) +
  geom_point(data = all_split_gz_mean,
             mapping = aes(x = tr, y = average, fill = group), position = position_dodge(width=0.75),
             pch = 21, color = "black",size = 1.5) +
  geom_line(data = all_split_gz_mean,
            mapping = aes(x = tr, y = average, group = group, color = group), position = position_dodge(width=0.75)) +
  theme_bw() +
  scale_y_continuous(breaks = seq(0, 14000, by = 1000))
ggsave(filename = "all_split_gz_df_plot.pdf", plot = all_split_gz_df_plot, height = 6, width = 8)
write.table(x = all_split_df, file = 'SRA_time_split.txt', sep = "\t", row.names = F, quote = F)
