barplot.missing <- function(obj){
        info <- t(
                info_table(obj, type = "missing", plot = F, plotlab = F) )

info <- as.data.frame(info[,-ncol(info)])

info$locus <- rownames(info)

info <- gather(info, key = pop, value = missing_data_freq, -locus)

mis <- info %>% 
        ggplot(aes(x=locus, y=missing_data_freq)) +
        geom_bar(stat = "identity") +
        geom_hline(yintercept = 0.10, colour = "red") +
        theme(axis.text.x = element_text(colour = NA),
              axis.ticks = element_line(colour = NA)) +
        ggtitle("Missing data by cohort after filtering") +
        facet_wrap(~pop, nrow=2 )

mis

}