library(here)
library(tidyverse)
library(latex2exp)
# library(cowplot)

source(here("expr", "utils.R"))


# fdp and power of the methods
draw_fdp_power_curve <- function(experiment, X_types, sample_size = 1,
                                 method_names, method_colors, method_shapes,
                                 error_bar = F){
    load(here("data", paste0(experiment, ".Rdata")))
    
    fdr_power <- lapply(X_types, function(X_type){
        results[[X_type]]$FDR_Power %>% mutate(design_mat = str_replace(X_type, "_", "-"))
    })
    fdr_power <- do.call(rbind, fdr_power)
    X_types <- str_replace(X_types, "_", "-")
    rm(results)
    
    fdr_power <- fdr_power %>% 
        filter(design_mat %in% X_types,
               method_name %in% method_names) %>%
        arrange(methods)

    methods_level <- unique(fdr_power$methods)
    # alphas <- unique(fdr_power$alpha)

    method_names <- parse_name(method_names)

    ref_prototype <- data.frame(fig_x = c(fig_x_var$value, fig_x_var$value),
                                threshold = c(unlist(alphas), rep(NA, length(fig_x_var$value))),
                                type = rep(c("FDR", "Power"), each = length(fig_x_var$value)))
    reference <- lapply(X_types, function(X_type){
        ref_prototype %>% mutate(design_mat = X_type)
    })
    reference <- do.call(rbind, reference)
    
    if(error_bar){
        add_error_bar <- quote(
            geom_errorbar(aes(x = fig_x, y = mean,
                              ymin = mean-2*std/sqrt(sample_size),
                              ymax = mean+2*std/sqrt(sample_size),
                              color = methods), width=0.05,
                          position = position_dodge(width=0.01))
        )
    } else{
        add_error_bar <- NULL
    }
    
    if(fig_x_var$value[1] <= fig_x_var$value[length(fig_x_var$value)]){
        set_x_axis <- scale_x_continuous
    } else{
        set_x_axis <- scale_x_reverse
    }

    plot <- ggplot(fdr_power) +
        geom_line(aes(x = fig_x, y = mean, color = methods)) +
        geom_point(aes(x = fig_x, y = mean, color = methods, shape = methods), size = 2) +
        eval(add_error_bar) +
        geom_line(data = reference, aes(x = fig_x, y = threshold),
                  linetype = "longdash", alpha = 0.6, na.rm = T) +
        facet_grid(vars(factor(design_mat, levels = X_types)),
                   vars(factor(type, levels = c("FDR", "Power"))), scales="free") +
        # facet_wrap(vars(factor(design_mat, levels = X_types), factor(type, levels = c("FDR", "Power"))),
        #            ncol = 2, scales="free_y") +
        set_x_axis(breaks = fig_x_var$value, labels = fig_x_var$value) +
        scale_color_manual(values = method_colors, labels = method_names, breaks = methods_level) +
        scale_shape_manual(values = method_shapes, labels = method_names, breaks = methods_level) +
        theme_bw() +
        theme(aspect.ratio = 1,
              panel.grid = element_blank(),
              strip.text = element_text(size = 15),
              axis.title = element_text(size = 13),
              axis.text = element_text(size = 10),
              legend.position = "right",
              legend.title=element_text(size=9),
              legend.text=element_text(size=9)) +
        labs(x = fig_x_var$name, y = "Estimated FDR/Power")

    ggsave(filename = here("latex", "figs", paste0("simu-", experiment, ".pdf")),
           plot, width = 7, height = 2*(length(X_types)+1))
    
    return(plot)
}





