## AP. Woodward, University of Georgia, 2025.
# Implementation of statistical analysis for a dose-response experiment for a single agent (vinblastine) using a 4-parameter log-logistic model.
# The statistical model is a multilevel nonlinear model with a mixture of hierarchical and non-hierarchical grouping features.
# The analysis is implemented in Stan (https://doi.org/10.18637/jss.v076.i01) with post-processing via packages 'ggplot2' (https://ggplot2.tidyverse.org) and 'ggdist' (https://doi.org/10.1109/TVCG.2023.3327195).

# Load the required packages.
library(rstan)
library(ggplot2)
library(rlist)
library(matrixStats)
library(viridis)
library(ggdist)

# Prepare the working data.
vinblastine_data <- single_dose_chemo_data_AW[single_dose_chemo_data_AW$DRUG == 'Vinblastine',]
vinblastine_data   <- rbind(vinblastine_data, single_dose_dead_controls[(single_dose_dead_controls$PLATE_NUM %in% vinblastine_data$PLATE_NUM),])
vinblastine_data$PLATE_IN <- as.numeric(as.factor(vinblastine_data$PLATE_NUM))
vinblastine_data$SUBJECT_IN <- as.numeric(factor(vinblastine_data$SUBJECT, levels = unique(vinblastine_data$SUBJECT)))
vinblastine_FIRST <- numeric(length = length(unique(vinblastine_data$PLATE_IN[!(vinblastine_data$DRUG == 'DMSO')])))
vinblastine_LAST <- numeric(length = length(unique(vinblastine_data$PLATE_IN[!(vinblastine_data$DRUG == 'DMSO')])))
for (i in 1:length(unique(vinblastine_data$PLATE_IN))){
  PLATE_iter <- which(vinblastine_data$PLATE_IN[!(vinblastine_data$DRUG == 'DMSO')] == (unique(vinblastine_data$PLATE_IN)[i]))
  vinblastine_FIRST[i] <- PLATE_iter[1]
  vinblastine_LAST[i] <- tail(PLATE_iter, n = 1)
}

# Initialize the input data for the Stan model, set the intial values, and conduct the parameter estimation.
vinblastine_list <- list(Nsub = length(unique(vinblastine_data$SUBJECT_IN)), Nplate = length(unique(vinblastine_data$PLATE_IN)), P = length(which(!(vinblastine_data$DRUG == 'DMSO'))), response_exp = vinblastine_data$Difference[!(vinblastine_data$DRUG == 'DMSO')],
                         plate_exp = vinblastine_data$PLATE_IN[!(vinblastine_data$DRUG == 'DMSO')], subject_exp = vinblastine_data$SUBJECT_IN[!(vinblastine_data$DRUG == 'DMSO')], dose = vinblastine_data$CONC_uM[!(vinblastine_data$DRUG == 'DMSO')],
                         plate_ind = vinblastine_data$PLATE_IN[!duplicated(vinblastine_data$PLATE_IN)], subj_ind = vinblastine_data$SUBJECT_IN[!duplicated(vinblastine_data$PLATE_IN)], D = length(which((vinblastine_data$DRUG == 'DMSO'))),
                         response_dead = vinblastine_data$Difference[(vinblastine_data$DRUG == 'DMSO')], subject_dead = vinblastine_data$SUBJECT_IN[(vinblastine_data$DRUG == 'DMSO')], plate_dead = vinblastine_data$PLATE_IN[(vinblastine_data$DRUG == 'DMSO')],
                         first_ind = vinblastine_FIRST, last_ind = vinblastine_LAST)
vinblastine_init <- function(){list(dead_resp_mu = runif(4,-2,-1), dead_resp_sd = runif(4,0.1,0.3))}
vinblastine_mod1 <- stan('doseresp_1d.stan', data = vinblastine_list, init = vinblastine_init, chains = 4, cores = 4, iter = 4000, warmup = 2000, control = list(adapt_delta = 0.98, max_treedepth = 12))
expose_stan_functions(vinblastine_mod1)

# Visualize the posterior probability distributions for the individual predictions at the plate level, with the observations.
vinblastine_labels <- c(`1` = 'Patient 6 Tissue', `2` = 'Patient 6 Urine', `3` = 'Patient 7 Tissue', `4` = 'Patient 7 Urine')
vinblastine_mod1_ind_samp <- extract(vinblastine_mod1, pars = c('Emin_VEC', 'Emax_VEC', 'EC50_VEC', 'H_VEC'))
vinblastine_mod1_ind_list <- vector(mode = 'list', length = length(unique(vinblastine_data$PLATE_IN)))
for (j in 1:length(vinblastine_mod1_ind_list)){
  iter_mat <- matrix(nrow = length(seq(-4,-1,0.1)), ncol = (dim(vinblastine_mod1_ind_samp$Emin_VEC)[1]))
  for (k in 1:(dim(vinblastine_mod1_ind_samp$Emin_VEC)[1])){
    iter_mat[,k] <- DR_model(dose = 10^seq(-4,-1,0.1), vinblastine_mod1_ind_samp$Emin_VEC[k,j], vinblastine_mod1_ind_samp$Emax_VEC[k,j], vinblastine_mod1_ind_samp$EC50_VEC[k,j], vinblastine_mod1_ind_samp$H_VEC[k,j])
  }
  vinblastine_mod1_ind_list[[j]] <- cbind(data.frame(SUBJECT_IN = (vinblastine_data$SUBJECT_IN[!duplicated(vinblastine_data$PLATE_IN)])[j], PLATE_IN = (vinblastine_data$PLATE_IN[!duplicated(vinblastine_data$PLATE_IN)])[j], dose = 10^seq(-4,-1,0.1)), as.data.frame(rowQuantiles(iter_mat, probs = c(0.05,0.25,0.50,0.75,0.95))))
}
vinblastine_mod1_ind_pred <- list.rbind(vinblastine_mod1_ind_list)
vinblastine_sub_plot <- ggplot(vinblastine_data, aes(x = CONC_uM*1000, y = Difference, shape = DRUG), fill = 'white') + geom_point(aes(color = as.factor(SUBJECT_IN)), alpha = 0.5) + facet_wrap(~PLATE_IN) + scale_shape_manual(values = c(21,19)) + theme_bw() + geom_line(data = vinblastine_mod1_ind_pred, aes(x = dose*1000, y = `50%`, color = as.factor(SUBJECT_IN)), inherit.aes = FALSE, alpha = 0.6) + geom_ribbon(data = vinblastine_mod1_ind_pred, aes(x = dose*1000, ymin = `5%`, ymax = `95%`, fill = as.factor(SUBJECT_IN)), inherit.aes = FALSE, alpha = 0.2) + scale_x_log10() + theme(legend.position = 'bottom', legend.title = element_blank()) + xlab('Vinblastine concentration (nM)') + ylab('Difference in absorption (arbitrary unit)') + coord_cartesian(xlim = c(1,100)) + scale_color_viridis(discrete = TRUE, end = 0.9, aesthetics = c('fill','color'), labels = vinblastine_labels) + guides(shape = 'none')
ggsave(vinblastine_sub_plot, file = 'vinblastine_sub_plot.svg', units = 'mm', width = 200, height = 150)

# Extract and export the individual parameters.
vinblastine_theta_samp <- extract(vinblastine_mod1, pars = c('Emin_pop','Emax_pop','EC50_pop','H_pop','Emax_logit'))
vinblastine_theta_list <- vector(mode = 'list', length = length(vinblastine_theta_samp))
for(k in 1:4){
  vinblastine_theta_list[[k]] <- cbind(data.frame(parameter = (names(vinblastine_theta_samp)[k]), subject = seq(1,4,1)), as.data.frame((colQuantiles(vinblastine_theta_samp[[k]], probs = c(0.05,0.25,0.5,0.75,0.95)))))
}
vinblastine_theta_table <- list.rbind(vinblastine_theta_list)
write.csv(vinblastine_theta_table, file = 'vinblastine_theta_table.csv', row.names = FALSE)

# Visualize the posterior probability distributions for the individual predictions.
vinblastine_pred_list <- vector(mode = 'list', length = ((dim(vinblastine_theta_samp$Emin_pop))[2]))
for(k in 1:((dim(vinblastine_theta_samp$Emin_pop))[2])){
  iter_pred <- matrix(nrow = length(10^seq(-4,-1,0.1)), ncol = ((dim(vinblastine_theta_samp$Emin_pop))[1]))
  for(i in 1:((dim(vinblastine_theta_samp$Emin_pop))[1])){
    iter_pred[,i] <-  DR_model(dose = 10^seq(-4,-1,0.1), vinblastine_theta_samp$Emin_pop[i,k], vinblastine_theta_samp$Emax_pop[i,k], vinblastine_theta_samp$EC50_pop[i,k], vinblastine_theta_samp$H_pop[i,k])
  }
  vinblastine_pred_list[[k]] <- cbind(data.frame(subject = k, dose = 10^seq(-4,-1,0.1)), as.data.frame(rowQuantiles(iter_pred, probs = c(0.05,0.25,0.50,0.75,0.95))))
}
vinblastine_pred_data <- list.rbind(vinblastine_pred_list)
vinblastine_pred_plot <- ggplot(vinblastine_pred_data, aes(x = dose*1000, y = `50%`)) + geom_line(aes(color = as.factor(subject))) + geom_ribbon(aes(ymin = `5%`, ymax = `95%`, fill = as.factor(subject)), alpha = 0.2) + geom_ribbon(aes(ymin = `25%`, ymax = `75%`, fill = as.factor(subject)), alpha = 0.3) + facet_wrap(~subject, labeller = labeller(subject = vinblastine_labels)) + scale_x_log10() + theme_bw() + scale_color_viridis(discrete = TRUE, end = 0.9, aesthetics = c('fill','color')) + theme(panel.grid.minor = element_blank(), legend.position = 'none') + xlab('Vinblastine concentration (nM)') + ylab('Difference in absorption (arbitrary unit)')
ggsave(vinblastine_pred_plot, file = 'vinblastine_pred_plot.svg', units = 'mm', width = 200, height = 200)

# Visualize the posterior probability distributions for the individual pharmacodynamic parameters.
vinblastine_theta_ind_list <- vector(mode = 'list', length = length(vinblastine_theta_samp))
for(i in 1:(length(vinblastine_theta_samp))){
  vinblastine_theta_ind_list[[i]] <- data.frame(subject = rep(c(1,2,3,4), each = (dim((vinblastine_theta_samp[[i]])))[1]), parameter = (names(vinblastine_theta_samp))[i], sample = as.vector(vinblastine_theta_samp[[i]]))
}
vinblastine_theta_ind <- list.rbind(vinblastine_theta_ind_list)
vinblastine_theta_ind$sample[vinblastine_theta_ind$parameter == 'EC50_pop'] <- vinblastine_theta_ind$sample[vinblastine_theta_ind$parameter == 'EC50_pop']*1000
vinblastine_theta_labeller <- as_labeller(c(EC50_pop = "EC[50]~(nM)", Emax_logit = "E[MAX]~(0-1~scale)", Emin_pop = "E[MIN]~(arbitrary~unit)", H_pop = "Slope~italic(h)"), default = label_parsed)
vinblastine_theta_plot <- ggplot(data = vinblastine_theta_ind[!(vinblastine_theta_ind$parameter == 'Emax_pop'),], aes(x = sample)) + stat_slab(aes(fill = as.factor(subject)), normalize = 'panels', alpha = 0.5) + facet_wrap(~parameter, scales = 'free', labeller = labeller(parameter = vinblastine_theta_labeller)) + scale_x_log10() + theme_bw() + scale_y_continuous(breaks = NULL) + scale_fill_viridis(discrete = TRUE, end = 0.9, labels = vinblastine_labels) + theme(legend.position = 'bottom', legend.title = element_blank()) + xlab('') + ylab('')
ggsave(vinblastine_theta_plot, file = 'vinblastine_theta_plot.svg', units = 'mm', width = 200, height = 200)
