GetPrior <- function(alpha = NULL, sigma = NULL) {

     currentTis <- c("human pbmc",
                     "human brain",
                     "human pancreas",
                     "human liver",
                     "human skin")

     if (is.character(alpha)) {
          l_alpha <- tolower(alpha)
          if (l_alpha %in% currentTis) {
               if (l_alpha == "human pbmc"){
                    alpha_prior <- c(0.09475, 0.23471, 0.33232, 0.0969, 0.24132)
                    sigma_prior <- c(0.09963, 0.14418, 0.16024, 0.10064, 0.14556)
                    names(alpha_prior) <- c("B cells", "CD8T", "CD4T",
                                            "NK cells", "Monocytes")
                    names(sigma_prior) <- names(alpha_prior)
               } else if (l_alpha == "human brain") {
                    alpha_prior <- c(0.04956, 0.00181, 0.49782, 0.13116, 0.02906,
                                     0.025208, 0.03572, 0.00279)
                    sigma_prior <- c(0.03113, 0.002348, 0.12279, 0.04686, 0.02136,
                                     0.11595, 0.01447, 0.00336)
                    names(alpha_prior) <- c("Astrocyte", "Endothelial",
                                            "Excitatory neuron",
                                            "Inhibitory neuron",
                                            "Microglia", "Oligodendrocyte",
                                            "Progenitor", "Pericyte")
                    names(sigma_prior) <- names(alpha_prior)
               } else if (l_alpha == "human pancreas") {
                    alpha_prior <- c(0.35888, 0.01054, 0.08185, 0.18519, 0.12705,
                                     0.13218, 0.04619, 0.04187, 0.00122, 0.0107)
                    sigma_prior <- c(0.08864, 0.00585, 0.04031, 0.07681, 0.07122,
                                     0.12007, 0.04223, 0.02737, 0.00165, 0.01734)
                    names(alpha_prior) <- c("Alpha", "Endothelial",
                                            "Delta", "Beta", "Duct",
                                            "Acinar", "PP", "Mesenchymal",
                                            "Epsilon", "Unclear")
                    names(sigma_prior) <- names(alpha_prior)
               } else if (l_alpha == "human liver") {
                    alpha_prior <- c(0.11457, 0.01384, 0.16225, 0.29166, 0.30844,
                                     0.09445, 0.00369, 0.0111)
                    sigma_prior <- c(0.059, 0.00547, 0.07945, 0.26653, 0.355869,
                                     0.1141, 0.00362, 0.01599)
                    names(alpha_prior) <- c("Endothelial", "Cholangiocyte",
                                            "Kupffer", "T/NK cell","Hepatocyte",
                                            "B cells", "Stellate", "Erythroid cell")
                    names(sigma_prior) <- names(alpha_prior)
               } else if (l_alpha == "human skin") {
                    alpha_prior <- c(0.1753, 0.1175, 0.3122, 0.04494, 0.01165,
                                     0.07953, 0.03754, 0.2213)
                    sigma_prior <- c(0.05884, 0.05909, 0.1408, 0.01456, 0.005914,
                                     0.04389, 0.00699, 0.08916)
                    names(alpha_prior) <- c("Basal keratinocytes", "Endothelial",
                                            "Fibroblasts", "Macrophages/DC",
                                            "Melanocytes", "Pericytes",
                                            "Smooth muscle", "Suprabasal Keratinocytes")
                    names(sigma_prior) <- names(alpha_prior)
               }
          } else {
               message("Supported tissue types: ")
               message(paste(currentTis, collapse="; "))
               stop("Entered tissue is not supported. Please manually specify alpha and sigma!")
          }
     } else {
          alpha_prior = alpha
          sigma_prior = sigma
     }
     return(list(alpha_prior = alpha_prior,
                 sigma_prior = sigma_prior))
}
