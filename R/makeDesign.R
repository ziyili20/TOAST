makeDesign <- function(design, Prop) {

    # design is a N by P matrix,
    #       with rows as samples and columns as phenotypes
    # Prop is a N by K matrix,
    #     with rows as samples and columns as cell types

    if (!is.vector(Prop)) {
        K <- ncol(Prop)
        
        if (is.null(colnames(Prop))) {
            colnames(Prop) <- paste0("celltype", seq_len(K))
        }
    } else {
        K = 1
    }
    
    dd <- cbind(Prop, design)

    ## make a formula
    formul <- paste("~", paste(colnames(dd)[seq_len(K)], collapse="+"))
    for(i in (K+1):ncol(dd)) {
       ## interaction terms for this factor
       tmp <- paste(colnames(dd)[seq_len(K)], colnames(dd)[i], sep=":")
       intterms <- paste(tmp, collapse="+")
       formul <- paste(formul, intterms, sep="+")
    }
    if(all(design==0)) {
        design_matrix <- Prop[,-1]
    } else {
        design_matrix <- model.matrix(as.formula(formul), dd)[,-1]
    }

    ## remove columns that are all zero
    design_matrix <- design_matrix[,!colSums(design_matrix==0) == nrow(Prop)]
    
    design_matrix <- unique(design_matrix, MARGIN = 2)
    
    formul <- paste0("~ ", paste(colnames(design_matrix), sep="+", collapse = "+"))

    return(list(design_matrix = design_matrix,
             Prop = Prop,
             design = design,
             all_coefs = colnames(design),
             all_cell_types = colnames(Prop),
             formula = formul))
}
