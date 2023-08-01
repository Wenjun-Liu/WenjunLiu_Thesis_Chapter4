spia_my <- function(de=NULL,BminsI = NULL){

  # remove pathway with 0 DEG in it
  kgToRemove <- lapply(names(BminsI), function(x){
    pathwayGene <- rownames(BminsI[[x]])
    length(intersect(names(de),
                     pathwayGene)) != 0
  }) %>%
    unlist()
  BminsI <- BminsI[kgToRemove]
  if(length(BminsI) == 0){
    tA <- NULL
    tA
  } else {
    # Extract the logFCs of DEGs contained in each pathway and assign nonDEGs in each pathway logFC of 0.
    delE <- sapply(names(BminsI), function(x){
      pathwayGene <- rownames(BminsI[[x]])
      de[pathwayGene]%>%
        replace(is.na(.),0)
    }, simplify = FALSE)
    PF <- sapply(names(delE), function(x){
      # if BminsI is invertible, solve for (B-I)PF = -delE
      if (det(BminsI[[x]]) != 0){
        pf <- solve(BminsI[[x]], -delE[[x]])
      } else {
        pf <- NULL
      }
    }, simplify = FALSE) %>%
      .[lapply(.,length) > 0]

    tA <- sapply(names(PF), function(x){
      sum(PF[[x]] - delE[[x]])
    })
    if (length(tA) == 0){
      tA <- NULL
      tA
    } else{
      tA %>%
        as.data.frame() %>%
        set_colnames("tA") %>%
        rownames_to_column("gs_name")
    }
  }

}
