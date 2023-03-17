#' Calculates the differentially expressed multivariate gene (DEMVG) score
#' @export
#' @param x DESeq normalized and vst transformed matrix
#' @param t_idx The column index of treatment group
#' @param c_idx The column index of control group

DEMVG<- function(x,t_idx,c_idx) {
  stopifnot("Input x must be DESeq normalized and vst treated matrix" = is.matrix(x), ncol(x)>=2, nrow(x)>=2)
  stopifnot("Column index should be a integer vector" = is.integer(t_idx),is.integer(c_idx),all(c(t_idx,c_idx)>0))
  col_space_x<-x[,c_idx]  # Control group
  col_space_y<-x[,t_idx]  # Treatment group

  # projection (columnar space)
  proj_fit_y_to_x<-col_space_x%*%MASS::ginv(t(col_space_x)%*%col_space_x)%*%t(col_space_x)%*%col_space_y
  proj_fit_x_to_y<-col_space_y%*%MASS::ginv(t(col_space_y)%*%col_space_y)%*%t(col_space_y)%*%col_space_x

  # orthogonal space
  proj_residual_y_to_x<-col_space_y - proj_fit_y_to_x
  proj_residual_x_to_y<-col_space_x - proj_fit_x_to_y

  proj_residual_merged_ortho_space<-cbind( ncol(proj_residual_x_to_y)*proj_residual_y_to_x, (-1)*ncol(proj_residual_y_to_x)*proj_residual_x_to_y)
  proj_residual_merged_ortho_space_prcomp<-stats::prcomp(proj_residual_merged_ortho_space,center=FALSE,scale. = FALSE)
  proj_residual_merged_ortho_space_DEMVG<-t(t(proj_residual_merged_ortho_space_prcomp$x[,1,drop=F]))* sign(mean(proj_residual_merged_ortho_space_prcomp$rotation[,1]))

  colnames(proj_residual_merged_ortho_space_DEMVG)<-"DEMVG_score"
  rownames(proj_residual_merged_ortho_space_DEMVG)<-rownames(x)
  proj_residual_merged_ortho_space_DEMVG
}
