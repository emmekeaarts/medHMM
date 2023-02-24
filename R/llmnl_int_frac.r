#' @keywords internal
# fractional log likelihood for multinomial intercept only model
	llmnl_int_frac <- function(int, Obs, pooled.likel, w, wgt){
		return((1-w)*llmnl_int(int = int, Obs = Obs) + w* wgt* pooled.likel)
	}
