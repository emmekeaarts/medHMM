# fractional log likelihood for multinomial intercept only model
	llmnl_int_frac <- function(int, Obs, n_cat, pooled.likel, w, wgt){
		return((1-w)*llmnl_int(int = int, Obs = Obs, n_cat = n_cat) + w* wgt* pooled.likel)		
	}
