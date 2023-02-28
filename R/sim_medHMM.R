#' Simulate data using a multilevel explicit-duration hidden (semi) Markov model
#'
#' \code{sim_medHMM} simulates data for multiple subjects, for which the data
#' have either categorical or continuous (i.e., normally distributed)
#' observations that follow an explicit-duration hidden (semi-) Markov model
#' (HMM) with a multilevel structure. The multilevel structure implies that
#' each subject is allowed to have its own set of parameters, and that the
#' parameters at the subject level (level 1) are tied together by a population
#' distribution at level 2 for each of the corresponding parameters. The shape
#' of the population distribution for each of the parameters is a normal
#' distribution. In addition to (natural and/or unexplained) heterogeneity
#' between subjects, the subjects parameters can also depend on a covariate.
#'
#' In simulating the data, having a multilevel structure means that the
#' parameters for each subject are sampled from the population level
#' distribution of the corresponding parameter. The user specifies the
#' population distribution for each parameter: the average population transition
#' probability matrix and its variance, and the average population emission
#' distribution and its variance. For now, the variance of the mean population
#' parameters is assumed fixed for all components of the transition probability
#' matrix and for all components of the emission distribution.
#'
#' One can simulate multivariate data. That is, the hidden states depend on more
#' than 1 observed variable simultaneously. The distributions of multiple
#' dependent variables for multivariate data are assumed to be independent, and
#' all distributions for one data set have to be of the same type (either
#' categorical or continuous).
#'
#' Note: the subject specific) initial state distributions (i.e., the
#' probability of each of the states at the first time point) needed to simulate
#' the data are obtained from the stationary distributions of the subject
#' specific transition probability matrices gamma.
#'
#'
#' \code{beta}: As the first element in each row of \code{gamma} is used as
#' reference category in the multinomial logistic regression, the first matrix
#' in the list \code{beta} used to predict transition probability matrix
#' \code{gamma} has a number of rows equal to \code{m} and the number of columns
#' equal to \code{m} - 1. The first element in the first row corresponds to the
#' probability of switching from state one to state two. The second element in
#' the first row corresponds to the probability of switching from state one to
#' state three, and so on. The last element in the first row corresponds to the
#' probability of switching from state one to the last state. The same principle
#' holds for the second matrix in the list \code{beta} used to predict
#' categorical emission distribution(s) \code{emiss_distr}: the first element in
#' the first row corresponds to the probability of observing category two in
#' state one. The second element in the first row corresponds to the probability
#' of observing category three is state one, and so on. The last element in the
#' first row corresponds to the probability of observing the last category in
#' state one.
#'
#'
#' @param n_t Numeric vector with length 1 denoting the length of the observed
#'   sequence to be simulated for each subject. To only simulate subject
#'   specific transition probability matrices gamma and emission distributions
#'   (and no data), set \code{t} to 0.
#' @param n Numeric vector with length 1 denoting the number of subjects for
#'   which data is simulated.
#' @param data_distr String vector with length 1 denoting the observation type
#'   of the data to be simulated. Should be set to either \code{'categorical'}
#'   or \code{'continuous'}. Note that when simulating multivariate data, all
#'   dependent variables are assumed to be of the same observation type. The
#'   default equals to \code{data_distr = 'categorical'}.
#' @param m Numeric vector with length 1 denoting the number
#'   of hidden states.
#' @param n_dep Numeric vector with length 1 denoting the
#'   number of dependent variables
#' @param start_state Optional numeric vector with length 1 denoting in which
#'   state the simulated state sequence should start. If left unspecified, the
#'   simulated state for time point 1 is sampled from the initial state
#'   distribution (which is derived from the transition probability matrix
#'   gamma).
#' @param q_emiss Only to be specified if the data to be simulated represents
#'   categorical data. Numeric vector with length \code{n_dep} denoting the
#'   number of observed categories for the categorical emission distribution for
#'   each of the dependent variables.
#' @param gamma A matrix with \code{m} rows and \code{m} columns containing the
#'   average population transition probability matrix used for simulating the
#'   data. That is, the probability to switch from hidden state \emph{i} (row
#'   \emph{i}) to hidden state \emph{j} (column  \emph{j}).
#' @param emiss_distr A list with \code{n_dep} elements containing the average
#'   population emission distribution(s) of the observations given the hidden
#'   states for each of the dependent variables. If \code{data_distr =
#'   'categorical'}, each element is a matrix with \code{m} rows and
#'   \code{q_emiss[k]} columns for each of the \code{k} in \code{n_dep} emission
#'   distribution(s). That is, the probability of observing category \emph{q}
#'   (column \emph{q}) in state \emph{i} (row \emph{i}). If \code{data_distr =
#'   'continuous'}, each element is a matrix with \code{m} rows and 2 columns;
#'   the first column denoting the mean of state \emph{i} (row \emph{i}) and the
#'   second column denoting the variance of state \emph{i} (row \emph{i}) of the
#'   Normal distribution.
#' @param dwell_type Family of distribution to sample state dwelling times.
#'   Currently it only takes de value \code{dwell_type = 'logNormal'} which
#'   specifies an approximation to the discrete dwelling times through
#'   a continuous log-Normal distribution.
#' @param dwell_distr A matrix containing the average population parameters for
#'   the dwelling distribution of the states. If \code{dwell_type =
#'   'logNormal'}, the matrix has a dimension of \code{m} rows and 2 columns;
#'   the first column denoting the log-mean of state \emph{i} (row \emph{i})
#'   and the second column denoting the log-variance of state \emph{i} (row
#'   \emph{i}) of the log-Normal distribution.
#' @param xx_vec List of 1 + \code{n_dep} vectors containing the covariate(s) to
#'   predict the transition probability matrix \code{gamma} and/or (specific)
#'   emission distribution(s) \code{emiss_distr} using the regression parameters
#'   specified in \code{beta} (see below). The first element in the list
#'   \code{xx_vec} is used to predict the transition matrix. Subsequent elements
#'   in the list are used to predict the emission distribution of (each of) the
#'   dependent variable(s). This means that the covariate used to predict
#'   \code{gamma} and \code{emiss_distr} can either be the same covariate,
#'   different covariates, or a covariate for certain elements and none for the
#'   other. At this point, it is only possible to use one covariate for both
#'   \code{gamma} and \code{emiss_distr}. The first vector of the list
#'   \code{xx_vec} is used to predict the transition matrix. The subsequent
#'   vectors in the list \code{xx_vec} are used to predict the emission
#'   distribution(s) of the dependent variable(s). For all elements in the list,
#'   the number of observations in the vectors should be  equal to the number of
#'   subjects to be simulated \code{n}. If \code{xx_vec} is omitted completely,
#'   \code{xx_vec} defaults to NULL, resembling no covariates at all. Specific
#'   elements in the list can also be left empty (i.e., set to \code{NULL}) to
#'   signify that either the transition probability matrix or (one of) the
#'   emission distribution(s) is not predicted by covariates.
#' @param beta List of 1 + \code{n_dep} matrices containing the regression
#'   parameters to predict \code{gamma} and/or \code{emiss_distr} in combination
#'   with \code{xx_vec} using (multinomial logistic) regression. The first
#'   matrix is used to predict the transition probability matrix \code{gamma}.
#'   The subsequent matrices are used to predict the emission distribution(s)
#'   \code{emiss_distr} of the dependent variable(s). For \code{gamma} and
#'   categorical emission distributions, one regression parameter is specified
#'   for each element in \code{gamma} and \code{emiss_distr}, with the following
#'   exception. The first element in each row of \code{gamma} and/or
#'   \code{emiss_distr} is used as reference category in the multinomial
#'   logistic regression. As such, no regression parameters can be specified for
#'   these parameters. Hence, the first element in the list \code{beta} to
#'   predict \code{gamma} consist of a matrix with the number of rows equal to
#'   \code{m} and the number of columns equal to \code{m} - 1. For categorical
#'   emission distributions, the subsequent elements in the list \code{beta} to
#'   predict \code{emiss_distr} consist of a matrix with the number of rows
#'   equal to \code{m} and the number of columns equal to \code{q_emiss[k]} - 1
#'   for each of the \code{k} in \code{n_dep} emission distribution(s). See
#'   \emph{details} for more information. For continuous emission distribuitons,
#'   the subsequent elements in the list \code{beta} consist of a matrix with
#'   the number of rows equal to \code{m} and 1 column.
#'
#'   Note that if \code{beta} is specified, \code{xx_vec} has to be specified as
#'   well. If \code{beta} is omitted completely, \code{beta} defaults to NULL,
#'   resembling no prediction of \code{gamma} and \code{emiss_distr} using
#'   covariates. One of the elements in the list can also be left empty
#'   (i.e., set to \code{NULL}) to signify that either the transition
#'   probability matrix or a specific emission distribution is not predicted by
#'   covariates.
#' @param var_gamma A numeric vector with length 1 denoting the amount of
#'   variance between subjects in the transition probability matrix. Note that
#'   this value corresponds to the variance of the parameters of the
#'   multinomial distribution (i.e., the intercepts of the regression equation
#'   of the multinomial distribution used to sample the transition probability
#'   matrix), see details below. In addition, only one variance value can be
#'   specified for the complete transition probability matrix, hence the
#'   variance is assumed fixed across all components. The default equals 0.1,
#'   which corresponds to little variation between subjects. If one wants to
#'   simulate data from exactly the same HMM for all subjects, var_gamma should
#'   be set to 0. Note that if data for only 1 subject is simulated (i.e., n =
#'   1), \code{var_gamma} is set to 0.
#' @param var_emiss A numeric vector with length \code{n_dep} denoting the
#'   amount of variance between subjects in the emission distribution(s). For
#'   categorical data, this value corresponds to the variance of the parameters
#'   of the multinomial distribution (i.e., the intercepts of the regression
#'   equation of the multinomial distribution used to sample the components of
#'   the emission distribution), see details below.  For continuous data, this
#'   value corresponds to the variance in the mean of the emission
#'   distribution(s) across subjects. Note that only one variance value can be
#'   specified each emission distribution, hence the variance is assumed fixed
#'   across states (and, for the categorical distribution, categories within
#'   a state) within an emission distribution. The default equals 0.1, which
#'   corresponds to little variation between subjects given categorical
#'   observations. If one wants to simulate data from exactly the same HMM for
#'   all subjects, \code{var_emiss} should be set to a vector of 0's. Note that if data
#'   for only 1 subject is simulated (i.e., n = 1), \code{var_emiss} is set to a
#'   vector of 0's.
#' @param var_dwell A numeric vector with length 1 denoting the amount of
#'   variance between subjects in the state dwelling distribution. Note that
#'   in the case of \code{dwell_type = 'logNormal'} this value corresponds to
#'   the variance of the log-mean of the log-Normal distribution in the
#'   logarithmic scale. In addition, only one variance value can be specified
#'   for the complete transition probability matrix, hence the
#'   variance is assumed fixed across all components. The default equals 0.01,
#'   which corresponds to little variation between subjects. If one wants to
#'   simulate data from exactly the same HMM for all subjects, \code{var_dwell}
#'   should be set to 0. Note that if data for only 1 subject is simulated
#'   (i.e., n = 1), \code{var_dwell} is set to 0.
#' @param return_ind_par A logical scalar. Should the subject specific
#'   transition probability matrix \code{gamma} and emission probability matrix
#'   \code{emiss_distr} be returned by the function (\code{return_ind_par =
#'   TRUE}) or not (\code{return_ind_par = FALSE}). The default equals
#'   \code{return_ind_par = FALSE}.
#'
#'
#' @return The following components are returned by the function
#'   \code{sim_medHMM}:
#' \describe{
#'   \item{\code{states}}{A matrix containing the simulated hidden state
#'   sequences, with one row per hidden state per subject. The first column
#'   indicates subject id number. The second column contains the simulated
#'   hidden state sequence, consecutively for all subjects. Hence, the id number
#'   is repeated over the rows (with the number of repeats equal to the length
#'   of the simulated hidden state sequence \code{T} for each subject).}
#'   \item{\code{obs}}{A matrix containing the simulated observed outputs, with
#'   one row per simulated observation per subject. The first column indicates
#'   subject id number. The second column contains the simulated observation
#'   sequence, consecutively for all subjects. Hence, the id number is repeated
#'   over rows (with the number of repeats equal to the length of the simulated
#'   observation sequence \code{T} for each subject).}
#'   \item{\code{subject_gamma}}{A list containing \code{n} elements with the simulated
#'   subject specific transition probability matrices \code{gamma}. Only
#'   returned if \code{return_ind_par} is set to \code{TRUE}.}
#'   \item{\code{subject_emiss}}{A list containing \code{n} elements with the
#'   simulated subject specific emission parameter matrices
#'   \code{emiss_distr}. Only returned if \code{return_ind_par} is set to
#'   \code{TRUE}.}
#'   \item{\code{subject_dwell}}{A list containing \code{n} elements with the
#'   simulated subject specific state dwelling parameter matrices
#'   \code{subject_dwell}. Only returned if \code{return_ind_par} is set to
#'   \code{TRUE}.}
#' }
#'
#'
#'
#' @examples
#' ###### Example on simulated data
#' Simulating multivariate continuous data
#' Define model parameters:
#' n_t <- 200
#' n <- 20
#' m <- 3
#' n_dep <- 2
#'
#'
#' gamma <- matrix(c(0, 0.7, 0.3,
#'                   0.5, 0, 0.5,
#'                   0.6, 0.4, 0), nrow = m, ncol = m, byrow = TRUE)
#'
#' emiss_distr <- list(matrix(c(10,2,
#'                              50,2,
#'                              2,2), nrow = m, ncol = 2, byrow = TRUE),
#'                     matrix(c(-5,2,
#'                              -20,2,
#'                              5,2), nrow = m, ncol = 2, byrow = TRUE))
#'
#' dwell_distr <- dwell_start <- matrix(log(c(20,2,
#'                                            10,2,
#'                                            2,2)), nrow = m, ncol = 2, byrow = TRUE)
#'
#' # Simulate data
#' sim_data <- sim_medHMM(n_t, n, data_distr = 'continuous', m, n_dep = n_dep,
#'                        dwell_distr = dwell_distr, dwell_type = 'logNormal',
#'                        start_state = NULL, q_emiss = NULL, gamma = gamma, emiss_distr = emiss_distr, xx_vec = NULL, beta = NULL,
#'                        var_gamma = 0.1, var_emiss = c(0.1,0.1), var_dwell = 0.01, return_ind_par = TRUE)
#'
#' @export

sim_medHMM <- function(n_t, n, data_distr = 'continuous', m, n_dep = 1,
                      start_state = NULL, q_emiss = NULL, gamma, emiss_distr,
                      dwell_distr, dwell_type = 'logNormal',
                      xx_vec = NULL, beta = NULL,
                      var_gamma = 0.1, var_emiss = NULL, var_dwell = NULL, return_ind_par = FALSE){

    #############
    # Inbuild checks for correct specification of parameters ---------------------
    #############

    if (dim(gamma)[1] != m){
        stop(paste("The transiton probability matrix gamma should be a", m, "by", m, "matrix."))
    }
    if (dim(gamma)[2] != m){
        stop(paste("The transiton probability matrix gamma should be a", m, "by", m, "matrix."))
    }
    if(!isTRUE(all.equal(apply(gamma,1,sum), rep(1,m)))){
        stop("The elements in each row of the transition probability matrix gamma should sum up to 1")
    }
    if(!is.list(emiss_distr)){
        stop("The format of emiss_distr should be a list with", n_dep, "elements.")
    }
    if(length(emiss_distr) != n_dep){
        stop("The number of dependent variables specified in n_dep and the number of elements specified in the list emiss_distr should be equal")
    }
    if(data_distr == "categorical" & length(q_emiss) != n_dep){
        stop("The lenght of q_emiss specifying the number of output categories for each of the number of dependent variables should equal the number of dependent variables specified in n_dep")
    }
    for(i in 1:n_dep){
        if (dim(emiss_distr[[i]])[1] != m){
            stop(paste("The number of rows of emission distribution matrix in element", i, "should be
                       equal to the number of states, which is", m, "."))
        }
        if(data_distr == 'categorical'){
            if (dim(emiss_distr[[i]])[2] != q_emiss[i]){
                stop(paste("The number of columns of the emission distribution matrix should be
                           equal to the number of observable categories, which is", q_emiss[i], ". See emission distribution in element", i, "."))
            }
            if(!isTRUE(all.equal(apply(emiss_distr[[i]], 1, sum), rep(1, m)))){
                stop("The elements in each row of the emission distribution matrix should sum up to 1, see emission distribution in element", i, ".")
            }
        }
        # if(data_distr == 'continuous'){
        #
        # }
    }

    if((is.null(xx_vec) & !is.null(beta)) | (!is.null(xx_vec) & is.null(beta))){
        stop("Either only xx_vec or only beta is specified. Please specify both 1) the values for the covariate
             in xx_vec and 2) the values of the regression parameters in beta, to allow correct simulation of the
             data.")
    }
    if(!is.null(xx_vec)){
        if((!is.null(xx_vec[[1]]) & is.null(beta[[1]])) |
           (!is.null(xx_vec[[2]]) & is.null(beta[[2]]))){
            stop("Either only xx_vec or only beta is specified in one of the elements.
                 Please specify both 1) the values for the covariate in xx_vec and 2)
                 the values of the regression parameters in beta if either one is not
                 empty, to allow correct simulation of the data.")
        }
    }
    if(!is.null(beta)){
        # extend to all 1 + n_dep
        if((!is.null(beta[[1]]) & is.null(xx_vec[[1]])) |
           (!is.null(beta[[2]]) & is.null(xx_vec[[2]]))){
            stop("Either only xx_vec or only beta is specified in one of the elements.
                 Please specify both 1) the values for the covariate in xx_vec and 2)
                 the values of the regression parameters in beta if either one is not
                 empty, to allow correct simulation of the data.")
        }
    }
    if(!is.null(xx_vec)){
        # extend to all 1 + n_dep
        if((!is.null(xx_vec[[1]]) & length(xx_vec[[1]]) != n) |
           (!is.null(xx_vec[[2]]) & length(xx_vec[[2]]) != n)){
            stop("The length of the vectors in xx_vec should be equal to the number of subjects to be simulated,
                 set in n, if (the element in) xx_vec is not set to NULL.")
        }
    }
    if (!is.null(beta)){
        if (!is.null(beta[[1]])){
            if ((dim(beta[[1]])[1] != (m)) | (dim(beta[[1]])[2] != (m-1))){
                stop(paste("The first element of beta to predict the transiton probability matrix gamma should be a m (", m, " ) by m - 1 (", m - 1, ") matrix."))
            }
        }
        if (!is.null(beta[[2]]) & data_distr == 'categorical'){
            # extend to all 1 + n_dep and continuous
            if((dim(beta[[2]])[1] != (m)) | (dim(beta[[2]])[2] != (q_emiss[1]-1))){
                stop(paste("The second element of beta to predict the emission distribution should be a m (", m, ") by q_emiss - 1 (", q_emiss[1] - 1, ") matrix."))
            }
        }
    }
    if(is.null(xx_vec)){
        xx_vec <- rep(list(NULL), n_dep + 1)
        for(i in 1:(n_dep + 1)){
            xx_vec[[i]] <- rep(1,n)
        }
    } else {
        for(i in 1:(n_dep + 1)){
            if(is.null(xx_vec[[i]])) {
                xx_vec[[i]] <- rep(1,n)
            }
        }
    }
    if(is.null(beta)){
        beta <- rep(list(NULL), n_dep + 1)
        beta[[1]] <- matrix(0, ncol = m - 1, nrow = m)
        for(i in 2:(n_dep + 1)){
            if(data_distr == 'categorical'){
                beta[[i]] <- matrix(0, ncol = q_emiss[i-1] - 1, nrow = m)
            } else if (data_distr == 'continuous'){
                beta[[i]] <- matrix(0, ncol = 1, nrow = m)
            }
        }
    } else {
        if(is.null(beta[[1]])) {
            beta[[1]] <- matrix(0, ncol = m - 1, nrow = m)
        }
        for (i in 2:(n_dep + 1)){
            if (is.null(beta[[i]])) {
                if(data_distr == 'categorical'){
                    beta[[i]] <- matrix(0, ncol = q_emiss[i-1] - 1, nrow = m)
                } else if (data_distr == 'continuous'){
                    beta[[i]] <- matrix(0, ncol = 1, nrow = m)
                }
            }
        }
    }

    if(n == 1){
        var_gamma <- 0
        var_emiss <- rep(0, n_dep)
    }
    if(is.null(var_emiss)){
        var_emiss <- rep(0.1, n_dep)
    } else if(length(var_emiss) != n_dep){
        stop("The lenght of var_emiss specifying variance between subjects in each of the number of dependent variables should equal the number of dependent variables specified in n_dep")
    }
    #----------------------------------------------------------------------#
    if(is.null(var_dwell)){
        var_dwell <- 0.01
    }
    #----------------------------------------------------------------------#

    #############
    # Simulating the data ---------------------
    #############

    states <- matrix(ncol = 2, nrow = n_t*n)
    states[,1] <- rep(1:n, each = n_t)
    obs <- matrix(ncol = 1 + n_dep, nrow = n_t*n)
    obs[,1] <- rep(1:n, each = n_t)
    sub_gamma <- rep(list(NULL), n)
    sub_emiss <- rep(list(vector("list", n_dep)), n)
    #----------------------------------------------------------------------#
    sub_durat <- rep(list(NULL), n)
    sub_gamma_full <- rep(list(NULL), n)
    #----------------------------------------------------------------------#
    mnl_gamma <- prob_to_int(gamma)
    if(data_distr == "categorical"){
        mnl_emiss <- rep(list(NULL), n_dep)
        for(i in 1:n_dep){
            mnl_emiss[[i]] <- prob_to_int(emiss_distr[[i]])
        }
    }
    for(j in 1:n){
        sub_gamma[[j]] <- int_to_prob(mnl_gamma + xx_vec[[1]][j] * beta[[1]] +
                                          rnorm(n = m * (m-1), mean = 0, sd = sqrt(var_gamma)))
        #----------------------------------------------------------------------#
        # Modified gamma
        gamma <- t(sub_gamma[[j]])
        sub_gamma[[j]] <- matrix(gamma[upper.tri(gamma) | lower.tri(gamma)], nrow = m, ncol = m-1, byrow = TRUE)

        if(dwell_type == "logNormal"){
            sub_durat[[j]] <- dwell_distr
            sub_durat[[j]][,1] <- dwell_distr[,1] + rnorm(n = m, mean = 0, sd = sqrt(var_dwell))
            # sub_durat[[j]][,1] <- max(dwell_distr[,1] + rnorm(n = m, mean = 0, sd = sqrt(var_dwell)),1)
        }

        #----------------------------------------------------------------------#

        for(i in 1:n_dep){
            if(data_distr == "categorical"){
                sub_emiss[[j]][[i]] <- int_to_prob(mnl_emiss[[i]] + xx_vec[[1+i]][j] * beta[[1+i]] +
                                                       rnorm(n = m * (q_emiss[i]-1), mean = 0, sd = sqrt(var_emiss[i])))
            } else if(data_distr == "continuous"){
                sub_emiss[[j]][[i]] <- emiss_distr[[i]]
                sub_emiss[[j]][[i]][,1] <- emiss_distr[[i]][,1] +  xx_vec[[1+i]][j] * beta[[1+i]] +
                    rnorm(n = m, mean = 0, sd = sqrt(var_emiss[i]))
            }
        }

        if(n_t != 0){

            sub_gamma_full[[j]] <- matrix(0, nrow = m, ncol = m)
            sub_gamma_full[[j]][upper.tri(sub_gamma_full[[j]]) | lower.tri(sub_gamma_full[[j]])] <- t(sub_gamma[[j]])
            sub_gamma_full[[j]] <- t(sub_gamma_full[[j]])

            init <- solve(t(diag(m) - sub_gamma_full[[j]] + 1), rep(1, m))
            if (is.null(start_state)){
                states[((j-1) * n_t + 1), 2] <- sample(x = 1:m, size = 1, prob = init)
            } else {
                states[((j-1) * n_t + 1), 2] <- start_state
            }
            if(data_distr == "categorical"){
                for(i in 1:n_dep){
                    obs[((j-1) * n_t + 1), (1+i)] <- sample(x = 1:q_emiss[i], size = 1, prob = sub_emiss[[j]][[i]][states[((j-1) * n_t + 1), 2],])
                }
            } else if (data_distr == "continuous"){
                for(i in 1:n_dep){
                    obs[((j-1) * n_t + 1), (1+i)] <- rnorm(1, mean = sub_emiss[[j]][[i]][states[((j-1) * n_t + 1), 2],1], sd = sqrt(sub_emiss[[j]][[i]][states[((j-1) * n_t + 1), 2],2]))
                }
            }
            for(t in 2:n_t){

                #--------------------------------------------------------------#
                # If the position has not been covered yet because the state finished, sample a state; else, if the state continues, keep it
                if(is.na(states[((j-1) * n_t + t), 2])){
                    # Sample state
                    allowed_states <- (1:m)[-states[((j-1) * n_t + t - 1), 2]]
                    states[((j-1) * n_t + t), 2] <- sample(x = allowed_states, size = 1, prob = sub_gamma[[j]][states[((j-1) * n_t + t - 1), 2],])

                    # Sample state duration
                    durat <- round(rlnorm(1, meanlog = max(sub_durat[[j]][states[((j-1) * n_t + t), 2],1],1), sdlog = sqrt(sub_durat[[j]][states[((j-1) * n_t + t), 2],2])),0)

                    # Paste the state over the min(duration sampled, remaining timesteps)
                    states[((j-1) * n_t + t):((j-1) * n_t + t + min(durat,(n_t-t))), 2] <- states[((j-1) * n_t + t), 2]

                }
                #--------------------------------------------------------------#

                if(data_distr == "categorical"){
                    for(i in 1:n_dep){
                        obs[((j-1) * n_t + t), (1+i)] <- sample(x = 1:q_emiss[i], size = 1, prob = sub_emiss[[j]][[i]][states[((j-1) * n_t + t), 2],])
                    }
                } else if (data_distr == "continuous"){
                    for(i in 1:n_dep){
                        obs[((j-1) * n_t + t), (1+i)] <- rnorm(1, mean = sub_emiss[[j]][[i]][states[((j-1) * n_t + t), 2],1], sd = sqrt(sub_emiss[[j]][[i]][states[((j-1) * n_t + t), 2],2]))
                    }
                }
            }
        }
    }

    #############
    # Returning output  ---------------------
    #############
    colnames(states) <- c("subj", "state")
    colnames(obs)    <- c("subj", paste("observation", 1:n_dep))
    if (return_ind_par == FALSE & n_t != 0){
        return(list(states = states, obs = obs))
    } else if (return_ind_par == TRUE & n_t != 0){
        return(list(states = states, obs = obs, subject_gamma = sub_gamma_full, subject_emiss = sub_emiss, subject_dwell = sub_durat))
    } else if (n_t == 0){
        return(list(subject_gamma = sub_gamma, subject_emiss = sub_emiss, subject_dwell = sub_durat))
    }
}
