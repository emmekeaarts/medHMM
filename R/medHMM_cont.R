#' Multilevel explicit duration hidden Markov model using Bayesian estimation
#' for continuous observations
#'
#' \code{medHMM_cont} fits a multilevel (also known as mixed or random effects)
#' explicit duration
#' hidden (semi-) Markov model (HMM) to intense longitudinal data with continuous
#' observations (i.e., normally distributed) of multiple subjects using Bayesian
#' estimation, and creates an object of class mHMM_cont. By using a multilevel
#' framework, we allow for heterogeneity in the model parameters between
#' subjects, while estimating one overall HMM. The function includes the
#' possibility to add covariates at level 2 (i.e., at the subject level) and
#' have varying observation lengths over subjects. For a short description of
#' the package see \link{mHMMbayes}. See \code{vignette("tutorial-mhmm")} for an
#' introduction to multilevel hidden Markov models and the package, and see
#' \code{vignette("estimation-mhmm")} for an overview of the used estimation
#' algorithms.
#'
#' Covariates specified in \code{xx} can either be dichotomous or continuous
#' variables. Dichotomous variables have to be coded as 0/1 variables.
#' Categorical or factor variables can as yet not be used as predictor
#' covariates. The user can however break up the categorical variable in
#' multiple dummy variables (i.e., dichotomous variables), which can be used
#' simultaneously in the analysis. Continuous predictors are automatically
#' centered. That is, the mean value of the covariate is subtracted from all
#' values of the covariate such that the new mean equals zero. This is done such
#' that the presented probabilities in the output (i.e., for the population
#' transition probability matrix and population emission probabilities)
#' correspond to the predicted probabilities at the average value of the
#' covariate(s).
#'
#' If covariates are specified and the user wants to set the values for the
#' parameters of the hyper-prior distributions manually, the specification of
#' the elements in the arguments of the hyper-prior parameter values of the
#' normal distribution on the means change as follows: the number of rows in the
#' matrices \code{gamma_mu0} and \code{emiss_mu0} are equal to 1 + the number of
#' covariates used to predict the transition probability matrix for
#' \code{gamma_mu0} and the emission distribution for \code{emiss_mu0} (i.e.,
#' the first row correspond to the hyper-prior mean values of the intercepts,
#' the subsequent rows correspond to the hyper-prior mean values of the
#' regression coefficients connected to each of the covariates), and
#' \code{gamma_K0} and \code{emiss_K0} are now a matrix with the number of
#' hypothetical prior subjects on which the vectors of means (i.e., the rows in
#' \code{gamma_mu0} or \code{emiss_mu0}) are based on the diagonal, and
#' off-diagonal elements equal to 0. Note that the hyper-prior parameter values
#' of the inverse Wishart distribution on the covariance matrix remains
#' unchanged, as the estimates of the regression coefficients for the covariates
#' are fixed over subjects.
#'
#' @param s_data A matrix containing the observations to be modelled, where the
#'   rows represent the observations over time. In \code{s_data}, the first
#'   column indicates subject id number. Hence, the id number is repeated over
#'   rows equal to the number of observations for that subject. The subsequent
#'   columns contain the dependent variable(s). Note that the dependent
#'   variables are assumed to be continuous (i.e., normally distributed
#'   depending on the hidden states). The total number of rows are equal to the
#'   sum over the number of observations of each subject, and the number of
#'   columns are equal to the number of dependent variables (\code{n_dep}) + 1.
#'   The number of observations can vary over subjects.
#' @param gen List containing the following elements denoting the general model
#'   properties:
#'   \itemize{\item{\code{m}: numeric vector with length 1 denoting the number
#'   of hidden states}
#'   \item{\code{n_dep}: numeric vector with length 1 denoting the
#'   number of dependent variables}}
#' @param xx An optional list of (level 2) covariates to predict the transition
#'   matrix and/or the emission probabilities. Level 2 covariate(s) means that
#'   there is one observation per subject of each covariate. The first element
#'   in the list \code{xx} is used to predict the transition matrix. Subsequent
#'   elements in the list are used to predict the emission distribution of (each
#'   of) the dependent variable(s). Each element in the list is a matrix, with
#'   the number of rows equal to the number of subjects. The first column of
#'   each matrix represents the intercept, that is, a column only consisting of
#'   ones. Subsequent columns correspond to covariates used to predict the
#'   transition matrix / emission distribution. See \emph{Details} for more
#'   information on the use of covariates.
#'
#'   If \code{xx} is omitted completely, \code{xx} defaults to \code{NULL},
#'   resembling no covariates. Specific elements in the list can also be left
#'   empty (i.e., set to \code{NULL}) to signify that either the transition
#'   probability matrix or a specific emission distribution is not predicted by
#'   covariates.
#' @param start_val List containing the start values for the transition
#'   probability matrix gamma, the emission distribution(s), and the dwell time
#'   distributions. The first
#'   element of the list contains a \code{m} by \code{m} matrix with the start
#'   values for gamma. The subsequent \code{n_dep} elements are matrices with
#'   \code{m} rows
#'   and 2 columns; the first column denoting the mean of state \emph{i} (row
#'   \emph{i}) and the second column denoting the variance of state \emph{i}
#'   (row \emph{i}) of the Normal distribution. The last element contains a
#'   matrix of \code{m} rows and 2 columns with the start values for the dwell time
#'   distribution. Note that \code{start_val} should not contain nested lists
#'   (i.e., lists within lists).
#' @param mcmc List of Markov chain Monte Carlo (MCMC) arguments, containing the
#'   following elements:
#'   \itemize{\item{\code{J}: numeric vector with length 1 denoting the number
#'   of iterations of the MCMC algorithm}
#'   \item{\code{burn_in}: numeric vector with length 1 denoting the
#'   burn-in period for the MCMC algorithm.}}
#' @param return_path A logical scalar. Should the sampled state sequence
#'   obtained at each iteration and for each subject be returned by the function
#'   (\code{sample_path = TRUE}) or not (\code{sample_path = FALSE}). Note that
#'   the sampled state sequence is quite a large object, hence the default
#'   setting is \code{sample_path = FALSE}. Can be used for local decoding
#'   purposes.
#' @param print_iter The argument print_iter is depricated; please use
#'   show_progress instead to show the progress of the algorithm.
#' @param show_progress A logical scaler. Should the function show a text
#'   progress bar in the \code{R} console to represent the progress of the
#'   algorithm (\code{show_progress = TRUE}) or not (\code{show_progress =
#'   FALSE}). Defaults to \code{show_progress = TRUE}.
#' @param gamma_hyp_prior An optional list containing user specified parameters
#'  of the hyper-prior distribution on the multivariate normal distribution
#'  of the intercepts (and regression coefficients given that covariates are
#'  used) of the multinomial regression model of the transition probability
#'  matrix gamma. The hyper-prior of the mean intercepts is a multivariate
#'  Normal distribution, the hyper-prior of the covariance matrix between the
#'  set of (state specific) intercepts is an Inverse Wishart distribution.
#'
#'  Hence, the list \code{gamma_hyp_prior} contains the following elements:
#'  \itemize{\item{\code{gamma_mu0}: a list containing m matrices; one matrix
#'  for each row of the transition probability matrix gamma. Each matrix
#'  contains the hypothesized mean values of the intercepts. Hence, each matrix
#'  consists of one row (when not including covariates in the model) and
#'  \code{m} - 1 columns}
#'  \item{\code{gamma_K0}: numeric vector with length 1 denoting the number of
#'  hypothetical prior subjects on which the vector of means \code{gamma_mu0} is
#'  based}
#'  \item{\code{gamma_nu}: numeric vector with length 1 denoting the degrees of
#'  freedom of the Inverse Wishart distribution}
#'  \item{\code{gamma_V}: matrix of \code{m} - 1 by \code{m} - 1 containing the
#'  hypothesized variance-covariance matrix between the set of intercepts.}}
#'  Note that \code{gamma_K0}, \code{gamma_nu} and \code{gamma_V} are assumed
#'  equal over the states. The mean values of the intercepts (and regression
#'  coefficients of the covariates) denoted by \code{gamma_mu0} are allowed to
#'  vary over the states.
#'
#'  The default values for the hyper-prior on gamma are: all elements of the
#'  matrices contained in \code{gamma_mu0} set to 0, \code{gamma_K0} set to 1,
#'  \code{gamma_nu} set to 3 + m - 1, and the diagonal of \code{gamma_V} (i.e.,
#'  the variance) set to 3 + m - 1 and the off-diagonal elements (i.e., the
#'  covariance) set to 0.
#'
#'  See \emph{Details} below if covariates are used for changes in the settings
#'  of the arguments of \code{gamma_hyp_prior}.
#' @param emiss_hyp_prior A list containing user specified parameters
#'   of the hyper-prior distribution on the Normal (i.e., Gaussian) emission
#'   distributions (and regression coefficients given that covariates are used)
#'   for each of the states. The hyper-prior connected to the means of the
#'   Normal emission distribution(s) is a Normal-Inverse-Gamma distribution
#'   (i.e., assuming both unknown populaltion mean and variance between subject
#'   level means). The hyper-prior on each of fixed variances of the Normal
#'   emission distribuitons is an Inverse gamma distribution (i.e., assuming a
#'   known mean).
#'
#'  Hence, the list \code{emiss_hyp_prior} contains the following elements:
#'  \itemize{\item{\code{emiss_mu0}: a list containing \code{n_dep} matrices
#'  with one row (when not including covariates in the model) and \code{m}
#'  columns denoting the hypothesized mean values of the Normal emission
#'  distributions in each of the states for each dependent variable \code{k}}.
#'  \item{\code{emiss_K0}: a list containing \code{n_dep} elements corresponding
#'  to each of the dependent variables, where each element is an integer denoting
#'  the number of hypothetical prior subjects on which the vector of means
#'  \code{emiss_mu0} is based}.
#'  \item{\code{emiss_nu}: a list containing \code{n_dep} elements corresponding
#'  to each of the dependent variables, where each element is an integer
#'  denoting the degrees of freedom of the Inverse Gamma hyper-prior
#'  distribution connected to the emission distribution means (note: here, the
#'  Inverse Gamma hyper-prior distribution is parametrized as a scaled inverse
#'  chi-squared distribution).}
#'  \item{\code{emiss_V}: a list containing \code{n_dep} elements corresponding
#'  to each of the dependent variables \code{k}, where each element is a vector
#'  with lenght \code{m} containing the hypothesized variances between the
#'  between subject level means of the Inverse Gamma hyper-prior distribution
#'  connected to the emission distribution means (note: here, the Inverse Gamma
#'  hyper-prior distribution is parametrized as a scaled inverse chi-squared
#'  distribution).}
#'  \item{\code{emiss_a0}: a list containing \code{n_dep} elements corresponding
#'  to each of the dependent variables \code{k}, where each element is a vector
#'  with lenght \code{m} containing the shape values of the Inverse Gamma
#'  hyper-prior on each of fixed variances of the Normal emission distribuitons
#'  (note: here the standard Inverse Gamma parametrization is used).}
#'  \item{\code{emiss_b0}: a list containing \code{n_dep} elements corresponding
#'  to each of the dependent variables \code{k}, where each element is a vector
#'  with lenght \code{m} containing the scale values of the Inverse Gamma
#'  hyper-prior on each of fixed variances of the Normal emission distribuitons
#'  (note: here the standard Inverse Gamma parametrization is used).}}
#'  Note that \code{emiss_K0} and \code{emiss_nu} are assumed
#'  equal over the states.
#' @param dwell_hyp_prior A list containing user specified parameters
#'   of the hyper-prior distribution on the log-Normal (i.e., log-Gaussian)
#'   dwell time distribution for each of the states. The hyper-prior connected
#'   to the log-means of the log-Normal state dwell time distribution(s) is a
#'   Normal-Inverse-Gamma distribution (i.e., assuming both unknown population
#'   mean and variance between subject level means). The hyper-prior on each of
#'   fixed variances of the Normal emission distributions is an Inverse gamma
#'   distribution (i.e., assuming a known mean). Note that the values used as
#'   log-mean priors are to be entered in the logarithmic scale.
#'
#'  Hence, the list \code{dwell_hyp_prior} contains the following elements:
#'  \itemize{\item{\code{d_mu0}: a vector of \code{m} elements denoting the
#'  hypothesized mean values of the log-Normal dwell distribution for each of
#'  the states in the logarithmic scale}.
#'  \item{\code{s2_0}: a vector with lenght \code{m} containing the
#'  hypothesized variances between the between subject level means of the
#'  Inverse Gamma hyper-prior distribution connected to the dwell distribution
#'  log-means (note: here, the Inverse Gamma hyper-prior distribution is
#'  parametrized as a scaled inverse chi-squared distribution).}
#'  \item{\code{alpha.sigma20}: a vector with lenght \code{m} containing the
#'  shape values of the Inverse Gamma hyper-prior on each of fixed variances of
#'  the log-Normal dwell distribution (note: here the standard Inverse Gamma
#'  parametrization is used).}
#'  \item{\code{beta.sigma20}: a vector with lenght \code{m} containing the
#'  scale values of the Inverse Gamma hyper-prior on each of fixed variances of
#'  the log-Normal dwell distribution (note: here the standard Inverse Gamma
#'  parametrization is used).}
#'  \item{\code{alpha.tau20}: a vector with lenght \code{m} containing the
#'  shape values of the Inverse Gamma hyper-prior on each of the between subject
#'  variances of the Normal prior on the dwell distribution (note: here the
#'  standard Inverse Gamma parametrization is used).}
#'  \item{\code{beta.tau20}: a vector with lenght \code{m} containing the
#'  scale values of the Inverse Gamma hyper-prior on each of the between subject
#'  variances of the Normal prior on the dwell distribution (note: here the
#'  standard Inverse Gamma parametrization is used).}}
#'  Note that \code{emiss_K0} and \code{emiss_nu} are assumed
#'  equal over the states.
#'
#'  See \emph{Details} below if covariates are used for changes in the settings
#'  of the arguments of \code{emiss_hyp_prior}.
#' @param gamma_sampler An optional list containing user specified settings for
#'   the proposal distribution of the random walk (RW) Metropolis sampler for
#'   the subject level parameter estimates of the intercepts modeling the
#'   transition probability matrix. The list \code{gamma_sampler} contains the
#'   following elements:
#'  \itemize{\item{\code{gamma_int_mle0}: a numeric vector with length \code{m}
#'  - 1 denoting the start values for the maximum likelihood estimates of the
#'  intercepts for the transition probability matrix gamma, based on the pooled
#'  data (data over all subjects)}
#'  \item{\code{gamma_scalar}: a numeric vector with length 1 denoting the scale
#'  factor \code{s}. That is, The scale of the proposal distribution is composed
#'  of a covariance matrix Sigma, which is then tuned by multiplying it by a
#'  scaling factor \code{s}^2}
#'  \item{\code{gamma_w}: a numeric vector with length 1 denoting the weight for
#'  the overall log likelihood (i.e., log likelihood based on the pooled data
#'  over all subjects) in the fractional likelihood.}}
#'   Default settings are: all elements in \code{gamma_int_mle0} set to 0,
#'   \code{gamma_scalar} set to 2.93 / sqrt(\code{m} - 1), and \code{gamma_w} set to
#'   0.1. See the section \emph{Scaling the proposal distribution of the RW
#'   Metropolis sampler} in \code{vignette("estimation-mhmm")} for details.
#' @param max_dwell An optional value for the maximal dwell time in discrete
#'   time steps allowed for the states. Lower values in \code{max_dwell} improve
#'   the speed of the algorithm at the risk of underestimating the duration of
#'   states. When omitted, it defaults to the number of occasions on each
#'   subject (that os, all possible durations are assessed).
#'
#' @return \code{mHMM_cont} returns an object of class \code{mHMM_cont}, which has
#'   \code{print} and \code{summary} methods to see the results.
#'   The object contains the following components:
#'   \describe{
#'   \item{\code{PD_subj}}{A list containing one matrix per subject with the
#'   subject level parameter estimates and the log likelihood over the
#'   iterations of the hybrid Metropolis within Gibbs sampler. The iterations of
#'   the sampler are contained in the rows, and the columns contain the subject
#'   level estimates of subsequently the emission means, the (fixed over subjects)
#'   emission variances, the transition probabilities, the dwell time log-means,
#'   the (fixed over subjects) dwell log-variances, and the log likelihood.}
#'   \item{\code{gamma_prob_bar}}{A matrix containing the group level parameter
#'   estimates of the transition probabilities over the iterations of the hybrid
#'   Metropolis within Gibbs sampler. The iterations of the sampler are
#'   contained in the rows, and the columns contain the group level parameter
#'   estimates. If covariates were included in the analysis, the group level
#'   probabilities represent the predicted probability given that the covariate
#'   is at the average value for continuous covariates, or given that the
#'   covariate equals zero for dichotomous covariates.}
#'   \item{\code{gamma_int_bar}}{A matrix containing the group level intercepts
#'   of the multinomial logistic regression modeling the transition
#'   probabilities over the iterations of the hybrid Metropolis within Gibbs
#'   sampler. The iterations of the sampler are contained in the rows, and the
#'   columns contain the group level intercepts.}
#'   \item{\code{gamma_V_int_bar}}{A matrix containing the (co-)variance
#'   components for the subject-level intercepts
#'   of the multinomial logistic regression modeling the transition
#'   probabilities over the iterations of the hybrid Metropolis within Gibbs
#'   sampler. The iterations of the sampler are contained in the rows, and the
#'   columns contain the variance components for the subject level intercepts.}
#'   \item{\code{gamma_cov_bar}}{A matrix containing the group level regression
#'   coefficients of the multinomial logistic regression predicting the
#'   transition probabilities over the iterations of the hybrid Metropolis within
#'   Gibbs sampler. The iterations of the sampler are contained in the rows, and
#'   the columns contain the group level regression coefficients.}
#'   \item{\code{gamma_int_subj}}{A list containing one matrix per subject
#'   denoting the subject level intercepts of the multinomial logistic
#'   regression modeling the transition probabilities over the iterations of the
#'   hybrid Metropolis within Gibbs sampler. The iterations of the sampler are
#'   contained in the rows, and the columns contain the subject level
#'   intercepts.}
#'   \item{\code{gamma_naccept}}{A matrix containing the number of accepted
#'   draws at the subject level RW Metropolis step for each set of parameters of
#'   the transition probabilities. The subjects are contained in the rows, and
#'   the columns contain the sets of parameters.}
#'   \item{\code{emiss_mu_bar}}{A list containing one matrix per dependent
#'   variable, denoting the group level means of the Normal emission
#'   distribution of each dependent variable over the iterations of the Gibbs
#'   sampler. The iterations of the sampler are contained in the rows of the
#'   matrix, and the columns contain the group level emission means. If
#'   covariates were included in the analysis, the group level means represent
#'   the predicted mean given that the covariate is at the average value for
#'   continuous covariates, or given that the covariate equals zero for
#'   dichotomous covariates.}
#'   \item{\code{emiss_varmu_bar}}{A list containing one matrix per dependent
#'   variable, denoting the variance between the subject level means of the
#'   Normal emision distributions over the iterations of the Gibbs sampler. The
#'   iterations of the sampler are contained in the rows of the matrix, and the
#'   columns contain the group level variance in the mean.}
#'   \item{\code{emiss_var_bar}}{A list containing one matrix per dependent
#'   variable, denoting the (fixed over subjects) variance of the Normal
#'   emission distributions over the iterations of the Gibbs sampler. The
#'   iterations of the sampler are contained in the rows of the matrix, and the
#'   columns contain the group level emission variances.}
#'   \item{\code{emiss_cov_bar}}{A list containing one matrix per dependent
#'   variable, denoting the group level regression coefficients predicting the
#'   emission means within each of the dependent variables over the iterations
#'   of the Gibbs sampler. The iterations of the sampler are contained in the
#'   rows  of the matrix, and the columns contain the group level regression
#'   coefficients.}
#'   \item{\code{dwell_mu_bar}}{A matrix denoting the group level means of the
#'   log-Normal dwell distribution over the iterations of the Gibbs
#'   sampler. The iterations of the sampler are contained in the rows of the
#'   matrix, and the columns contain the group level dwell log-means. Note that
#'   the log-means are returned in the logarithmic scale.}
#'   \item{\code{dwell_varmu_bar}}{A matrix denoting the variance between the
#'   subject level log-means of the log-Normal dwell distribution over the
#'   iterations of the Gibbs sampler. The iterations of the sampler are
#'   contained in the rows of the matrix, and the columns contain the group
#'   level variance of the dwell mean.}
#'   \item{\code{dwell_var_bar}}{A matrix denoting the (fixed over subjects)
#'   variance of the log-Normal dwell distribution over the
#'   iterations of the Gibbs sampler. The iterations of the sampler are
#'   contained in the rows of the matrix, and the columns contain the group
#'   level variance of the subject level dwell distribution.}
#'   \item{\code{input}}{Overview of used input specifications: the number of
#'   states \code{m}, the number of used dependent variables \code{n_dep}, the
#'   number of iterations \code{J} and the specified burn in period
#'   \code{burn_in} of the hybrid Metropolis within Gibbs sampler, the number of
#'   subjects \code{n_subj}, the observation length for each subject
#'   \code{n_vary}, and the column names of the dependent variables
#'   \code{dep_labels}.}
#'   \item{\code{sample_path}}{A list containing one matrix per subject with the
#'   sampled hidden state sequence over the hybrid Metropolis within Gibbs
#'   sampler. The time points of the dataset are contained in the rows, and the
#'   sampled paths over the iterations are contained in the columns. Only
#'   returned if \code{return_path = TRUE}. }
#' }
#'
#' @seealso \code{\link{sim_medHMM}} for simulating multilevel hidden Markov
#'   data, and \code{\link{vit_medHMM}} for obtaining the most likely hidden
#'   state sequence for each subject using the Viterbi algorithm.
#'
#' @references
#' \insertRef{rabiner1989}{mHMMbayes}
#'
#' \insertRef{scott2002}{mHMMbayes}
#'
#' \insertRef{altman2007}{mHMMbayes}
#'
#' \insertRef{rossi2012}{mHMMbayes}
#'
#' \insertRef{zucchini2017}{mHMMbayes}
#'
#'
#'
#'
#' @examples
#' ###### Example on simulated data
#' # simulating multivariate continuous data
#'
#' # 3 states
#' n_t <- 500
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
#' dwell_distr <- dwell_start <- matrix(log(c(50,2,
#'                                            10,2,
#'                                            2,2)), nrow = m, ncol = 2, byrow = TRUE)
#'
#' # Simulate data
#' sim_data <- sim_medHMM(n_t, n, data_distr = 'continuous', m, n_dep = n_dep,
#'                        dwell_distr = dwell_distr, dwell_type = 'logNormal',
#'                        start_state = NULL, q_emiss = NULL, gamma = gamma, emiss_distr = emiss_distr, xx_vec = NULL, beta = NULL,
#'                        var_gamma = 0.1, var_emiss = c(0.1,0.1), var_dwell = 0.01, return_ind_par = TRUE)
#'
#' # Specify hyper-prior for the continuous emission distribution
#' emiss_hyp_pr <- list(
#'     emiss_mu0 = list(matrix(c(10,50,2), nrow = 1),
#'                      matrix(c(-5, -20, 5), nrow = 1)),
#'     emiss_K0  = list(1, 1),
#'     emiss_nu  = list(1, 1),
#'     emiss_V   = list(rep(100, m), rep(100, m)),
#'     emiss_a0  = list(rep(0.01, m), rep(0.01, m)),
#'     emiss_b0  = list(rep(0.01, m), rep(0.01, m))
#' )
#'
#' # Specify hyper-prior for the dwelling time (log-normal) distribution
#' dwell_hyp_pr <- list(
#'     d_mu0			= log(c(20,10,2)),
#'     s2_0			= log(rep(100, m)),
#'     alpha.sigma20	= rep(0.01, m),
#'     beta.sigma20	= rep(0.01, m),
#'     alpha.tau20		= rep(0.01, m),
#'     beta.tau20		= rep(0.01, m)
#' )
#'
#'
#' # Run the model on the simulated data:
#' out <- medHMM_cont(s_data = sim_data$obs,
#'                    gen = list(m = m, n_dep = n_dep),
#'                    start_val = c(list(gamma), emiss_distr, list(dwell_start)),
#'                    emiss_hyp_prior = emiss_hyp_pr,
#'                    dwell_hyp_prior = dwell_hyp_pr,
#'                    show_progress = TRUE,
#'                    mcmc = list(J = 200, burn_in = 100), return_path = TRUE, max_dwell = 120)
#'
#' # To do:
#' #   - Improve function structure and output for 2 states
#' #   - Make sure function works for covariates
#' #   - Improve stability of function implementing control checks
#' #   - Make model more memory efficient (rm unnecessary objects)
#'
#' @export

medHMM_cont <- function(s_data, gen, xx = NULL, start_val, emiss_hyp_prior, dwell_hyp_prior,
                        mcmc, return_path = FALSE, show_progress = TRUE,
                        gamma_hyp_prior = NULL, gamma_sampler = NULL, max_dwell = NULL
){

    # Initialize data -----------------------------------
    # dependent variable(s), sample size, dimensions gamma and conditional distribuiton
    n_dep			 <- gen$n_dep
    dep_labels <- colnames(s_data[,2:(n_dep+1)])
    id         <- unique(s_data[,1])
    n_subj     <- length(id)
    subj_data  <- rep(list(NULL), n_subj)
    if(sum(sapply(s_data, is.factor)) > 0 ){
        stop("Your data contains factorial variables, which cannot be used as input in the function mHMM. All variables have to be numerical.")
    }
    for(s in 1:n_subj){
        subj_data[[s]]$y <- as.matrix(s_data[s_data[,1] == id[s],][,-1], ncol = n_dep)
    }
    ypooled    <- n <- NULL
    n_vary     <- numeric(n_subj)
    m          <- gen$m

    for(s in 1:n_subj){
        ypooled   <- rbind(ypooled, subj_data[[s]]$y)
        n         <- dim(subj_data[[s]]$y)[1]

        switched 			<- rep(0, n)
        switched[1] 		<- 1
        for (t in 2:n) {
            if(any(subj_data[[s]]$y[t] != subj_data[[s]]$y[t-1])) {
                switched[t] <- 1}
        }
        switched2 		<- c(switched[-1],1)

        if (is.null(max_dwell)){
            max_dwell_n <- n
        } else if (max_dwell > n) {
            max_dwell_n <- n
        } else {
            max_dwell_n <- max_dwell
        }

        n_vary[s] <- n
        subj_data[[s]]	<- c(
            subj_data[[s]], n = n, list(gamma_converge = numeric(m), gamma_int_mle = matrix(NA_real_, m, (m - 2)),
                                        gamma_mhess = matrix(NA_real_, (m - 2) * m, (m - 2))),
            list(Mx = max_dwell_n),
            list(Mx2 = rep(max_dwell_n, m)),
            list(switch = switched),
            list(switch2 = switched2)
        )
    }
    n_total 		<- dim(ypooled)[1]

    # Intialize covariates
    n_dep1 <- 1 + n_dep
    nx <- numeric(n_dep1)
    if (is.null(xx)){
        xx <- rep(list(matrix(1, ncol = 1, nrow = n_subj)), n_dep1)
        nx[] <- 1
    } else {
        if(!is.list(xx) | length(xx) != n_dep1){
            stop("If xx is specified, xx should be a list, with the number of elements equal to the number of dependent variables + 1")
        }
        for(i in 1:n_dep1){
            if (is.null(xx[[i]])){
                xx[[i]] <- matrix(1, ncol = 1, nrow = n_subj)
                nx[i] <- 1
            } else {
                nx[i] <- ncol(xx[[i]])
                if (sum(xx[[i]][,1] != 1)){
                    stop("If xx is specified, the first column in each element of xx has to represent the intercept. That is, a column that only consists of the value 1")
                }
                if(nx[i] > 1){
                    for(j in 2:nx[i]){
                        if(is.factor(xx[[i]][,j])){
                            stop("Factors currently cannot be used as covariates, see help file for alternatives")
                        }
                        if((length(unique(xx[[i]][,j])) == 2) & (sum(xx[[i]][,j] != 0 & xx[[i]][,j] !=1) > 0)){
                            stop("Dichotomous covariates in xx need to be coded as 0 / 1 variables. That is, only conisting of the values 0 and 1")
                        }
                        if(length(unique(xx[[i]][,j])) > 2){
                            xx[[i]][,j] <- xx[[i]][,j] - mean(xx[[i]][,j])
                        }
                    }
                }
            }
        }
    }

    # Initialize mcmc argumetns
    J 				<- mcmc$J
    burn_in			<- mcmc$burn_in


    # Initalize priors and hyper priors --------------------------------
    # Initialize gamma sampler
    if(is.null(gamma_sampler)) {
        gamma_int_mle0  <- rep(0, m - 2)
        gamma_scalar    <- 2.93 / sqrt(m - 2)
        gamma_w         <- .1
    } else {
        gamma_int_mle0  <- gamma_sampler$gamma_int_mle0
        gamma_scalar    <- gamma_sampler$gamma_scalar
        gamma_w         <- gamma_sampler$gamma_w
    }



    # Initialize Gamma hyper prior
    if(is.null(gamma_hyp_prior)){
        gamma_mu0	        <- rep(list(matrix(0,nrow = nx[1], ncol = m - 2)), m)
        gamma_K0			<- diag(1, nx[1])
        gamma_nu			<- 3 + m - 2
        gamma_V			    <- gamma_nu * diag(m - 2)
    } else {
        ###### BUILD in a warning / check if gamma_mu0 is a matrix when given, with  nrows equal to the number of covariates
        gamma_mu0			<- gamma_hyp_prior$gamma_mu0
        gamma_K0			<- diag(gamma_hyp_prior$gamma_K0, nx[1])
        gamma_nu			<- gamma_hyp_prior$gamma_nu
        gamma_V			    <- gamma_hyp_prior$gamma_V
    }


    # Initialize emiss hyper prior
    if(missing(emiss_hyp_prior)){
        stop("The hyper-prior values for the Normal emission distribution(s) denoted by emiss_hyp_prior needs to be specified")
    }

    # emiss_mu0: a list containing n_dep matrices with in the first row the hypothesized mean values of the Normal emission
    # distributions in each of the states over the m coloumns. Subsequent rows contain the hypothesised regression
    # coefficients for covariates influencing the state dependent mean value of the normal distribution
    emiss_mu0	  <- rep(list(NULL), n_dep)
    emiss_a0	  <- rep(list(NULL), n_dep)
    emiss_b0	  <- rep(list(NULL), n_dep)
    emiss_V	  <- rep(list(NULL), n_dep)
    emiss_nu	    <- rep(list(NULL), n_dep)
    emiss_K0     <- rep(list(NULL), n_dep)
    for(q in 1:n_dep){
        # emiss_hyp_prior[[q]]$emiss_mu0 has to contain a list with lenght equal to m, and each list contains matrix with number of rows equal to number of covariates for that dep. var.
        # stil build in a CHECK, with warning / stop / switch to default prior
        emiss_mu0[[q]]	 <- emiss_hyp_prior$emiss_mu0[[q]]
        emiss_nu[[q]]	 <- emiss_hyp_prior$emiss_nu[[q]]
        emiss_V[[q]]		 <- emiss_hyp_prior$emiss_V[[q]]
        emiss_K0[[q]]	 <- diag(emiss_hyp_prior$emiss_K0, nx[1 + q])
        emiss_a0[[q]] <- emiss_hyp_prior$emiss_a0[[q]]
        emiss_b0[[q]] <- emiss_hyp_prior$emiss_b0[[q]]
    }


    # Initialize dwell hyper prior =====================================================================================================
    d_mu0			<- dwell_hyp_prior$d_mu0
    s2_0			<- dwell_hyp_prior$s2_0
    alpha.sigma20	<- dwell_hyp_prior$alpha.sigma20
    beta.sigma20	<- dwell_hyp_prior$beta.sigma20
    alpha.tau20		<- dwell_hyp_prior$alpha.tau20
    beta.tau20		<- dwell_hyp_prior$beta.tau20



    # Define objects used to store data in mcmc algorithm, not returned ----------------------------
    # overall
    c <- llk <- numeric(1)
    sample_path <- lapply(n_vary, dif_matrix, cols = J)
    trans <- rep(list(vector("list", m)), n_subj)

    # gamma
    gamma_int_mle_pooled <- gamma_pooled_ll <- vector("list", m)
    gamma_c_int <- rep(list(matrix(NA_real_, n_subj, (m-2))), m)
    gamma_mu_int_bar <- gamma_V_int <- vector("list", m)
    gamma_mu_prob_bar <- rep(list(numeric(m)), m)
    gamma_naccept <- matrix(0, n_subj, m)

    # emiss
    cond_y <- lapply(rep(n_dep, n_subj), nested_list, m = m)
    cond_y_pooled <- rep(list(rep(list(NULL),n_dep)), m)
    emiss_c_mu <- rep(list(rep(list(matrix(NA_real_,ncol = 1, nrow = n_subj)),n_dep)), m)
    for(i in 1:m){
        for(q in 1:n_dep){
            emiss_c_mu[[i]][[q]][,1] <- start_val[[1 + q]][i,1]
        }
    }
    emiss_V_mu <- emiss_c_mu_bar <- emiss_c_V <- rep(list(rep(list(NULL),n_dep)), m)
    ss_subj <- n_cond_y <- numeric(n_subj)
    label_switch <- matrix(0, ncol = n_dep, nrow = m, dimnames = list(c(paste("mu_S", 1:m, sep = "")), dep_labels))

    # dwell ============================================================================================================================= Added dwell pars
    d <- B_star <- Occupancy <- N <- numeric()
    Dur <- vector("list", n_subj)
    # sample.path <- Dur <- vector("list", n.mice)
    n.Dur <- numeric(n_subj)
    tau2 <- logmu2 <- a1 <- b1 <- double(m)


    # Define objects that are returned from mcmc algorithm ----------------------------
    # Define object for subject specific posterior density, put start values on first row
    if(length(start_val) != n_dep + 2){
        stop("The number of elements in the list start_val should be equal to 2 + the number of dependent variables,
         and should not contain nested lists (i.e., lists within lists)")
    }
    PD 					  <- matrix(NA_real_, nrow = J, ncol = n_dep * m * 2 + m * m + m + m + 1) # ===================================================== Addedd dwell time pars and positions
    colnames(PD) 	<- c(paste("dep", rep(1:n_dep, each = m), "_mu", "_S", rep(1:m), sep = ""),
                       paste("dep", rep(1:n_dep, each = m), "_fixvar", "_S", rep(1:m), sep = ""),
                       paste("S", rep(1:m, each = m), "toS", rep(1:m, m), sep = ""),
                       paste("dur_logmu", 1:m, sep = ""),
                       paste("dur_logsigma2", 1:m, sep = ""),
                       # paste("dur_tau2_bar", 1:m, sep = ""),
                       "LL")

    PD[1, ((n_dep * m * 2 + 1)) :((n_dep * m * 2 + m * m))] <- unlist(sapply(start_val, t))[1:(m*m)]
    for(q in 1:n_dep){
        PD[1, ((q-1) * m + 1):(q * m)] <- start_val[[q + 1]][,1]
        PD[1, (n_dep * m + (q-1) * m + 1):(n_dep * m + q * m)] <- start_val[[q + 1]][,2]
    }
    PD[1, ((n_dep * m * 2 + m * m + 1)) :((n_dep * m * 2 + m * m + m))] <- start_val[[1 + n_dep + 1]][,1]
    PD[1, ((n_dep * m * 2 + m * m + m + 1)) :((n_dep * m * 2 + m * m + m + m))] <- start_val[[1 + n_dep + 1]][,2]
    # PD[1, ((n_dep * m * 2 + m * m + 1)) :((n_dep * m * 2 + m * m + m))] <- unlist(sapply(start_val, t))[(n_dep*m*2+m*m+1):(n_dep*m*2+m*m+m*2)]

    PD_subj				<- rep(list(PD), n_subj)

    # Define object for population posterior density (probabilities and regression coefficients parameterization )
    gamma_prob_bar		<- matrix(NA_real_, nrow = J, ncol = (m * m))
    colnames(gamma_prob_bar) <- paste("S", rep(1:m, each = m), "toS", rep(1:m, m), sep = "")
    gamma_prob_bar[1,] <- PD[1,(n_dep * 2 * m + 1):(n_dep * 2 * m + m * m)]
    emiss_mu_bar			<- rep(list(matrix(NA_real_, ncol = m, nrow = J, dimnames = list(NULL, c(paste("mu_", 1:m, sep = ""))))), n_dep)
    names(emiss_mu_bar) <- dep_labels
    for(q in 1:n_dep){
        emiss_mu_bar[[q]][1,] <- PD[1, ((q-1) * m + 1):(q * m)]
    }

    gamma_int_bar				<- matrix(NA_real_, nrow = J, ncol = ((m-2) * m))
    temp_gamma <- matrix(gamma_prob_bar[1,], byrow = TRUE, ncol = m, nrow = m)
    diag(temp_gamma) <- NA
    gamma_int_bar[1,] <- as.vector(prob_to_int(t(apply(matrix(as.numeric(t(temp_gamma))[!is.na(as.numeric(t(temp_gamma)))], nrow = m, byrow = TRUE),1, function(r) r/sum(r)))))

    temp <- matrix(paste0("int_S", rep(1:m, each = m), "toS", rep(1:m, m)), nrow = m, byrow = TRUE)
    diag(temp) <- NA
    colnames(gamma_int_bar) <- as.vector(t(matrix(t(temp)[!is.na(as.vector(t(temp)))], nrow = m, byrow = TRUE)[,-1]))

    gamma_V_int_bar <- matrix(NA_real_, nrow = J, ncol = ((m-2) * (m-2) * m))
    gamma_V_int_bar[1,] <- unlist(lapply(gamma_V, function(e) as.vector(t(e))))

    int_names <- paste0("int_S", rep(1:m, each = m), "toS", rep(1:m, m))[!(paste0("int_S", rep(1:m, each = m), "toS", rep(1:m, m)) %in% as.vector(t(matrix(t(temp)[!is.na(as.vector(t(temp)))], nrow = m, byrow = TRUE)[,-1])))]
    var_names <- paste0("var_int_S", rep(1:m, each = (m-1)*(m-1)), "toS", rep(2:m, each=m-1), "_with_", "int_S", rep(1:m, each = (m-1)*(m-1)), "toS", rep(2:m, m))
    colnames(gamma_V_int_bar) <- var_names[!grepl(paste(int_names, collapse = "|"), var_names)]

    if(nx[1] > 1){
        gamma_cov_bar				<- matrix(NA_real_, nrow = J, ncol = ((m-2) * m) * (nx[1] - 1))
        # colnames(gamma_cov_bar) <- paste( paste("cov", rep(1 : (nx[1] - 1),each = (m-1)), "_", sep = ""), "S", rep(1:m, each = (m-1) * (nx[1] - 1)), "toS", rep(2:m, m * (nx[1] - 1)), sep = "") # CHECK NAMES
        gamma_cov_bar[1,] <- 0
    } else{
        gamma_cov_bar <- "No covariates where used to predict the transition probability matrix"
    }
    emiss_varmu_bar			<- rep(list(matrix(NA_real_, ncol = m, nrow = J, dimnames = list(NULL, c(paste("varmu_", 1:m, sep = ""))))), n_dep)
    names(emiss_varmu_bar) <- dep_labels
    emiss_var_bar			<- rep(list(matrix(NA_real_, ncol = m, nrow = J, dimnames = list(NULL, c(paste("var_", 1:m, sep = ""))))), n_dep)
    names(emiss_var_bar) <- dep_labels
    for(q in 1:n_dep){
        emiss_var_bar[[q]][1,] <- PD[1, (n_dep * m + (q-1) * m + 1):(n_dep * m + q * m)]
    }
    if(sum(nx[-1]) > n_dep){
        emiss_cov_bar			<- lapply(m * (nx[-1] - 1 ), dif_matrix, rows = J)
        names(emiss_cov_bar) <- dep_labels
        for(q in 1:n_dep){
            if(nx[1 + q] > 1){
                colnames(emiss_cov_bar[[q]]) <-  paste( paste("cov", 1 : (nx[1+q] - 1), "_", sep = ""), "mu_S", rep(1:m, each = (nx[1 + q] - 1)), sep = "")
                emiss_cov_bar[[q]][1,] <- 0
            } else {
                emiss_cov_bar[[q]] <- "No covariates where used to predict the emission probabilities for this outcome"
            }
        }
    } else{
        emiss_cov_bar <- "No covariates where used to predict the emission probabilities"
    }

    # dwell ============================================================================================================================= Added dwell pars
    # duration_bar		<- matrix(NA_real_, nrow = J, ncol = m * 3)
    # colnames(duration_bar)	<- c(paste("mu_d_bar", 1:m, sep = ""), paste("logsigma2", 1:m, sep = ""), paste("tau2_d_bar", 1:m, sep = ""))
    dwell_mu_bar <- matrix(NA_real_, nrow = J, ncol = m)
    colnames(dwell_mu_bar) <- paste("mu_d_bar", 1:m, sep = "")
    dwell_mu_bar[1,] <- matrix(PD[1, ((n_dep * m * 2 + m * m + 1)) :((n_dep * m * 2 + m * m + m))], nrow = 1, ncol = m, byrow = TRUE)

    dwell_varmu_bar <- dwell_var_bar <- matrix(NA_real_, nrow = J, ncol = m)
    colnames(dwell_varmu_bar) <- paste("tau2_d_bar", 1:m, sep = "")
    colnames(dwell_var_bar) <- paste("logsigma2", 1:m, sep = "")

    dwell_var_bar[1,] <- dwell_varmu_bar[1,] <- matrix(PD[1, ((n_dep * m * 2 + m * m + 1)) :((n_dep * m * 2 + m * m + m))], nrow = 1, ncol = m, byrow = TRUE)


    # Define object for subject specific posterior density (regression coefficients parameterization )
    gamma_int_subj			<- rep(list(gamma_int_bar), n_subj)

    # Put starting values in place for fist run forward algorithm
    emiss_sep 	<- rep(list(matrix(NA_real_, ncol = 2, nrow = m)), n_dep)
    for(q in 1:n_dep){
        emiss_sep[[q]][,1] <- PD[1, ((q-1) * m + 1):(q * m)]
        emiss_sep[[q]][,2] <- PD[1, (n_dep * m + (q-1) * m + 1):(n_dep * m + q * m)]
    }
    emiss				<- rep(list(emiss_sep), n_subj)
    gamma 			<- rep(list(matrix(PD[1,(m * 2 * n_dep + 1):(m * 2 * n_dep + m * m)], byrow = TRUE, ncol = m)), n_subj)
    delta 			<- rep(list(solve(t(diag(m) - gamma[[1]] + 1), rep(1, m))), n_subj)

    # dwell ============================================================================================================================= Added dwell starting values
    mu_d_bar        <- matrix(PD[1, ((n_dep * m * 2 + m * m + 1)) :((n_dep * m * 2 + m * m + m))], nrow = 1, ncol = m, byrow = TRUE)
    logmu           <- matrix(PD[1, ((n_dep * m * 2 + m * m + 1)) :((n_dep * m * 2 + m * m + m))], nrow = n_subj, ncol = m, byrow = TRUE)
    logsigma2		<- tau2_d_bar <- PD[1, ((n_dep * m * 2 + m * m + m + 1)) :((n_dep * m * 2 + m * m + m + m))]
    # logmu		  	<- mu_d_bar <- matrix(Start_Val[((m*q_pr + m*m) + 1) : ((m*q_pr + m*m) + m)], n.mice, m, byrow = TRUE)
    # logsigma2		<- tau2_d_bar <- Start_Val[((m*q_pr + m*m) + m + 1) : ((m*q_pr + m*m) + m*2)]


    # Start analysis --------------------------------------------
    # Run the MCMC algorithm

    if (m==2){cat("A model with only two hidden states has been specified. As a result, all transition probabilities will be fixed to 1, and no variation between individuals will be modelled.\n")}

    itime <- proc.time()[3]
    if(show_progress == TRUE){
        cat("Progress of the Bayesian mHMM algorithm:", "\n")
        pb <- utils::txtProgressBar(min = 2, max = J, style = 3)
    }

    for (iter in 2 : J){

        if (m == 2){

            # For each subject, obtain sampled state sequence with subject individual parameters ----------
            sample_path_state <- Dur <- vector("list", n_subj)

            for(s in 1:n_subj){

                # Idea: pre-compute emissions likelihood to pass to FBalgC() here:

                # Run forward backward algorithm in C++, using the runlength distribution d for each state ================================================================
                d 	<- get.d.lognorm(run.p = list(logmu = logmu[s,], logsd = sqrt(logsigma2)), Mx = subj_data[[s]]$Mx, m = m)

                delta[[s]] <- get_delta(gamma[[s]], m)

                allprobs <- get_all1(x = subj_data[[s]]$y, emiss = emiss[[s]], n_dep = n_dep, data_distr = "continuous")

                FB	<- mult_ed_fb_cpp(
                    # y2 = subj_data[[s]]$y,
                    m = m,
                    n = subj_data[[s]]$n,
                    allprobs = t(allprobs),
                    Mx = subj_data[[s]]$Mx,
                    Mx2 = subj_data[[s]]$Mx2,
                    gamma = gamma[[s]],
                    d = d,
                    S2 = subj_data[[s]]$switch2,
                    S = subj_data[[s]]$switch,
                    delta = delta[[s]]
                )

                B_star				<- FB[[2]]
                Occupancy			<- FB[[3]]
                N 					<- FB[[1]]
                PD_subj[[s]][iter-1, m*n_dep*2 + m*m + m*2 + 1] <- llk <- sum(N)	# adjust index; we may need to do log-sum-ex

                # Using the outcomes of the forward backward algorithm, sample the state sequence ==========================================================================
                trans[[s]]				                <- vector("list", m)
                sample_path_state[[s]][1] 	            <- sample(1:m, 1, prob = delta[[s]] * exp(B_star[,1]))
                Dur[[s]][1] 			                <- sample(1:subj_data[[s]]$Mx, 1, prob = (Occupancy[[1]][sample_path_state[[s]][1],] / exp(B_star[sample_path_state[[s]][1],1])))
                sample_path[[s]][1:Dur[[s]][1], iter]   <- sample_path_state[[s]][1]

                t <- 1
                while(sum(Dur[[s]]) < subj_data[[s]]$n){
                    t 						                                        <- t + 1
                    Mx.l 					                                        <- min(subj_data[[s]]$Mx, subj_data[[s]]$n-sum(Dur[[s]]))
                    sample_path_state[[s]][t] 	                                    <- sample(1:m, 1, prob = gamma[[s]][sample_path_state[[s]][t-1],] * exp(B_star[,sum(Dur[[s]])+1]))
                    trans[[s]][[sample_path_state[[s]][t-1]]]                       <- c(trans[[s]][[sample_path_state[[s]][t-1]]], sample_path_state[[s]][t])
                    Dur[[s]][t]			                                            <- sample(1:Mx.l, 1, prob = (Occupancy[[sum(Dur[[s]])+1]][sample_path_state[[s]][t],] / exp(B_star[sample_path_state[[s]][t], sum(Dur[[s]])+1])))
                    sample_path[[s]][sum(Dur[[s]][1:t-1],1):sum(Dur[[s]]), iter]    <- sample_path_state[[s]][t]
                }

                n.Dur[s] <- length(Dur[[s]])
                for (i in 1:m){
                    # trans[[s]][[i]]         <- c(trans[[s]][[i]], 1:m) # to avoid errors, check if we can drop
                    for (q in 1:n_dep) {
                        # cond_y[[s]][[i]]    <- subj_data[[s]]$y[sample_path[[s]][, iter] == i, q]
                        if(iter == 2){
                            cond_y[[s]][[i]][[q]] <- c(subj_data[[s]]$y[sample_path[[s]][, iter] == i, q][!is.na(subj_data[[s]]$y[sample_path[[s]][, iter] == i, q])],emiss_mu0[[q]][1,i])
                        } else {
                            cond_y[[s]][[i]][[q]] <- c(subj_data[[s]]$y[sample_path[[s]][, iter] == i, q][!is.na(subj_data[[s]]$y[sample_path[[s]][, iter] == i, q])],emiss_c_mu_bar[[i]][[q]][1])
                        }
                    }

                }
            }

            # The remainder of the mcmc algorithm is state specific
            for(i in 1:m){

                # Sample populaton values for gamma and conditional probabilities using Gibbs sampler -----------
                # gamma_mu0_n and gamma_mu_int_bar are matrices, with the number of rows equal to the number of covariates, and ncol equal to number of intercepts estimated
                gamma_V_int[[i]]      <- 0
                gamma_mu_int_bar[[i]] <- matrix(matrix(c(Inf,
                                                             -Inf), nrow = 2, byrow = TRUE)[i,], nrow = 1)
                gamma_mu_prob_bar[[i]] 	<- as.vector(matrix(c(0,1,
                                                              1,0), nrow = 2, byrow = TRUE)[i,])

                # sample population mean (and regression parameters if covariates) of the Normal emission distribution, and it's variance (so the variance between the subject specific means)
                # note: the posterior is thus one of a Bayesian linear regression because of the optional regression parameters
                for(q in 1:n_dep){
                    emiss_mu0_n                    <- solve(t(xx[[1 + q]]) %*% xx[[1 + q]] + emiss_K0[[q]]) %*% (t(xx[[1 + q]]) %*% emiss_c_mu[[i]][[q]] + emiss_K0[[q]] %*% emiss_mu0[[q]][,i])
                    emiss_a_mu_n                   <- (emiss_K0[[q]] + n_subj) / 2
                    emiss_b_mu_n                   <- (emiss_nu[[q]] * emiss_V[[q]][i]) / 2 + (t(emiss_c_mu[[i]][[q]]) %*% emiss_c_mu[[i]][[q]] +
                                                                                                   t(emiss_mu0[[q]][,i]) %*% emiss_K0[[q]] %*% emiss_mu0[[q]][,i] -
                                                                                                   t(emiss_mu0_n) %*% (t(xx[[1 + q]]) %*% xx[[1 + q]] + emiss_K0[[q]]) %*% emiss_mu0_n) / 2
                    emiss_V_mu[[i]][[q]]       <- solve(stats::rgamma(1, shape = emiss_a_mu_n, rate = emiss_b_mu_n))
                    if(all(dim(emiss_V_mu[[i]][[q]]) == c(1,1))){
                        emiss_c_mu_bar[[i]][[q]]	  <- emiss_mu0_n + rnorm(1 + nx[1 + q] - 1, mean = 0, sd = sqrt(diag(as.numeric(emiss_V_mu[[i]][[q]]) * solve(t(xx[[1 + q]]) %*% xx[[1 + q]] + emiss_K0[[q]]))))
                    } else {
                        emiss_c_mu_bar[[i]][[q]]	  <- emiss_mu0_n + rnorm(1 + nx[1 + q] - 1, mean = 0, sd = sqrt(diag(emiss_V_mu[[i]][[q]] * solve(t(xx[[1 + q]]) %*% xx[[1 + q]] + emiss_K0[[q]]))))
                    }
                }

                # Sample subject values  -----------
                for (s in 1:n_subj){

                    gamma[[s]][i,]  	    <- PD_subj[[s]][iter, c((n_dep * 2 * m + 1 + (i - 1) * m):(n_dep * 2 * m + (i - 1) * m + m))] <- as.vector(matrix(c(0,1,
                                                                                                                                                             1,0), nrow = 2, byrow = TRUE)[i,])
                    gamma_naccept[s, i]		<- gamma_naccept[s, i] + 1
                    gamma_c_int[[i]][s,]	<- matrix(c(Inf,
                                                                -Inf), nrow = 2, byrow = TRUE)[i,]
                    # gamma_int_subj[[s]][iter, c((1 + (i - 1) * (m - 2)):((m - 2) + (i - 1) * (m - 2)))[-i]] <- gamma_c_int[[i]][s,] # CHECK
                    gamma_int_subj[[s]][iter, c((1 + (i - 1) * (m - 2)):((m - 2) + (i - 1) * (m - 2)))] <- gamma_c_int[[i]][s,]

                    if(i == m){
                        delta[[s]] 		<- solve(t(diag(m) - gamma[[s]] + 1), rep(1, m))
                    }
                }
                # Sample subject values for normal emission distribution using Gibbs sampler   ---------

                # population level, conditional probabilities, seperate for each dependent variable
                for(q in 1:n_dep){
                    for (s in 1:n_subj){
                        ss_subj[s] <- t(matrix(cond_y[[s]][[i]][[q]] - emiss_c_mu[[i]][[q]][s,1], nrow = 1) %*%
                                            matrix(cond_y[[s]][[i]][[q]] - emiss_c_mu[[i]][[q]][s,1], ncol = 1))
                        n_cond_y[s]       <- length(cond_y[[s]][[i]][[q]])
                    }
                    emiss_a_resvar_n <- sum(n_cond_y) / 2 + emiss_a0[[q]][i]
                    emiss_b_resvar_n <- (sum(ss_subj) + 2 * emiss_b0[[q]][i]) / 2
                    emiss_c_V[[i]][[q]] <- emiss_var_bar[[q]][iter, i] <- solve(stats::rgamma(1, shape = emiss_a_resvar_n, rate = emiss_b_resvar_n))
                }

                ### sampling subject specific means for the emission distributions, assuming known mean and var, see Lynch p. 244
                for(q in 1:n_dep){
                    emiss_c_V_subj    <- (emiss_V_mu[[i]][[q]] * emiss_c_V[[i]][[q]]) / (2 * emiss_V_mu[[i]][[q]] + emiss_c_V[[i]][[q]])
                    for (s in 1:n_subj){
                        emiss_mu0_subj_n  <- (emiss_V_mu[[i]][[q]] * sum(cond_y[[s]][[i]][[q]]) +  emiss_c_V[[i]][[q]] * c(t(emiss_c_mu_bar[[i]][[q]]) %*% xx[[q+1]][s,])) /
                            (n_cond_y[s] * emiss_V_mu[[i]][[q]] + emiss_c_V[[i]][[q]])
                        emiss[[s]][[q]][i,1] <- PD_subj[[s]][iter, ((q - 1) * m + i)] <- emiss_c_mu[[i]][[q]][s,1] <- rnorm(1, emiss_mu0_subj_n, sqrt(emiss_c_V_subj))
                        emiss[[s]][[q]][i,2] <- PD_subj[[s]][iter, (n_dep * m + (q - 1) * m + i)] <- emiss_c_V[[i]][[q]]
                    }
                }


                #################
                # Obtain hierarchical and mouse specific parameters for duration distribuiton using gibbs sampler ====================================
                #################

                #draw logmu's
                for(s in 1:n_subj){
                    tau2 		<- 1/ ((1/tau2_d_bar[i]) + (1/logsigma2[i]) * sum(sample_path_state[[s]][-n.Dur[s]] == i))
                    logmu[s,i]	<- rnorm(1,
                                        mean = tau2 * ((1/tau2_d_bar[i]) * mu_d_bar[i] + (1/logsigma2[i]) * sum(log(Dur[[s]][-n.Dur[s]][sample_path_state[[s]][-n.Dur[s]] == i]))),
                                        sd = sqrt(tau2))
                    PD_subj[[s]][iter, n_dep*m*2 + m*m + i]			<- logmu[s, i] # adjust index
                }

                # draw mu_d_bar
                s2				<- 1 / ((1/s2_0[i]) + (1/tau2_d_bar[i]) * n_subj)
                mu_d_bar[i] 	<- rnorm(1,
                                      mean = s2 * ((1/s2_0[i]) * d_mu0[i] + (1/tau2_d_bar[i]) * sum(logmu[,i])),
                                      sd = sqrt(s2))

                # draw tau2_d_bar
                a1 <- alpha.tau20[i] + n_subj/2
                b1 <- beta.tau20[i] + (1/2) * (sum((logmu[,i] - mu_d_bar[i])^2))
                tau2_d_bar[i] <- 1 / rgamma(1, shape = a1, rate = b1)

                #draw logsigma
                ss <- numeric(1)
                n.ss <- numeric(1)
                for (s in 1:n_subj){
                    ss <- ss + sum((log(Dur[[s]][-n.Dur[s]][sample_path_state[[s]][-n.Dur[s]] == i]) - logmu[s,i])^2)
                    n.ss <- n.ss + sum(sample_path_state[[s]][-n.Dur[s]] == i)
                }
                c1 <- alpha.sigma20[i] + n.ss/2
                d1 <- beta.sigma20[i] + (1/2) * ss
                logsigma2[i] <- 1 / rgamma(1, shape = c1, rate = d1)
                for (s in 1:n_subj) {
                    PD_subj[[s]][iter, n_dep*m*2 + m*m + m + i] <- logsigma2[i]
                }

            }

        } else if (m >= 3){

            # For each subject, obtain sampled state sequence with subject individual parameters ----------
            sample_path_state <- Dur <- vector("list", n_subj)

            for(s in 1:n_subj){

                # Idea: pre-compute emissions likelihood to pass to FBalgC() here:

                # Run forward backward algorithm in C++, using the runlength distribution d for each state ================================================================
                d 	<- get.d.lognorm(run.p = list(logmu = logmu[s,], logsd = sqrt(logsigma2)), Mx = subj_data[[s]]$Mx, m = m)

                delta[[s]] <- get_delta(gamma[[s]], m)

                allprobs <- get_all1(x = subj_data[[s]]$y, emiss = emiss[[s]], n_dep = n_dep, data_distr = "continuous")

                FB	<- mult_ed_fb_cpp(
                    # y2 = subj_data[[s]]$y,
                    m = m,
                    n = subj_data[[s]]$n,
                    allprobs = t(allprobs),
                    Mx = subj_data[[s]]$Mx,
                    Mx2 = subj_data[[s]]$Mx2,
                    gamma = gamma[[s]],
                    d = d,
                    S2 = subj_data[[s]]$switch2,
                    S = subj_data[[s]]$switch,
                    delta = delta[[s]]
                )

                B_star				<- FB[[2]]
                Occupancy			<- FB[[3]]
                N 					<- FB[[1]]
                PD_subj[[s]][iter-1, m*n_dep*2 + m*m + m*2 + 1] <- llk <- sum(N)	# adjust index; we may need to do log-sum-ex

                # Using the outcomes of the forward backward algorithm, sample the state sequence ==========================================================================
                trans[[s]]				                <- vector("list", m)
                sample_path_state[[s]][1] 	            <- sample(1:m, 1, prob = delta[[s]] * exp(B_star[,1]))
                Dur[[s]][1] 			                <- sample(1:subj_data[[s]]$Mx, 1, prob = (Occupancy[[1]][sample_path_state[[s]][1],] / exp(B_star[sample_path_state[[s]][1],1])))
                sample_path[[s]][1:Dur[[s]][1], iter]   <- sample_path_state[[s]][1]

                t <- 1
                while(sum(Dur[[s]]) < subj_data[[s]]$n){
                    t 						                                        <- t + 1
                    Mx.l 					                                        <- min(subj_data[[s]]$Mx, subj_data[[s]]$n-sum(Dur[[s]]))
                    sample_path_state[[s]][t] 	                                    <- sample(1:m, 1, prob = gamma[[s]][sample_path_state[[s]][t-1],] * exp(B_star[,sum(Dur[[s]])+1]))
                    trans[[s]][[sample_path_state[[s]][t-1]]]                       <- c(trans[[s]][[sample_path_state[[s]][t-1]]], sample_path_state[[s]][t])
                    Dur[[s]][t]			                                            <- sample(1:Mx.l, 1, prob = (Occupancy[[sum(Dur[[s]])+1]][sample_path_state[[s]][t],] / exp(B_star[sample_path_state[[s]][t], sum(Dur[[s]])+1])))
                    sample_path[[s]][sum(Dur[[s]][1:t-1],1):sum(Dur[[s]]), iter]    <- sample_path_state[[s]][t]
                }

                n.Dur[s] <- length(Dur[[s]])
                for (i in 1:m){
                    # trans[[s]][[i]]         <- c(trans[[s]][[i]], 1:m) # to avoid errors, check if we can drop
                    for (q in 1:n_dep) {
                        # cond_y[[s]][[i]]    <- subj_data[[s]]$y[sample_path[[s]][, iter] == i, q]
                        if(iter == 2){
                            cond_y[[s]][[i]][[q]] <- c(subj_data[[s]]$y[sample_path[[s]][, iter] == i, q][!is.na(subj_data[[s]]$y[sample_path[[s]][, iter] == i, q])],emiss_mu0[[q]][1,i])
                        } else {
                            cond_y[[s]][[i]][[q]] <- c(subj_data[[s]]$y[sample_path[[s]][, iter] == i, q][!is.na(subj_data[[s]]$y[sample_path[[s]][, iter] == i, q])],emiss_c_mu_bar[[i]][[q]][1])
                        }
                    }

                }
            }

            # The remainder of the mcmc algorithm is state specific
            for(i in 1:m){

                # Obtain MLE of the covariance matrices and log likelihood of gamma and emiss at subject and population level -----------------
                # used to scale the propasal distribution of the RW Metropolis sampler

                # population level, transition matrix ============================================================ Check number of categories in gamma_int_mle0
                trans_pooled			<- factor(c(unlist(sapply(trans, "[[", i)), c(1:m)[-i]))
                trans_pooled            <- factor(trans_pooled, levels = paste(c(1:m)[-i]))
                trans_pooled            <- as.numeric(trans_pooled)
                gamma_mle_pooled		<- optim(gamma_int_mle0, llmnl_int, Obs = trans_pooled,
                                           method = "BFGS", hessian = FALSE,
                                           control = list(fnscale = -1))

                gamma_int_mle_pooled[[i]]   <- gamma_mle_pooled$par
                gamma_pooled_ll[[i]]        <- gamma_mle_pooled$value

                # subject level
                for (s in 1:n_subj){
                    wgt 		        <- subj_data[[s]]$n / n_total

                    trans_subj          <- factor(c(trans[[s]][[i]], c(1:m)[-i]))
                    trans_subj          <- factor(trans_subj, levels = paste(c(1:m)[-i]))
                    trans_subj          <- as.numeric(trans_subj)

                    # subject level, transition matrix ============================================================ Check number of categories in gamma_int_mle_pooled
                    gamma_out			<- optim(gamma_int_mle_pooled[[i]], llmnl_int_frac, Obs = trans_subj,
                                         # n_cat = m,
                                         pooled_likel = gamma_pooled_ll[[i]], w = gamma_w, wgt = wgt,
                                         method="BFGS", hessian = FALSE, control = list(fnscale = -1))
                    if(gamma_out$convergence == 0){
                        subj_data[[s]]$gamma_converge[i] <- 1
                        subj_data[[s]]$gamma_int_mle[i,] <- gamma_out$par
                        subj_data[[s]]$gamma_mhess[(1 + (i - 1) * (m - 2)):((m - 2) + (i - 1) * (m - 2)), ]	<-
                            mnlHess_int(int = gamma_out$par, Obs = c(trans[[s]][[i]], c(1:m)), n_cat =  (m-1) )
                    } else {
                        subj_data[[s]]$gamma_converge[i] <- 0
                        subj_data[[s]]$gamma_int_mle[i,] <- rep(0, m - 2)
                        subj_data[[s]]$gamma_mhess[(1 + (i - 1) * (m - 2)):((m - 2) + (i - 1) * (m - 2)), ]	<- diag(m-2)
                    }

                    # if this is first iteration, use MLE for current values RW metropolis sampler
                    if (iter == 2){
                        gamma_c_int[[i]][s,]		<- gamma_out$par
                    }
                }

                # Sample pouplaton values for gamma and conditional probabilities using Gibbs sampler -----------
                # gamma_mu0_n and gamma_mu_int_bar are matrices, with the number of rows equal to the number of covariates, and ncol equal to number of intercepts estimated
                gamma_mu0_n           <- solve(t(xx[[1]]) %*% xx[[1]] + gamma_K0)  %*% (t(xx[[1]]) %*% gamma_c_int[[i]] + gamma_K0 %*% gamma_mu0[[i]])
                gamma_V_n             <- gamma_V + t(gamma_c_int[[i]] - xx[[1]] %*% gamma_mu0_n) %*% (gamma_c_int[[i]] - xx[[1]] %*% gamma_mu0_n) + t(gamma_mu0_n - gamma_mu0[[i]]) %*% gamma_K0 %*% (gamma_mu0_n - gamma_mu0[[i]])
                gamma_V_int[[i]]      <- solve(rwish(S = solve(gamma_V_n), v = gamma_nu + n_subj))
                gamma_mu_int_bar[[i]] <- gamma_mu0_n + solve(chol(t(xx[[1]]) %*% xx[[1]] + gamma_K0)) %*% matrix(rnorm((m - 2) * nx[1]), nrow = nx[1]) %*% t(solve(chol(solve(gamma_V_int[[i]]))))
                gamma_exp_int				  <- matrix(exp(c(0, gamma_mu_int_bar[[i]][1,] )), nrow  = 1)
                gamma_mu_prob_bar[[i]] 	<- gamma_exp_int / as.vector(gamma_exp_int %*% c(rep(1,(m-1))))

                # sample population mean (and regression parameters if covariates) of the Normal emission distribution, and it's variance (so the variance between the subject specific means)
                # note: the posterior is thus one of a Bayesian linear regression because of the optional regression parameters
                for(q in 1:n_dep){
                    emiss_mu0_n                    <- solve(t(xx[[1 + q]]) %*% xx[[1 + q]] + emiss_K0[[q]]) %*% (t(xx[[1 + q]]) %*% emiss_c_mu[[i]][[q]] + emiss_K0[[q]] %*% emiss_mu0[[q]][,i])
                    emiss_a_mu_n                   <- (emiss_K0[[q]] + n_subj) / 2
                    emiss_b_mu_n                   <- (emiss_nu[[q]] * emiss_V[[q]][i]) / 2 + (t(emiss_c_mu[[i]][[q]]) %*% emiss_c_mu[[i]][[q]] +
                                                                                                   t(emiss_mu0[[q]][,i]) %*% emiss_K0[[q]] %*% emiss_mu0[[q]][,i] -
                                                                                                   t(emiss_mu0_n) %*% (t(xx[[1 + q]]) %*% xx[[1 + q]] + emiss_K0[[q]]) %*% emiss_mu0_n) / 2
                    emiss_V_mu[[i]][[q]]       <- solve(stats::rgamma(1, shape = emiss_a_mu_n, rate = emiss_b_mu_n))
                    if(all(dim(emiss_V_mu[[i]][[q]]) == c(1,1))){
                        emiss_c_mu_bar[[i]][[q]]	  <- emiss_mu0_n + rnorm(1 + nx[1 + q] - 1, mean = 0, sd = sqrt(diag(as.numeric(emiss_V_mu[[i]][[q]]) * solve(t(xx[[1 + q]]) %*% xx[[1 + q]] + emiss_K0[[q]]))))
                    } else {
                        emiss_c_mu_bar[[i]][[q]]	  <- emiss_mu0_n + rnorm(1 + nx[1 + q] - 1, mean = 0, sd = sqrt(diag(emiss_V_mu[[i]][[q]] * solve(t(xx[[1 + q]]) %*% xx[[1 + q]] + emiss_K0[[q]]))))
                    }
                }

                # Sample subject values  -----------
                for (s in 1:n_subj){

                    trans_subj          <- factor(c(trans[[s]][[i]]))
                    trans_subj          <- factor(trans_subj, levels = paste(c(1:m)[-i]))
                    trans_subj          <- as.numeric(trans_subj)

                    # Sample subject values for gamma using RW Metropolis sampler   ---------
                    gamma_candcov_comb 		<- chol2inv(chol(subj_data[[s]]$gamma_mhess[(1 + (i - 1) * (m - 2)):((m - 2) + (i - 1) * (m - 2)), ] + chol2inv(chol(gamma_V_int[[i]]))))
                    gamma_RWout				<- mnl_RW_once(int1 = gamma_c_int[[i]][s,],
                                                  Obs = trans_subj,
                                                  n_cat = m-1,
                                                  mu_int_bar1 = c(t(gamma_mu_int_bar[[i]]) %*% xx[[1]][s,]),
                                                  V_int1 = gamma_V_int[[i]],
                                                  scalar = gamma_scalar,
                                                  candcov1 = gamma_candcov_comb)
                    gamma[[s]][i,-i]  	    <- PD_subj[[s]][iter, c((n_dep * 2 * m + 1 + (i - 1) * m):(n_dep * 2 * m + (i - 1) * m + m))[-i]] <- gamma_RWout$prob
                    PD_subj[[s]][iter, c((n_dep * 2 * m + 1 + (i - 1) * m):(n_dep * 2 * m + (i - 1) * m + m))[i]] <- 0
                    gamma_naccept[s, i]		<- gamma_naccept[s, i] + gamma_RWout$accept
                    gamma_c_int[[i]][s,]	<- gamma_RWout$draw_int
                    # gamma_int_subj[[s]][iter, c((1 + (i - 1) * (m - 2)):((m - 2) + (i - 1) * (m - 2)))[-i]] <- gamma_c_int[[i]][s,] # CHECK
                    gamma_int_subj[[s]][iter, c((1 + (i - 1) * (m - 2)):((m - 2) + (i - 1) * (m - 2)))] <- gamma_c_int[[i]][s,]

                    if(i == m){
                        delta[[s]] 		<- solve(t(diag(m) - gamma[[s]] + 1), rep(1, m))
                    }
                }
                # Sample subject values for normal emission distribution using Gibbs sampler   ---------

                # population level, conditional probabilities, seperate for each dependent variable
                for(q in 1:n_dep){
                    for (s in 1:n_subj){
                        ss_subj[s] <- t(matrix(cond_y[[s]][[i]][[q]] - emiss_c_mu[[i]][[q]][s,1], nrow = 1) %*%
                                            matrix(cond_y[[s]][[i]][[q]] - emiss_c_mu[[i]][[q]][s,1], ncol = 1))
                        n_cond_y[s]       <- length(cond_y[[s]][[i]][[q]])
                    }
                    emiss_a_resvar_n <- sum(n_cond_y) / 2 + emiss_a0[[q]][i]
                    emiss_b_resvar_n <- (sum(ss_subj) + 2 * emiss_b0[[q]][i]) / 2
                    emiss_c_V[[i]][[q]] <- emiss_var_bar[[q]][iter, i] <- solve(stats::rgamma(1, shape = emiss_a_resvar_n, rate = emiss_b_resvar_n))
                }

                ### sampling subject specific means for the emission distributions, assuming known mean and var, see Lynch p. 244
                for(q in 1:n_dep){
                    emiss_c_V_subj    <- (emiss_V_mu[[i]][[q]] * emiss_c_V[[i]][[q]]) / (2 * emiss_V_mu[[i]][[q]] + emiss_c_V[[i]][[q]])
                    for (s in 1:n_subj){
                        emiss_mu0_subj_n  <- (emiss_V_mu[[i]][[q]] * sum(cond_y[[s]][[i]][[q]]) +  emiss_c_V[[i]][[q]] * c(t(emiss_c_mu_bar[[i]][[q]]) %*% xx[[q+1]][s,])) /
                            (n_cond_y[s] * emiss_V_mu[[i]][[q]] + emiss_c_V[[i]][[q]])
                        emiss[[s]][[q]][i,1] <- PD_subj[[s]][iter, ((q - 1) * m + i)] <- emiss_c_mu[[i]][[q]][s,1] <- rnorm(1, emiss_mu0_subj_n, sqrt(emiss_c_V_subj))
                        emiss[[s]][[q]][i,2] <- PD_subj[[s]][iter, (n_dep * m + (q - 1) * m + i)] <- emiss_c_V[[i]][[q]]
                    }
                }


                #################
                # Obtain hierarchical and mouse specific parameters for duration distribuiton using gibbs sampler ====================================
                #################

                #draw logmu's
                for(s in 1:n_subj){
                    tau2 		<- 1/ ((1/tau2_d_bar[i]) + (1/logsigma2[i]) * sum(sample_path_state[[s]][-n.Dur[s]] == i))
                    logmu[s,i]	<- rnorm(1,
                                        mean = tau2 * ((1/tau2_d_bar[i]) * mu_d_bar[i] + (1/logsigma2[i]) * sum(log(Dur[[s]][-n.Dur[s]][sample_path_state[[s]][-n.Dur[s]] == i]))),
                                        sd = sqrt(tau2))
                    PD_subj[[s]][iter, n_dep*m*2 + m*m + i]			<- logmu[s, i] # adjust index
                }

                # draw mu_d_bar
                s2				<- 1 / ((1/s2_0[i]) + (1/tau2_d_bar[i]) * n_subj)
                mu_d_bar[i] 	<- rnorm(1,
                                      mean = s2 * ((1/s2_0[i]) * d_mu0[i] + (1/tau2_d_bar[i]) * sum(logmu[,i])),
                                      sd = sqrt(s2))

                # draw tau2_d_bar
                a1 <- alpha.tau20[i] + n_subj/2
                b1 <- beta.tau20[i] + (1/2) * (sum((logmu[,i] - mu_d_bar[i])^2))
                tau2_d_bar[i] <- 1 / rgamma(1, shape = a1, rate = b1)

                #draw logsigma
                ss <- numeric(1)
                n.ss <- numeric(1)
                for (s in 1:n_subj){
                    ss <- ss + sum((log(Dur[[s]][-n.Dur[s]][sample_path_state[[s]][-n.Dur[s]] == i]) - logmu[s,i])^2)
                    n.ss <- n.ss + sum(sample_path_state[[s]][-n.Dur[s]] == i)
                }
                c1 <- alpha.sigma20[i] + n.ss/2
                d1 <- beta.sigma20[i] + (1/2) * ss
                logsigma2[i] <- 1 / rgamma(1, shape = c1, rate = d1)
                for (s in 1:n_subj) {
                    PD_subj[[s]][iter, n_dep*m*2 + m*m + m + i] <- logsigma2[i]
                }

            }

        }


        # End of MCMC iteration, save output values --------
        gamma_int_bar[iter, ]				   	<- unlist(lapply(gamma_mu_int_bar, "[",1,))
        gamma_V_int_bar[iter, ] <- unlist(lapply(gamma_V_int, function(e) as.vector(t(e))))
        if(nx[1] > 1){
            gamma_cov_bar[iter, ]      	<- unlist(lapply(gamma_mu_int_bar, "[",-1,))
        }
        if(m == 2) {
            gamma_prob_bar[iter,]                                   <- unlist(gamma_mu_prob_bar)
        } else {
            gamma_prob_bar[iter,(1:m^2)[-seq(1,m*m, m+1)]]			<- unlist(gamma_mu_prob_bar)
            gamma_prob_bar[iter,seq(1,m*m, m+1)]			        <- 0
        }
        for(q in 1:n_dep){
            emiss_mu_bar[[q]][iter, ]	<- as.vector(unlist(lapply(
                lapply(emiss_c_mu_bar, "[[", q), "[",1,)
            ))
            if(nx[1+q] > 1){
                emiss_cov_bar[[q]][iter, ]  <- as.vector(unlist(lapply(
                    lapply(emiss_c_mu_bar, "[[", q), "[",-1,)
                ))
            }
            emiss_varmu_bar[[q]][iter,]	<- as.vector(unlist(sapply(emiss_V_mu, "[[", q)))
        }

        # dwell time parameters
        # duration_bar[iter,] <- c(mu_d_bar, logsigma2, tau2_d_bar)
        dwell_mu_bar[iter, ] <- mu_d_bar # IMPROVE memory use
        dwell_varmu_bar[iter, ] <- tau2_d_bar
        dwell_var_bar[iter, ] <- logsigma2

        if(show_progress == TRUE){
            utils::setTxtProgressBar(pb, iter)
        }
    }
    if(show_progress == TRUE){
        close(pb)
    }

    # End of function, return output values --------
    ctime = proc.time()[3]
    message(paste("Total time elapsed (hh:mm:ss):", hms(ctime-itime)))
    if(return_path == TRUE){
        out <- list(input = list(m = m, n_dep = n_dep, J = J,
                                 burn_in = burn_in, n_subj = n_subj, n_vary = n_vary, dep_labels = dep_labels),
                    PD_subj = PD_subj,
                    gamma_prob_bar = gamma_prob_bar,
                    gamma_int_bar = gamma_int_bar,
                    gamma_cov_bar = gamma_cov_bar,
                    gamma_V_int_bar = gamma_V_int_bar,
                    gamma_int_subj = gamma_int_subj,
                    gamma_naccept = gamma_naccept,
                    emiss_mu_bar = emiss_mu_bar, emiss_varmu_bar = emiss_varmu_bar, emiss_cov_bar = emiss_cov_bar,
                    emiss_var_bar = emiss_var_bar,
                    dwell_mu_bar = dwell_mu_bar, dwell_varmu_bar = dwell_varmu_bar,
                    dwell_var_bar = dwell_var_bar,
                    sample_path = sample_path)
    } else {
        out <- list(input = list(m = m, n_dep = n_dep, J = J,
                                 burn_in = burn_in, n_subj = n_subj, n_vary = n_vary, dep_labels = dep_labels),
                    PD_subj = PD_subj,
                    gamma_prob_bar = gamma_prob_bar,
                    gamma_int_bar = gamma_int_bar,
                    gamma_cov_bar = gamma_cov_bar,
                    gamma_V_int_bar = gamma_V_int_bar,
                    gamma_int_subj = gamma_int_subj,
                    gamma_naccept = gamma_naccept,
                    emiss_mu_bar = emiss_mu_bar, emiss_varmu_bar = emiss_varmu_bar, emiss_cov_bar = emiss_cov_bar,
                    emiss_var_bar = emiss_var_bar,
                    dwell_mu_bar = dwell_mu_bar, dwell_varmu_bar = dwell_varmu_bar,
                    dwell_var_bar = dwell_var_bar)
    }
    class(out) <- append(class(out), "medHMM_cont")
    return(out)
}
