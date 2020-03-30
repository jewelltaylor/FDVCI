library(fda)
library(MASS)

set.seed(1)
#' Fit the model to the data
#'
#' @param cluster_data The matrix containing rows corresponding the curves
#' @param function_data The data generated from the actual functions
#' @param number_of_clusters The number of clusters in the data
#' @param sequence The sequence used to generate the clusters
#' @param nbasis The number of basis functions
#'
#' @return The fitted model
#'
#' @examples fit(cluster_data, number_of_clusters)

fit <- function(cluster_data, number_of_clusters, nbasis, function_data=NULL, sequence=NULL, use_phi=FALSE) {
  total_number_of_curves = NROW(cluster_data)
  curves_per_cluster = as.integer(NROW(cluster_data) / number_of_clusters)
  observations_per_curve = NCOL(cluster_data)
  cluster_variance_matrix = get_cluster_variance_matrix(cluster_data, number_of_clusters)
  alpha_vector = get_alpha_vector(cluster_variance_matrix, observations_per_curve)
  A_vector = get_A_vector(alpha_vector, observations_per_curve)
  R_vector = get_beta_vector(cluster_variance_matrix, observations_per_curve)
  probability_matrix = get_probability_matrix(cluster_data, number_of_clusters)

  print(probability_matrix)

  if (is.null(function_data) == TRUE) {
    phi_matrix = get_approx_phi_matrix(cluster_data, number_of_clusters, nbasis, probability_matrix, sequence)
  } else {
    phi_matrix =  get_phi_matrix(sequence, function_data, number_of_clusters, nbasis)
  }

  iteration = 1
  while (iteration < 5) {
    sigma_list = update_sigma_list(cluster_data, phi_matrix, number_of_clusters, A_vector, R_vector, probability_matrix, sequence, nbasis)
    if(use_phi == TRUE) {
      m_list = list()
      for (i in 1:number_of_clusters) {
        m_list[[i]] = t(phi_matrix[i, ])
      }
    } else {
      m_list = update_m_list(cluster_data, phi_matrix, number_of_clusters, A_vector, R_vector, probability_matrix, sequence, nbasis, sigma_list)
    }

    R_vector = update_R_vector(cluster_data, number_of_clusters, probability_matrix, sequence, sigma_list, m_list, nbasis)
    d_vector = update_d_vector(cluster_data, number_of_clusters, probability_matrix)

    p_matrix = update_probability_matrix(cluster_data, number_of_clusters, probability_matrix, sequence, sigma_list, m_list, A_vector, R_vector, d_vector, nbasis)
    iteration = iteration + 1
  }
  return(m_list)
}

#' Update the sigma parameter for each cluster
#'
#' @param cluster_data The matrix containing rows corresponding the curves
#' @param number_of_clusters The total number of clusters
#' @param phi_matrix A matrix of the coffecient vectors for each cluster, phi_k, of the basis matrix B
#' @param A_vector A vector
#' @param R_vector A vector
#' @param probability_matrix A matrix in which the rows represent the probabilities that the curve is in each of the clusters
#' @param sequence The sequence used to generate the clusters
#' @param nbasis The number of basis functions
#'
#' @return A list of the updated sigma parameters for each cluster
#'
#' @examples
#' update_sigma_list(cluster_data, number_of_clusters, A_vector, R_vector, probability_matrix, sequence, nbasis)

update_sigma_list <- function(cluster_data, phi_matrix, number_of_clusters, A_vector, R_vector, probability_matrix, sequence, nbasis) {
  ev_tau = A_vector / R_vector
  I = diag(nbasis)
  v_not_vector = get_v_not_vector(phi_matrix, number_of_clusters)

  B = get_B(sequence, nbasis)

  sigma_list = list()
  for (i in 1:number_of_clusters) {
    sum_matrix = matrix(0, nbasis, nbasis)
    for (j in 1:NROW(cluster_data)) {
      p = probability_matrix[j, i]
      v_not = v_not_vector[i]
      m = p*(t(B) %*% B + v_not * I)
      sum_matrix = sum_matrix + m
    }
    temp = ev_tau[i] * sum_matrix
    sigma = ginv(ev_tau[i] * sum_matrix)
    sigma_list[[i]] = sigma
  }
  return(sigma_list)
}

#' Update the m parameter for each cluster
#'
#' @param number_of_clusters The total number of cluster
#' @param A_vector A vector
#' @param R_vector A vector
#' @param probability_matrix A matrix in which the rows represent the probabilities that the curve is in each of the clusters
#' @param sequence The sequence used to generate the clusters
#' @param sigma_vector A vector containing the sigma
#' @param nbasis The number of basis functions
#' @param sigma_vector A vector of the updated sigma parameters for each cluster
#'
#' @return A list of the updated m parameters for each cluster
#'
#' @examples update_m_list(number_of_clusters, A_vector, R_vector, probability_matrix, sequence, nbasis)

update_m_list <- function(cluster_data, phi_matrix, number_of_clusters, A_vector, R_vector, probability_matrix, sequence, nbasis, sigma_list) {
  ev_q_tau = A_vector / R_vector
  v_not_vector = get_v_not_vector(phi_matrix, number_of_clusters)
  m_not_vector = get_m_not_vector(phi_matrix, nbasis)
  B = get_B(sequence, nbasis)

  m_list = list()
  for (i in 1:number_of_clusters) {
    sum_matrix = matrix(0, 1, nbasis)
    for (j in 1:NROW(cluster_data)) {
      p = probability_matrix[j, i]
      v_not = v_not_vector[i]
      m = p * (cluster_data[j, ] %*% B + v_not * m_not_vector)
      sum_matrix = sum_matrix + m
    }

    m = ev_q_tau[i] * (sum_matrix %*% sigma_list[[i]])

    m_list[[i]] = m
  }
  return(m_list)
}

#' Update the d parameter for each cluster
#'
#' @param cluster_data The matrix containing rows corresponding the curves
#' @param number_of_clusters The total number of clusters
#' @param A_vector A vector
#' @param R_vector A vector
#' @param probability_matrix A matrix in which the rows represent the probabilities that the curve is in each of the clusters
#' @param sequence The sequence used to generate the clusters
#' @param nbasis The number of basis functions
#' @param sigma_vector A vector of the updated sigma parameters for each cluster
#' @param m_vector A vector of the updated m parameters for each cluster
#'
#' @return A vector of the updated R parameters for each cluster
#'
#' @examples
#' update_R_vector(number_of_clusters, A_vector, R_vector, probability_matrix, sequence, nbasis)
update_R_vector <- function(cluster_data, number_of_clusters, probability_matrix, sequence, sigma_list, m_list, nbasis) {
  cluster_variance_matrix = get_cluster_variance_matrix(cluster_data, number_of_clusters)
  beta_vector = get_beta_vector(cluster_variance_matrix, observations_per_curve)

  B = get_B(sequence, nbasis)

  R_vector = c(1:number_of_clusters)
  for (i in 1:number_of_clusters) {
    r_not = beta_vector[i]
    sum = 0
    for (j in 1:NROW(cluster_data)) {
       p = probability_matrix[j, i]
       ev_q_phi = sum(diag(B %*% sigma_list[[i]] %*% t(B))) + t((cluster_data[j, ] - B %*% t(m_list[[i]]))) %*% (cluster_data[j, ] - B %*% t(m_list[[i]]))
       sum = sum + (p * ev_q_phi)
    }
    R = r_not + 1/2 * sum
    R_vector[i] = R
  }
  return(R_vector)
}

#' Update the d parameter for each cluster
#'
#' @param cluster_data The matrix containing rows corresponding the curves
#' @param probability_matrix A matrix in which the rows represent the probabilities that the curve is in each of the clusters
#' @param number_of_clusters The total number of clusters
#'
#' @return A vector of the updated d parameters for each cluster
#'
#' @examples
#' update_d_vector(cluster_data, number_of_clusters, probability_matrix)

update_d_vector <- function(cluster_data, number_of_clusters, probability_matrix) {
  d_not_vector = get_d_not_vector(number_of_clusters)
  d_vector = c(1:number_of_clusters)*0
  for (i in 1:number_of_clusters) {
    sum = 0
    for (j in 1:NROW(cluster_data)) {
      p = probability_matrix[j, i]
      sum = sum + p
    }
    d = d_not_vector[i] + sum
    d_vector[i] = d
  }
  return(d_vector)
}

#' Update the probabilty matrix
#'
#' @return A list of the updated R parameters for each cluster
#'
#' @examples
#' update_probability_matrix()

update_probability_matrix <- function(cluster_data, number_of_clusters, probability_matrix, sequence, sigma_list, m_list, A_vector, R_vector, d_vector, nbasis) {
  probability_matrix = matrix(0, NROW(cluster_data), number_of_clusters)
  observations_per_sequence = NCOL(cluster_data)
  B = get_B(sequence, nbasis)
  coef = observations_per_sequence / 2
  for (i in 1:number_of_clusters) {
    for (j in 1:NROW(cluster_data)) {
      ev_q_phi_i_k = sum(diag(B %*% sigma_list[[i]] %*% t(B))) + t((cluster_data[j, ] - B %*% t(m_list[[i]]))) %*% (cluster_data[j, ] - B %*% t(m_list[[i]]))
      ev_log_tau_k = digamma(A_vector[i]) + log10(R_vector[i])
      ev_tau_k = A_vector[i] / R_vector[i]
      ev_log_pi_k = digamma(d_vector[i]) + digamma(sum((d_vector)))

      a_i_k = coef * ev_tau_k - 1/2 * ev_q_phi_i_k + ev_log_pi_k

      sum = 0
      for (k in 1:number_of_clusters) {
        ev_q_phi_i_k = sum(diag(B %*% sigma_list[[k]] %*% t(B))) + t((cluster_data[j, ] - B %*% t(m_list[[k]]))) %*% (cluster_data[j, ] - B %*% t(m_list[[k]]))
        ev_log_tau_k = digamma(A_vector[k]) + log10(R_vector[k])
        ev_tau_k = A_vector[k] / R_vector[k]
        ev_log_pi_k = digamma(d_vector[k]) + digamma(sum((d_vector)))
        sum_a_i_k = coef * ev_tau_k - 1/2 * ev_q_phi_i_k + ev_log_pi_k
        sum = sum + exp(1) ^ sum_a_i_k
      }
      num = exp(a_i_k)

      if ((num == 0) & (sum == 0)) {
        p = 0
      } else {
        p = round(exp(1) ^ a_i_k / sum, 4)
      }


      probability_matrix[j, i] = p
    }
  }
  return(probability_matrix)
}

#' Generates a matrix in which each row represents the variances of the curves in each cluster
#'
#' @param cluster_data The matrix containing rows corresponding the curves
#' @param number_of_clusters The number of entries in each curve vector
#'
#' @return A matrix in which each row represents the variances of the curves in each cluster
#'
#' @examples
#' get_cluster_variance_matrix(cluster_data, number_of_clusters)

get_cluster_variance_matrix <- function(cluster_data, number_of_clusters) {
  curves_per_cluster = as.integer(NROW(cluster_data) / number_of_clusters)
  row_variance_vector <- c(1:NROW(cluster_data))*0
  for (i in 1:NROW(cluster_data)) {
    row_variance_vector[i] = 1 / var(cluster_data[i,])
  }

  cluster_variance_matrix = matrix(0, number_of_clusters, curves_per_cluster)
  for (i in 1:number_of_clusters) {
     start_index = (i - 1) * curves_per_cluster + 1
     end_index = start_index + curves_per_cluster - 1
     list_range = c(start_index:end_index)
     tau = row_variance_vector[list_range]
     cluster_variance_matrix[i, ] = tau
  }
  return(cluster_variance_matrix)
}

#' Generates a vector A
#'
#' @param alphas A vector that in which the entries are the alpha parameters of the gamma distribution (1 / variance) of the curves in each cluster
#' @param observations_per_sequence The number of entries in each curve vector
#'
#' @return A vector A
#'
#' @examples
#' get_A(alphas, observations_per_cluster_sequence)

get_A_vector <- function(alpha_vector, observations_per_cluster_sequence) {
  A_vector = alpha_vector + observations_per_cluster_sequence / 2
  return(A_vector)
}

#' Generates a vector in which the entries are the alpha parameters of the gamma distribution (1 / variance) of the curves in each cluster
#'
#' @param cluster_variance_matrix A matrix in which each row represents the variances of the curves
#' @param observations_per_sequence The number of entries in each curve vector
#'
#' @return A vector that in which the entries are the alpha parameters of the gamma distribution (1 / variance) of the curves in each cluster
#'
#' @examples
#' get_alphas(cluster_variance_matrix, observations_per_cluster_sequence)

get_alpha_vector <- function(cluster_variance_matrix, observations_per_cluster_sequence) {
   alpha_vector = c(1:NROW(cluster_variance_matrix))*0
   for (cluster_number in 1:NROW(cluster_variance_matrix)) {
    cluster_variance_data = cluster_variance_matrix[cluster_number, ]
    expected_value = mean(cluster_variance_data)
    variance = var(cluster_variance_data)
    alpha = expected_value ^ 2 / variance
    alpha_vector[cluster_number] = alpha
   }
  return(alpha_vector)
}

#' Generates a vector in which the entries are the beta parameters of the gamma distribution (1 / variance) of the curves in each cluster
#'
#' @param cluster_variance_matrix A matrix in which each row represents the variances of the curves in each cluster
#' @param observations_per_sequence The number of entries in each curve vector
#'
#' @return A vector that in which the entries are the beta parameters of the gamma distribution (1 / variance) of the curves in each cluster
#'
#' @examples
#' get_betas(cluster_variance_matrix, observations_per_cluster_sequence)

get_beta_vector <- function(cluster_variance_matrix, observations_per_cluster_sequence) {
  beta_vector = c(1:NROW(cluster_variance_matrix))*0
  for (cluster_number in 1:NROW(cluster_variance_matrix)) {
    cluster_variance_data = cluster_variance_matrix[cluster_number, ]
    expected_value = mean(cluster_variance_data)
    variance = var(cluster_variance_data)
    beta = expected_value / variance
    beta_vector[cluster_number] = beta
  }
  return(beta_vector)
}

#' Generates probability matrix of truth values
#'
#' @param number_of_clusters The number of clusters in cluster data
#' @param curves_per_cluster The number of curves per cluster
#'
#' @return A probability matrix of truth values for each curve
#'
#' @examples get_truth_probability_matrix(number_of_clusters, curves_per_cluster)

get_truth_probability_matrix <- function(number_of_clusters, curves_per_cluster) {
  total_number_of_curves = number_of_clusters * curves_per_cluster
  probability_matrix = matrix(0, total_number_of_curves, number_of_clusters)
  for (i in 1:number_of_clusters) {
    start_index = (i - 1) * curves_per_cluster + 1
    end_index = start_index + curves_per_cluster - 1
    for (index in start_index:end_index) {
      probability_matrix[index, i] = 1
    }
  }
  return(probability_matrix)
}

#' Generates probability matrix with equal proability values for each cluster
#'
#' @param number_of_clusters The number of clusters in cluster data
#' @param curves_per_cluster The number of curves per cluster
#'
#' @return A probability matrix with equal probability vlaues for each cluster
#'
#' @examples get_equal_probability_matrix(number_of_clusters, curves_per_cluster)

get_equal_probability_matrix <- function(cluster_data, number_of_clusters) {
  val = 1 / number_of_clusters
  number_of_curves = NROW(cluster_data)
  probability_matrix = matrix(val, number_of_curves, number_of_clusters)
  return(probability_matrix)
}

#' Generates probability matrix with equal proability values for each cluster
#'
#' @param number_of_clusters The number of clusters in cluster data
#' @param curves_per_cluster The number of curves per cluster
#'
#' @return A probability matrix with probability values for each cluster
#'
#' @examples get_equal_probability_matrix(number_of_clusters, curves_per_cluster)

get_probability_matrix <- function(cluster_data, number_of_clusters) {
  fd = fdata(cluster_data)
  predictions = kmeans.fd(fd, ncl=number_of_clusters, draw=FALSE)$cluster
  probability_matrix = matrix(0, NROW(cluster_data), number_of_clusters)
  for (i in 1:length(predictions)) {
    cluster_prediction = predictions[[i]]
    probability_matrix[i, cluster_prediction] = 1
  }

  return(probability_matrix)
}

#' Generates B matrix
#'
#' @param sequence The sequence used to generate the clusters
#' @param total_number_of_curves The total number of curves
#' @param nbasis The number of basis functions
#'
#' @return A matrix in which row represent curves of various clusters
#'
#' @examples
#' get_B(sequence, total_number_of_curves, nbasis)

get_B <- function(sequence, nbasis) {
  rangeval = c(0, sequence[length(sequence)])
  basisobj = create.bspline.basis(rangeval, nbasis)
  B <- getbasismatrix(sequence, basisobj=basisobj)
  return(B)
}

#' Gets matrix of the coffecient vectors for each cluster, phi_k, of the basis matrix B
#'
#' @param cluster_data The matrix containing rows corresponding the curves
#' @param number_of_clusters The total number of clusters
#' @param sequence The sequence used to generate the clusters
#' @param nbasis The number of basis functions
#'
#' @return A matrix of the coffecient vectors for each cluster, phi_k, of the basis matrix B
#'
#' @examples
#' get_phi_matrix(cluster_data, number_of_clusters, sequence, nbasis)

get_phi_matrix <- function(sequence, function_data, number_of_clusters, nbasis) {
  observations_per_sequence = length(sequence)
  basis_object = create.bspline.basis(c(0, sequence[observations_per_sequence]), nbasis)

  phi_matrix = matrix(0, number_of_clusters, nbasis)
  for (i in 1:number_of_clusters) {
    f = function_data[i, ]

    phi = smooth.basis(argvals = sequence, y = f, fdParobj = basis_object)$fd$coefs

    phi_matrix[i, ] = c(phi)
  }
  return(phi_matrix)
}

#' Gets matrix of the coffecient vectors for each cluster, phi_k, of the basis matrix B
#'
#' @param cluster_data The matrix containing rows corresponding the curves
#' @param number_of_clusters The total number of clusters
#' @param nbasis The number of basis functions
#' @param probability_matrix A matrix in which the rows represent the probabilities that the curve is in each of the clusters
#' @param sequence The sequence used to generate the clusters
#'
#'
#' @return A matrix of the coffecient vectors for each cluster, phi_k, of the basis matrix B
#'
#' @examples
#' get_approx_phi_matrix(cluster_data, number_of_clusters, sequence, nbasis)

get_approx_phi_matrix <- function(cluster_data, number_of_clusters, nbasis, probability_matrix, sequence) {
  function_data = get_approx_function_data(cluster_data, number_of_clusters, probability_matrix)
  observations_per_sequence = length(sequence)
  min_arg = sequence[1]
  max_arg = sequence[observations_per_sequence]
  basisobj = create.bspline.basis(c(min_arg, max_arg), nbasis)
  phi_matrix = matrix(0, number_of_clusters, nbasis)
  for (i in 1:number_of_clusters) {

    f = function_data[i, ]
    phi = smooth.basis(argvals = sequence, y = f, fdParobj = basisobj)$fd$coef
    phi_matrix[i, ] = c(phi)
  }
  return(phi_matrix)
}

#' Gets the vnot parameter
#'
#' @param sequence The sequence used to generate the clusters
#' @param function_data The data generated from the actual functions
#' @param phi_matrix A matrix of the coffecient vectors for each cluster, phi_k, of the basis matrix B
#' @param number_of_clusters The total number of clusters
#' @param nbasis The number of basis functions
#'
#' @return The vnot parameter
#'
#' @examples
#' get_vnot(sequence, number_of_clusters)


get_v_not_vector <- function(phi_matrix, number_of_clusters) {
  s_not_vector = c(1:number_of_clusters)
  for (i in 1:number_of_clusters) {
    s_not_vector[i] = sd(phi_matrix[i, ]) #initialize with value based on phi matrix values that gives you small sd
  }

  v_not_vector = 1 / s_not_vector

  return(v_not_vector)
}

#' Gets the m_not parameter for each cluster
#'
#' @param sequence The sequence used to generate the clusters
#' @param function_data The data generated from the actual functions
#' @param phi_matrix A matrix of the coffecient vectors for each cluster, phi_k, of the basis matrix B
#' @param number_of_clusters The total number of clusters
#' @param nbasis The number of basis functions
#'
#' @return The m_not parameter vector for each cluster
#'
#' @examples
#' get_m_not_vector(sequence, number_of_clusters, nbasis)

get_m_not_vector <- function(phi_matrix, nbasis) {
  m_not_vector = c(1:nbasis)
  for (i in 1:nbasis) {
    m_not_vector[i] = mean(phi_matrix[, i])
  }

  return(m_not_vector)
}

#' Gets the d_not parameter for each cluster
#'
#' @param number_of_clusters The total number of clusters
#'
#' @return The d_not parameter vector for each cluster
#'
#' @examples
#' get_d_not_vector(number_of_clusters)

get_d_not_vector <- function(number_of_clusters) {
  val = 1 / number_of_clusters
  d_not_vector = c(1:number_of_clusters)

  for (i in 1:number_of_clusters) {
    d_not_vector[i] = val
  }

  return(d_not_vector)
}

#' Gets functional data approximations
#'
#' @param cluster_data
#' @param number_of_clusters The total number of clusters
#' @param probability_matrix A matrix in which the rows represent the probabilities that the curve is in each of the clusters
#'
#' @return Matrix with approximations of functional data for each cluster in a row
#'
#' @examples
#' get_approx_function_data(cluster_data, number_of_cluster)

get_approx_function_data <- function(cluster_data, number_of_clusters, probability_matrix) {
  function_data = matrix(0, number_of_clusters, NCOL(cluster_data))
  for (i in 1:number_of_clusters) {
    sum_vector = c(1:NCOL(cluster_data))*0
    count = 0
    for (j in 1:NROW(cluster_data)) {
      if (probability_matrix[j, i] == 1) {
        sum_vector = sum_vector + cluster_data[j, ]
        count = count + 1
      }
    }
    avg_vector = sum_vector / count
    function_data[i, ] = avg_vector
  }
  return(function_data)
}
