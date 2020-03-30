set.seed(1)

MAX_MEAN = 10
MAX_STD = 1
MAX_NOISE = .2

#' Generates a specified number of clusters based on sequence
#'
#' @param sequence The sequence used to generate the clusters
#' @param number_of_clusters The number of clusters to be generated
#' @param curves_per_cluster The number of curver per cluster
#'
#' @return A matrix in which row represent curves of various clusters
#'
#' @examples
#' generate_cluster_data(sequence, number_of_clusters, curves_per_cluster)

generate_cluster_data <- function(sequence, number_of_clusters, curves_per_cluster) {
  length_of_sequence = length(sequence)
  Y = matrix(0, number_of_clusters * curves_per_cluster, length_of_sequence)
  for (i in 1:number_of_clusters) {
    for (j in 1:curves_per_cluster) {
      b <- runif(1, -1/2, 1/2)
      f <- b + i*5 + sin(pi*sequence*i)
      y <- f + rnorm(length(sequence),0,0.4^2)
      index = (i - 1) * curves_per_cluster + j
      Y[index, ] = y
    }
  }
  return(Y)
}

#' Generates the data corresponding the actual functions
#'
#' @param sequence The sequence used to generate the clusters
#' @param number_of_clusters The number of clusters to be generated
#' @param curves_per_cluster The number of curver per cluster
#'
#' @return A matrix in which row represent curves of various clusters
#'
#' @examples
#' get_function_data(sequence, number_of_clusters, curves_per_cluster)

get_function_data <- function(sequence, number_of_clusters) {
  length_of_sequence = length(sequence)
  F = matrix(0, number_of_clusters, length_of_sequence)
  for (i in 1:number_of_clusters) {
    b <- runif(1, -1/2, 1/2)
    f <- b + i*2.5 + sin(pi*sequence*i)
    F[i, ] = f
  }
  return(F)
}

#' Generates a specified number of clusters based on sequence
#'
#' @param sequence The sequence used to generate the clusters
#' @param number_of_clusters The number of clusters to be generated
#' @param curves_per_cluster The number of curver per cluster
#'
#' @return A matrix in which row represent curves of various clusters
#'
#' @examples
#' generate_cluster_data_2(function_daya, curves_per_cluster)

generate_cluster_data_2 <- function(function_data, curves_per_cluster) {
  cluster_data = matrix(0, NROW(function_data) * curves_per_cluster, NCOL(function_data))
  count = 1
  for(i in 1:NROW(function_data)) {
    for(j in 1:curves_per_cluster) {
      noise = runif(1) * MAX_NOISE
      sequence = function_data[i, ] + rnorm(NCOL(function_data),0,noise)
      cluster_data[count, ] = sequence
      count = count + 1
    }
  }
  return(cluster_data)
}

#' Generates the data corresponding the actual functions
#'
#' @param sequence_length The length of each sequence (curve)
#' @param number_of_clusters The number of clusters to be generated
#'
#' @return A matrix in which row represent curves of various clusters
#'
#' @examples
#' get_function_data_2(sequence_length, number_of_clusters)

get_function_data_2 <- function(sequence_length, number_of_clusters) {
  function_data = matrix(0, number_of_clusters, sequence_length)
  for (i in 1:number_of_clusters) {
    sequence_mean = runif(1) * MAX_MEAN
    sequence_std = runif(1) * MAX_STD
    sequence = rnorm(n=sequence_length, mean=MAX_MEAN, sd=MAX_STD)
    function_data[i, ] = sequence
  }
  return(function_data)

}


