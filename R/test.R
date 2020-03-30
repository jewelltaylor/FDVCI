  source("main.R")
  source("generateData.R")
  library(fda.usc)
  library(weathercan)
  library(tidyverse)
  set.seed(1)

  test <- function() {
    sequence = seq(from=0,to=pi/3, length = 100)
    number_of_clusters = 2
    nbasis = 6
    curves_per_cluster = 10
    B = get_B(sequence, nbasis)

    function_data = get_function_data(sequence, number_of_clusters)
    cluster_data = generate_cluster_data_2(function_data, curves_per_cluster)


    p = fit(cluster_data, number_of_clusters, nbasis=6, function_data = NULL, sequence = sequence, use_phi=FALSE)

    print(p)
  }

  validation_test <- function() {
    number_of_clusters = 2
    nbasis = 6

    boys = growth[[1]]

    girls  = growth[[2]]
    age = growth$age
    cluster_data = matrix(0, NCOL(boys) + NCOL(girls), length(age))
    for (i in 1:NCOL(boys)) {
      cluster_data[i, ] = boys[, i]
    }

    for (i in 1:NCOL(girls)) {
      cluster_data[NCOL(boys) + i, ] = girls[, i]
    }

    p = fit(cluster_data, number_of_clusters, nbasis=6, function_data = NULL, sequence = age, use_phi=FALSE)

    print(p)

  }

# validation_test()

  test_viz <- function() {
    sequence = seq(from=0,to=pi/3, length = 100)
    number_of_clusters = 2
    nbasis = 6
    curves_per_cluster = 10
    B = get_B(sequence, nbasis)

    function_data = get_function_data(sequence, number_of_clusters)
    cluster_data = generate_cluster_data_2(function_data, curves_per_cluster)


    m_list = fit(cluster_data, number_of_clusters, nbasis=6, function_data = NULL, sequence = sequence, use_phi=FALSE)

    dev.new(width=5, height=4)
    matplot(sequence, t(cluster_data), type="l",col="gray",lty=3, xlab="x", ylab="f(x)", main="Simulation Results (Hard)")

    lines(sequence, function_data[1, ],col="red", xlab="x", ylab="f(x)")


    lines(sequence, function_data[2, ], col="green")

    lines(sequence, B %*% t(m_list[[1]]), col="blue")

    lines(sequence, B %*% t(m_list[[2]]), col="black")

    legend(x="topright", legend=c("True Function 1", "True Function 2", "Estimated Function 1", "Estimated Function 2", "Functional Data"), col=c("red", "green", "blue", "black", "grey"), lty=1, lwd=1)

  }

  validation_test2 <- function() {
    cluster_data = matrix(0, 4, 365)
    sequence = seq(1, 365)
    number_of_clusters = 4
    nbasis = 6
    cluster_data[1, ] = CanadianWeather$dailyAv[,,1][,26]
    cluster_data[2, ] = CanadianWeather$dailyAv[,,1][,35]

    p = fit(cluster_data, number_of_clusters, nbasis=6, function_data = NULL, sequence = sequence, use_phi=FALSE)

    print(p)
  }

  #validation_test2()


  validation_test3 <- function() {
    #print((stations_search("VANCOUVER", interval = "day")))

    res_temp = as.vector(weather_dl(station_ids = 1776, start = "1994-01-01", end = "2003-12-31", interval = "day")$max_temp)
    tor_temp = as.vector(weather_dl(station_ids = 5051, start = "1994-01-01", end = "2003-12-31", interval = "day")$max_temp)
    #print(van_temp[1:365])

   # print(van_temp[1:10])
    cl = matrix(0, 18, 31)
    count = 1
    for (i in seq(1, 3285, 365)) {
      upperbound1 = i + 364
      res_temp_curr = replace_na(res_temp[i:upperbound1], 0)[187:217]
      print(res_temp_curr)
      cl[count, ] = res_temp_curr

      count = count + 1
    }

    for (j in seq(1, 3285, 365)) {
      upperbound2 = j + 364
      tor_temp_curr = replace_na(tor_temp[j:upperbound2], 0)[187:217]
      cl[count, ] = tor_temp_curr
      count = count + 1
    }

    sequence = seq(1, 31)
    par(mar=c(1,1,1,1))

    matplot(sequence, t(cl), type="l",col="gray",lty=3)
    B = get_B(sequence, 10)

    view(cl)
    m = fit(cl, 2, nbasis=10, function_data = NULL, sequence = sequence, use_phi=FALSE)
    dev.new(width=5, height=4)
    matplot(sequence, t(cl), type="l",col="gray",lty=3, xlab="Day", ylab="Temperature", main="Validation Results")
    lines(sequence, B %*% t(m[[1]]), col="red")
    lines(sequence, B %*% t(m[[2]]), col="blue")
    #title(main="Valid")
    legend(x="center", legend=c("Estimated Function Vancouver", "Estimated Function Resolute", "Functional Data"), col=c("blue", "red", "grey"), lty=1, lwd=1)

  }

  validation_test3()

