# load necessary packages
library(tidyverse)
library(combinat)
library(mvtnorm)
library(gtools)

# function that computes the probability of persistence
# input: A = initial community structure (interaction matrix); B = initial community structure (interaction matrix);
#        exp = the scale of grid approximation algorithm; constraint = two extreme cases
# output: the probability of persistence from A to B given the constraint
compute_prob_persistence <- function(A, B, exp, constraint, persistence_type = "stability") {

  # function that computes the size of the domain of feasibility
  # input: alpha = community structure (interaction matrix)
  # output: the size of the domain of feasibility
  Omega <- function(alpha) {
    S <- nrow(alpha)
    Sigma <- tryCatch(solve(t(alpha) %*% alpha))
    m <- matrix(0, S, 1)
    a <- matrix(0, S, 1)
    b <- matrix(Inf, S, 1)
    d <- pmvnorm(lower = rep(0, S), upper = rep(Inf, S), mean = rep(0, S), sigma = Sigma)
    out <- d[1]^(1 / S)
    out
  }

  # function that computes the spanning vectors of the domain of feasibility
  # input: alpha = community structure (interaction matrix)
  # output: the spanning vectors of the domain of feasibility
  span_vectors <- function(alpha) {
    Span <- matrix(0, ncol = ncol(alpha), nrow = nrow(alpha))
    for (k in 1:ncol(alpha)) {
      Span[, alpha] <- -alpha[, k]
    }
    Span
  }

  # function that checks dynamical stability
  # input: A = community structure (interaction matrix); N = equilibrium species abundance
  # output: whether community structure A and speces abundance N makes a dynamically stable system
  criteria_stability_N <- function(A, N) {
    N <- as.vector(N)
    num <- length(N)
    sta_mat <- diag(N) %*% A
    numd <- sum(Re(eigen(sta_mat)$values) < -1e-6)
    if (numd == num) {
      state <- 1
    } else {
      state <- 0
    }

    state
  }

  criteria_permanence_N <- function(A, N, r) {
    N <- as.vector(N)
    num <- length(N)
    criteria_2 <- if_else(det(-A) > 0, 1, 0)
    for (i in 1:3) {
      sub_matrix <- A[-i, -i]
      sub_r <- r[-i]
      if (sub_matrix[1, 2] > 0 & sub_matrix[2, 1] > 0) {
        criteria_3 <- if_else(sub_matrix[1, 1] * sub_matrix[2, 2] > sub_matrix[1, 2] * sub_matrix[2, 1], 1, 0)
      } else if (sub_matrix[1, 2] < 0 & sub_matrix[2, 1] < 0) {
        criteria_3 <- if_else(sum(which(solve(sub_matrix, -sub_r) < 0)) == 0, 1, 0)
      } else {
        criteria_3 <- 1
      }
    }
    criteria_2 * criteria_3
  }

  criteria_persistence_N <- function(persistence_type, A, N, r) {
    if (persistence_type == "stability") {
      persistence <- criteria_stability_N(A, N)
    }
    if (persistence_type == "permanence") {
      persistence <- criteria_permanence_N(A, N, r)
    }
    persistence
  }

  # grid approximation algorithm to evenly sample the feasible equilibria inside the feasibility domain
  # input: S = the size of the ecological community; steps = the level of precision
  # output: the sampled feasible equilibria
  sample_points <- function(S, steps) {
    fixed <- 1

    points <- matrix(0, nrow = 1, ncol = S)
    colnames(points) <- paste("a_", 1:S)

    for (i in c(1:((steps - 2) / (S - 2)))) {
      fixed <- i
      up <- c(1:(steps - fixed * (S - 2) - 1))
      down <- c(steps - fixed * (S - 2) - up)
      len_row <- (length(up))

      make_points <- cbind(matrix(up, ncol = 1, nrow = len_row), matrix(down, ncol = 1, nrow = len_row), matrix(fixed, ncol = (S - 2), nrow = len_row))
      points <- rbind(points, make_points)
    }

    points <- points / steps

    initial_points <- nrow(points)

    orders <- permn(1:S)
    orders <- orders[-1]

    for (val in orders) {
      new_points <- points[1:initial_points, val]
      points <- rbind(points, new_points)
    }

    points <- unique(points,
      incomparables = FALSE, MARGIN = 1,
      fromLast = FALSE
    )
    points <- points[-1, ]
    return(points)
  }

  num <- nrow(B)
  test_rs <- sample_points(num, exp)

  if (constraint == "without") {
    point <- c()
    for (j in 1:nrow(test_rs)) {
      r_B <- -B %*% matrix(test_rs[j, ], ncol = 1)
      N_B <- solve(B, -r_B)
      point[j] <- if_else(criteria_persistence_N(persistence_type, B, N_B, r_B) == 1, 1, 0)
    }
    prob <- mean(point, na.rm = T) * Omega(B)
  }

  if (constraint == "with") {
    point <- c()
    for (j in 1:nrow(test_rs)) {
      r_A <- -A %*% matrix(test_rs[j, ], ncol = 1)
      N_A <- solve(A, -r_A)
      if (criteria_persistence_N(persistence_type, A, N_A, r_A) == 1) {
        N_B <- solve(B, -r_A)
        r_B <- -B %*% N_B
        if (sum(N_B > 0) == num) {
          point <- c(point, if_else(criteria_persistence_N(persistence_type, B, N_B, r_B) == 1, 1, 0))
        }
        else {
          point <- c(point, 0)
        }
      }
    }
    prob <- if_else(length(point) > 0, mean(point, na.rm = T), 0)
  }

  prob
}

# function that computes the probability of transition
# input: A = initial community structure (interaction matrix)
# output: the probability of transition of A
compute_prob_transition <- function(A) {
  labels_transitioned <- compute_transition_paths(A)
  1 / length(labels_transitioned)
}

# function that computes the transition paths from A given 1-point transition
# input: A = initial community structure (interaction matrix)
# output: the labels of network topology that A can transition to
compute_transition_paths <- function(A) {
  # get all isomorphism representations of the same topology of the network structure
  # input: AllType_Inte = the set of all possible topologies
  # output: all isomorphism representations for a topology
  generate_isomorphism_matrix <- function(AllType_Inte) {
    iso_class <- list()
    for (j in 1:length(AllType_Inte)) {
      iso_class_Inte <- list()
      perm <- gtools::permutations(n = 3, r = 3, v = c(1:3))
      A <- AllType_Inte[[j]]
      for (i in 1:6) {
        iso_class_Inte[[i]] <- A[perm[i, ], perm[i, ]]
      }
      iso_class[[j]] <- iso_class_Inte
    }
    iso_class
  }

  # identify the lebales of network structures from a given topology after 1-point transition
  # input: A = original network topology, iso_class = all isomorphism representations
  # output: all isomorphism representations for a topology
  search_iso_class <- function(A, iso_class) {
    state_class <- c()
    for (i in 1:length(iso_class)) {
      iso_class_Inte <- iso_class[[i]]
      state <- 0
      for (j in 1:6) {
        if (sum(abs(A - iso_class_Inte[[j]])) == 0) state <- 1
      }
      state_class[i] <- state
    }
    which(state_class == 1)
  }

  # get the number of transitions after 1-point transition
  # input: A = original network topology
  # output: all isomorphism representations for a topology
  generate_mutation <- function(A) {
    pos <- list()
    l <- 1
    diag(A) <- 0
    A[A > 0] <- 1
    A[A < 0] <- -1
    for (i in 1:3) {
      for (j in 1:2) {
        B <- A
        B[i, j] <- -A[i, j]
        pos[[l]] <- B
        l <- l + 1
        B <- A
        B[j, i] <- -A[j, i]
        pos[[l]] <- B
        l <- l + 1
        # }
      }
    }
    Allmutation <- c()
    for (i in 1:(l - 1)) {
      Allmutation[i] <- search_iso_class(pos[[i]], iso_class)
    }

    unique(Allmutation)
  }

  # get the number of transitions after 1-point transition
  # input: A = original network topology
  # output: all isomorphism representations for a topology
  identify_Interaction_type <- function(A) {
    pos <- list()
    l <- 1
    diag(A) <- 0
    A[A > 0] <- 1
    A[A < 0] <- -1
    for (i in 1:3) {
      for (j in 1:2) {
        B <- A
        B[i, j] <- A[i, j]
        pos[[l]] <- B
        l <- l + 1
        B <- A
        B[j, i] <- A[j, i]
        pos[[l]] <- B
        l <- l + 1
        # }
      }
    }
    Allmutation <- c()
    for (i in 1:(l - 1)) {
      Allmutation[i] <- search_iso_class(pos[[i]], iso_class)
    }

    unique(Allmutation)
  }

  AllType_Inte <- Generate_AllType_Inte()
  iso_class <- generate_isomorphism_matrix(AllType_Inte)
  label_original <- identify_Interaction_type(A)
  labels_transitioned <- setdiff(generate_mutation(AllType_Inte[[label_original]]), label_original)

  labels_transitioned
}

# Generate all possible network topologies
# input: none
# output: the list of all possible network topologies
# Note: Case label is from Milo et al., Science, 2003
Generate_AllType_Inte <- function() {
  AllType_Inte <- list()
  # Case 1
  AllType_Inte[[1]] <- matrix(c(
    0, -1, 0,
    1, 0, 1,
    0, -1, 0
  ), ncol = 3, byrow = TRUE) # Exploitative competition
  # Case 2
  AllType_Inte[[2]] <- matrix(c(
    0, -1, 0,
    1, 0, -1,
    0, 1, 0
  ), ncol = 3, byrow = TRUE) # Tri-trophic chain
  # Case 3
  AllType_Inte[[3]] <- matrix(c(
    0, -1, 0,
    1, 0, 1,
    0, 1, 0
  ), ncol = 3, byrow = TRUE)
  AllType_Inte[[4]] <- matrix(c(
    0, -1, 0,
    1, 0, -1,
    0, -1, 0
  ), ncol = 3, byrow = TRUE)
  # Case 4
  AllType_Inte[[5]] <- matrix(c(
    0, -1, -1,
    1, 0, 0,
    1, 0, 0
  ), ncol = 3, byrow = TRUE) # Apparent competition
  # Case 5
  AllType_Inte[[6]] <- matrix(c(
    0, -1, -1,
    1, 0, 1,
    1, -1, 0
  ), ncol = 3, byrow = TRUE) # Omnivory
  # Case 6
  AllType_Inte[[7]] <- matrix(c(
    0, -1, -1,
    1, 0, 1,
    1, 1, 0
  ), ncol = 3, byrow = TRUE)
  AllType_Inte[[8]] <- matrix(c(
    0, -1, -1,
    1, 0, -1,
    1, -1, 0
  ), ncol = 3, byrow = TRUE)
  # Case 7
  AllType_Inte[[9]] <- matrix(c(
    0, 1, 0,
    -1, 0, 1,
    0, 1, 0
  ), ncol = 3, byrow = TRUE)
  AllType_Inte[[10]] <- matrix(c(
    0, 1, 0,
    -1, 0, -1,
    0, -1, 0
  ), ncol = 3, byrow = TRUE)
  # Case 8
  AllType_Inte[[11]] <- matrix(c(
    0, 1, 0,
    1, 0, 1,
    0, 1, 0
  ), ncol = 3, byrow = TRUE)
  AllType_Inte[[12]] <- matrix(c(
    0, -1, 0,
    -1, 0, 1,
    0, 1, 0
  ), ncol = 3, byrow = TRUE)
  AllType_Inte[[13]] <- matrix(c(
    0, -1, 0,
    -1, 0, -1,
    0, -1, 0
  ), ncol = 3, byrow = TRUE)
  # Case 9
  AllType_Inte[[14]] <- matrix(c(
    0, 1, -1,
    -1, 0, 1,
    1, -1, 0
  ), ncol = 3, byrow = TRUE) # Rock-Scissor Game
  # Case 10
  AllType_Inte[[15]] <- matrix(c(
    0, 1, -1,
    1, 0, 1,
    1, -1, 0
  ), ncol = 3, byrow = TRUE)
  AllType_Inte[[16]] <- matrix(c(
    0, -1, -1,
    -1, 0, 1,
    1, -1, 0
  ), ncol = 3, byrow = TRUE)
  # Case 11
  AllType_Inte[[17]] <- matrix(c(
    0, 1, -1,
    1, 0, -1,
    1, 1, 0
  ), ncol = 3, byrow = TRUE)
  AllType_Inte[[18]] <- matrix(c(
    0, -1, -1,
    -1, 0, -1,
    1, 1, 0
  ), ncol = 3, byrow = TRUE)
  # Case 12
  AllType_Inte[[19]] <- matrix(c(
    0, 1, -1,
    1, 0, 1,
    1, 1, 0
  ), ncol = 3, byrow = TRUE)
  AllType_Inte[[20]] <- matrix(c(
    0, 1, -1,
    1, 0, -1,
    1, -1, 0
  ), ncol = 3, byrow = TRUE)
  AllType_Inte[[21]] <- matrix(c(
    0, -1, -1,
    -1, 0, 1,
    1, 1, 0
  ), ncol = 3, byrow = TRUE)
  AllType_Inte[[22]] <- matrix(c(
    0, -1, -1,
    -1, 0, -1,
    1, -1, 0
  ), ncol = 3, byrow = TRUE)
  # Case 13
  AllType_Inte[[23]] <- matrix(c(
    0, 1, 1,
    1, 0, 1,
    1, 1, 0
  ), ncol = 3, byrow = TRUE)
  AllType_Inte[[24]] <- matrix(c(
    0, 1, 1,
    1, 0, -1,
    1, -1, 0
  ), ncol = 3, byrow = TRUE)
  AllType_Inte[[25]] <- matrix(c(
    0, 1, -1,
    1, 0, -1,
    -1, -1, 0
  ), ncol = 3, byrow = TRUE)

  AllType_Inte[[26]] <- matrix(c(
    0, -1, -1,
    -1, 0, -1,
    -1, -1, 0
  ), ncol = 3, byrow = TRUE)

  AllType_Inte
}

# Sample the switch probability from all topologies
# input: strength = the standard deviation of the interspecific interaction strength
# output: the list of all possible network topologies
sample_switch_probability <- function(strength = 1, persistence_type = "stability") {
  pb$tick()$print()

  sample_switch_probability_single <- function(label_original, strength, persistence_type) {
    generate_random_matrix <- function(topology, strength) {
      nspp <- nrow(topology)
      alpha <- topology * abs(rnorm(nspp, 0, strength))
      diag(alpha) <- -1

      alpha
    }

    determine_pair_type <- function(int) {
      type <- c(0, 0, 0)
      names(type) <- c("(-,-)", "(+,-)", "(+,+)")
      if (int[1] < 0 & int[2] < 0) type[1] <- 1
      if (int[1] < 0 & int[2] > 0) type[2] <- 1
      if (int[1] > 0 & int[2] < 0) type[2] <- 1
      if (int[1] > 0 & int[2] > 0) type[3] <- 1
      type
    }

    detemine_type_of_transition <- function(A, B) {
      count_interaction_types <- function(A) {
        A[A > 0] <- 1
        A[A < 0] <- -1
        determine_pair_type(c(A[1, 2], A[2, 1])) + determine_pair_type(c(A[1, 3], A[3, 1])) + determine_pair_type(c(A[2, 3], A[3, 2]))
      }

      int_A <- count_interaction_types(A)
      int_B <- count_interaction_types(B)

      # type_from <- names(which(int_A-int_B==1))
      # type_to <- names(which(int_A-int_B==-1))
      type_from <- which(int_A - int_B == 1)
      type_to <- which(int_A - int_B == -1)

      paste(type_from, type_to, sep = "->")
    }

    generate_transformed_alpha <- function(alpha_original, type_of_transition) {
      determine_change <- function(pair, alpha) {
        type <- 0
        int <- c(alpha[pair[1], pair[2]], alpha[pair[2], pair[1]])
        if (int[1] < 0 & int[2] < 0) type <- 1
        if (int[1] < 0 & int[2] > 0) type <- 2
        if (int[1] > 0 & int[2] < 0) type <- 2
        if (int[1] > 0 & int[2] > 0) type <- 3
        if_else(str_sub(type_of_transition, 1, 1) == type, 1, 0)
      }
      make_change <- function(pair, alpha) {
        from <- str_sub(type_of_transition, 1, 1)
        to <- str_sub(type_of_transition, 4, 4)
        if (to == "1") {
          alpha[pair[1], pair[2]] <- -1 * abs(alpha[pair[1], pair[2]])
          alpha[pair[2], pair[1]] <- -1 * abs(alpha[pair[2], pair[1]])
        }
        if (to == "2") {
          if (runif(1) < .5) {
            alpha[pair[1], pair[2]] <- -1 * abs(alpha[pair[1], pair[2]])
            alpha[pair[2], pair[1]] <- 1 * abs(alpha[pair[2], pair[1]])
          } else {
            alpha[pair[1], pair[2]] <- 1 * abs(alpha[pair[1], pair[2]])
            alpha[pair[2], pair[1]] <- -1 * abs(alpha[pair[2], pair[1]])
          }
        }
        if (to == "3") {
          alpha[pair[1], pair[2]] <- 1 * abs(alpha[pair[1], pair[2]])
          alpha[pair[2], pair[1]] <- 1 * abs(alpha[pair[2], pair[1]])
        }
        alpha
      }

      pairs <- list(c(1, 2), c(1, 3), c(2, 3))
      pair_change <- pairs %>%
        map_dbl(~ determine_change(., alpha_original)) %>%
        {
          which(. == 1)
        } %>%
        sample(1) %>%
        {
          pairs[[.]]
        }

      make_change(pair_change, alpha_original)
    }

    topology_original <- AllType_Inte[[label_original]]
    alpha_original <- generate_random_matrix(topology_original, strength)

    label_transformed <- compute_transition_paths(AllType_Inte[[label_original]])
    types_of_transition <- label_transformed %>%
      map_chr(~ detemine_type_of_transition(alpha_original, AllType_Inte[[.x]]))

    alpha_transfomed <- types_of_transition %>%
      map(~ generate_transformed_alpha(alpha_original, .))

    prob_transition <- compute_prob_transition(alpha_original)

    prob_persistence_with <- alpha_transfomed %>%
      map_dbl(~ compute_prob_persistence(alpha_original, ., exp = 50, constraint = "with", persistence_type))

    prob_persistence_without <- alpha_transfomed %>%
      map_dbl(~ compute_prob_persistence(alpha_original, ., exp = 50, constraint = "without", persistence_type))

    tibble(
      from = label_original,
      to = label_transformed,
      type_of_transition = types_of_transition,
      prob_transition = prob_transition,
      prob_persistence_with = prob_persistence_with,
      prob_persistence_without = prob_persistence_without,
      persistence_type = persistence_type
    )
  }

  map_dfr(1:length(AllType_Inte), ~ sample_switch_probability_single(., strength, persistence_type))
}
