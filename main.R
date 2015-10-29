library(lpSolve)

# be sure to have NFL_Stadium_Locations.csv and load_data.R
# in same directory
source('load_data.R')

out <- load_distance_matrix()
cost_matrix <- matrix(unlist(out[1]), ncol=32, byrow=T)
team_names <- out[2]

num_teams <- 32
num_weeks <- 16
num_cities <- 32

#### for now ###
cost_matrix <- cost_matrix[1:8, 1:8]
num_teams <- 8
num_weeks <- 4
num_cities <- 8
################

nfl.obj <- vector(length = (num_cities^2 * num_teams * num_weeks +
                              num_teams * num_cities * num_weeks + 
                              num_cities^2 * num_teams * num_weeks))

# fill in objective function values
index <- 1
for (L in 1:num_cities) {
  for (M in 1:num_cities) {
    for (W in 1:num_weeks) {
      C <- cost_matrix[L, M]
      for (team in 1:num_teams) {
        nfl.obj[index] <- C
        index <- index + 1
      }
    }
  }
}

start_index_x <- num_cities^2 * num_teams * num_weeks + 1

start_index_y <- num_cities^2 * num_teams * num_weeks + 
  num_cities * num_teams * num_weeks + 1

# add x values to objective function 
# (adding zeros because x doesn't play into cost)
# start at number of v values + 1
for (i in start_index_x:(start_index_y - 1)) {
  nfl.obj[i] <- 0
}

# add y values to objective function 
# (adding zeros because y doesn't play into cost)
# start at number of v values + number of x values + 1
for (i in start_index_y:(start_index_y + num_cities * num_teams^2 * num_weeks - 1)) {
  nfl.obj[i] <- 0
}

# build function that maps (L, M, W, team) to an index.
# ranges: 1-32, 1-32, 1-16, 1-32
# 1, 1, 1, 1 --> 1,
# 1 ,1, 1, 2 --> 2,
# 1, 1, 2, 1 --> 32 + 1 = 33,
# 1, 2, 1, 1 --> 16*32 + 1 = 513,
# (L-1)*32*16*32 + (M-1)*16*32 + (W-1)*32 + team
v_index_map <- function(team, L, M, W) {
  index <- (team-1)*num_cities*num_teams*num_weeks + 
    (L-1)*num_teams*num_weeks + 
    (M-1)*num_weeks + 
    W
  return(index)
}

# Same map as before, but stick x values below v
# values in rhs vector. Hence offset by the number
# of v values.
x_index_map <- function(team, L, week) {
  index <- (team-1)*num_cities*num_weeks +
    (L-1)*num_weeks +
    week + 
    v_index_map(num_teams, num_cities, num_cities, num_weeks) # start where v's left off
  return(index)
}

# Again, same map as before, but stick y values below x
# values in rhs vector. Hence offset by the number of x
# and v values.
y_index_map <- function(i, j, L, week) {
  index <- (i-1)*num_teams*num_cities*num_weeks +
    (j-1)*num_cities*num_weeks +
    (L-1)*num_weeks +
    week + 
    x_index_map(num_teams, num_cities, num_weeks) # start where x's left off
  return(index)
}

rhs_length <- y_index_map(num_teams, num_teams, num_cities, num_weeks)

### constraint 1 ###
# V_tlmw - x_tlw - x_tm(w+1) >= -1
# For all teams, city from, city to, and all but last week.
# Number of rows is 32*32*32*15.
const_1.mat <- matrix(rep(0, times=num_teams*num_cities^2*(num_weeks-1)*rhs_length),
                      ncol=rhs_length)
row <- 1
for (team in 1:num_teams) {
  for (L in 1:num_cities) {
    for (M in 1:num_cities) {
      for (W in 1:(num_weeks - 1)) {
        const_1.mat[row, v_index_map(team, L, M, W)] <- 1
        const_1.mat[row, x_index_map(team, L, W)] <- -1
        const_1.mat[row, x_index_map(team, L, W + 1)] <- -1
        row <- row + 1
      }
    }
  }
}
const_1.dir <- rep('>=', times=nrow(const_1.mat))
const_1.rhs <- rep(-1, times=nrow(const_1.mat))

### constraint 2 ###
# V_tlmw - x_tlw <= 0
# For all teams, city from, city to, and all weeks.
# Number of rows is 32*32*32*16.
const_2.mat <- matrix(rep(0, times=num_teams*num_cities^2*num_weeks*rhs_length),
                      ncol=rhs_length)
row <- 1
for (team in 1:num_teams) {
  for (L in 1:num_cities) {
    for (M in 1:num_cities) {
      for (W in 1:num_weeks) {
        const_2.mat[row, v_index_map(team, L, M, W)] <- 1
        const_2.mat[row, x_index_map(team, L, W)] <- -1
        row <- row + 1
      }
    }
  }
}
const_2.dir <- rep('<=', times=nrow(const_2.mat))
const_2.rhs <- rep(0, times=nrow(const_2.mat))


### constraint 3 ###
# V_tlmw - x_tmw <= 0
# For all teams, city from, city to, and all weeks.
# Number of rows is 32*32*32*16.
const_3.mat <- matrix(rep(0, times=num_teams*num_cities^2*num_weeks*rhs_length),
                      ncol=rhs_length)
row <- 1
for (team in 1:num_teams) {
  for (L in 1:num_cities) {
    for (M in 1:num_cities) {
      for (W in 1:num_weeks) {
        const_3.mat[row, v_index_map(team, L, M, W)] <- 1
        const_3.mat[row, x_index_map(team, M, W)] <- -1
        row <- row + 1
      }
    }
  }
}
const_3.dir <- rep('<=', times=nrow(const_3.mat))
const_3.rhs <- rep(0, times=nrow(const_3.mat))


### constraint 4 ###
# 16 constraints
# for all 1<=W<=16
# sum over m of x_mmw = 16
const_4.mat <- matrix(rep(0, times=num_weeks*rhs_length),
                      ncol=rhs_length)
row <- 1
for (W in 1:num_weeks) {
  for (team in 1:num_teams) {
    const_4.mat[row, x_index_map(team, team, W)] <- 1
  }
  row <- row + 1
}

const_4.dir <- rep('=', times=num_weeks)
const_4.rhs <- rep(num_weeks, times=num_weeks)


### constraint 5 ###
# (1)
# y_ijmw - x_imw - x_jmw >= -1
const_5.mat_1 <- matrix(rep(0, times=num_teams*(num_teams-1)*num_cities*num_weeks*rhs_length),
                        ncol=rhs_length)
row <- 1
for (i in 1:num_teams) {
  for (j in 1:num_teams) {
    if (i != j) {
      for (M in 1:num_cities) {
        for (W in 1:num_weeks) {
          const_5.mat_1[row, y_index_map(i, j, M, W)] <- 1
          const_5.mat_1[row, x_index_map(i, M, W)] <- -1
          const_5.mat_1[row, x_index_map(j, M, W)] <- -1
          row <- row + 1
        }
      }
    }
  }
}

const_5.dir_1 <- rep('>=', times=nrow(const_5.mat_1))
const_5.rhs_1 <- rep(-1, times=nrow(const_5.mat_1))

# (2)
# y_ijmw - x_imw <= 0
const_5.mat_2 <- matrix(rep(0, times=num_teams*(num_teams-1)*num_cities*num_weeks*rhs_length),
                        ncol=rhs_length)
row <- 1
for (i in 1:num_teams) {
  for (j in 1:num_teams) {
    if (i != j) {
      for (M in 1:num_cities) {
        for (W in 1:num_weeks) {
          const_5.mat_2[row, y_index_map(i, j, M, W)] <- 1
          const_5.mat_2[row, x_index_map(i, M, W)] <- -1
          row <- row + 1
        }
      }
    }
  }
}

const_5.dir_2 <- rep('<=', times=nrow(const_5.mat_2))
const_5.rhs_2 <- rep(0, times=nrow(const_5.mat_2))

# (3)
# y_ijmw - x_jmw <= 0
const_5.mat_3 <- matrix(rep(0, times=num_teams*(num_teams-1)*num_cities*num_weeks*rhs_length),
                        ncol=rhs_length)
row <- 1
for (i in 1:num_teams) {
  for (j in 1:num_teams) {
    if (i != j) {
      for (M in 1:num_cities) {
        for (W in 1:num_weeks) {
          const_5.mat_3[row, y_index_map(i, j, M, W)] <- 1
          const_5.mat_3[row, x_index_map(j, M, W)] <- -1
          row <- row + 1
        }
      }
    }
  }
}

const_5.dir_3 <- rep('<=', times=nrow(const_5.mat_3))
const_5.rhs_3 <- rep(0, times=nrow(const_5.mat_3))

# (4)
# the constraint with the sum over x_imw && x_jmw
const_5.mat_4 <- matrix(rep(0, times=num_teams*(num_teams-1)*num_cities*num_weeks*rhs_length),
                        ncol=rhs_length)
row <- 1
for (i in 1:num_teams) {
  for (j in 1:num_teams) {
    if (i != j) {
      for (M in 1:num_cities) {
        for (W in 1:num_weeks) {
          const_5.mat_4[row, y_index_map(i, j, M, W)] <- 1
        }
      }
      row <- row + 1
    }
  }
}

const_5.dir_4 <- rep('<=', times=nrow(const_5.mat_4)) 
const_5.rhs_4 <- rep(1, times=nrow(const_5.mat_4))


# glue all of the constraint 5 matrices and vectors together
const_5.mat <- rbind(const_5.mat_1, const_5.mat_2, const_5.mat_3, const_5.mat_4)
const_5.dir <- c(const_5.dir_1, const_5.dir_2, const_5.dir_3, const_5.dir_4)
const_5.rhs <- c(const_5.rhs_1, const_5.rhs_2, const_5.rhs_3, const_5.rhs_4)


### constraint 6 ###
# 32 constraints
# for all 1<=i<=32
# sum over weeks of x_iiw = 16
const_6.mat <- matrix(rep(0, times=num_cities*rhs_length),
                      ncol=rhs_length)
row <- 1
for (team in 1:num_teams) {
  for (W in 1:num_weeks) {
    const_6.mat[row, x_index_map(team, team, W)] <- 1
  }
  row <- row + 1
}

const_6.dir <- rep('=', times=num_cities)
const_6.rhs <- rep(num_weeks / 2, times=num_cities)


### constraint 7 ###
# 32 constraints
# for each team 1<=T<=32
# some over months and weeks of x_tmw = 16
const_7.mat <- matrix(rep(0, times=num_teams*rhs_length),
                      ncol=rhs_length)
row <- 1
for (team in 1:num_teams) {
  for (M in 1:num_cities) {
    for (W in 1:num_weeks) {
      const_7.mat[row, x_index_map(team, M, W)] <- 1
    }
  }
  row <- row + 1
}

const_7.dir <- rep('=', times=num_teams)
const_7.rhs <- rep(num_weeks, times=num_teams)


### constraint 8 ###
# 32*16 constraints
# for team, and for each week
# sum over city m of x_tmw = 1
const_8.mat <- matrix(rep(0, times=num_teams*num_weeks*rhs_length),
                      ncol=rhs_length)
row <- 1
for (team in 1:num_teams) {
  for (W in 1:num_weeks) {
    for (M in 1:num_cities) {
      const_8.mat[row, x_index_map(team, M, W)] <- 1
    }
    row <- row + 1
  }
}

const_8.dir <- rep('=', times=num_teams*num_weeks)
const_8.rhs <- rep(1, times=num_teams*num_weeks)


# glue together results from each constraint
nfl.const_mat <- rbind(const_1.mat, const_2.mat, const_3.mat,
                       const_4.mat, const_5.mat, const_6.mat, 
                       const_7.mat, const_8.mat)

nfl.dir <- c(const_1.dir, const_2.dir, const_3.dir, const_4.dir,
             const_5.dir, const_6.dir, const_7.dir, const_8.dir)

nfl.rhs <- c(const_1.rhs, const_2.rhs, const_3.rhs, const_4.rhs,
             const_5.rhs, const_6.rhs, const_7.rhs, const_8.rhs)


# run the lp solver
nfl.lp = lp(direction='min', objective.in=nfl.obj, const.mat=nfl.const_mat,
            const.dir=nfl.dir, const.rhs=nfl.rhs, all.bin=TRUE)

nfl.lp
solution = nfl.lp[['solution']]


# stuff below doesn't work yet (index_to_x function is incorrect)

index_to_x <- function(index) {
  week <- (index - num_cities^2 * num_teams * num_weeks) %% num_weeks
  location <- (index - num_cities^2 * num_teams * num_weeks) %% (num_weeks*num_cities)
  team <- (index - num_cities^2 * num_teams * num_weeks) %% (num_weeks*num_cities*num_teams)
  return(list(week, location, team))
}

organize_results <- function(solution) {
  results_matrix <- matrix(, nrow=num_weeks*num_teams, ncol=3)
  row <- 1
  for (x in solution[(num_cities^2 * num_teams * num_weeks + 1):length(solution)]) {
    if (x != 0) {
      out <- index_to_x(row)
      results_matrix[row, 1] <- out[[1]]
      results_matrix[row, 2] <- out[[2]]
      results_matrix[row, 3] <- out[[3]]
      row <- row + 1
    }
  }
  #results_df <- data.frame(results_matrix, )
  return(results_matrix)
}

# results <- organize_results(solution)


library(lpSolveAPI)
# try printing lp to file
the.lp <- make.lp(nrow(nfl.const_mat), ncol(nfl.const_mat))
set.objfn(the.lp, nfl.obj)
for (i in 1:nrow(nfl.const_mat)) {
  add.constraint(the.lp, nfl.const_mat[i,], nfl.dir[i], nfl.rhs[i])
}
#write.lp(the.lp, 'nfl.lp')
