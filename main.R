source('toolbox.r')
## Code accompanying 
#Song, Von Ahn, Rohr, Saavedra. Towards a probabilistic understanding about the context-dependency of species interactions
#the file toolbox.r contains all the necessary functions
#Please note that the permanence condition is only applicable to 3-species communities.

#Initial community structure A
A <- matrix(c(-1, .4, .2, .2, -1, -.2, .2, -.2, -1), byrow = T, ncol=3)

#Switched interaction giving rise to the new community structure B
B <- A
B[2,1] <- -B[2,1]
B[1,2] <- -B[1,2]

#Probability of transition
compute_prob_transition(A)

#Probability of persistence under constraints
compute_prob_persistence(A, B, exp = 100, constraint = 'with', persistence_type = 'stability')
compute_prob_persistence(A, B, exp = 100, constraint = 'with', persistence_type = 'permanence')

#Probability of persistence without constraints
compute_prob_persistence(A, B, exp = 100, constraint = 'without', persistence_type = 'stability')
compute_prob_persistence(A, B, exp = 100, constraint = 'without', persistence_type = 'permanence')




