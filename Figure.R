source('toolbox.r')
# Code to reproduce the results of Song, Von Ahn, Rohr, Saavedra. Towards a probabilistic understanding about the context-dependency of species interactions
# Generate all possible network topologies for 3-species communities
# For each topology, we generate a weighted interaction matrix drawing values from a half normal distribution.
# Then, we assume 1-point transitions and 1 transition at a time
# For each posible switch we calculate its swiching probability
# This switching probability is given by the product of the transition probability and the persistence probability (see Methods)
# The transition probability (in this exercise) is given by the inverse of the number of possible transitions of the community
# The persistence probability is given by the re-normalized size of the persistence domain (see Methods)
# Because we are generating random interaction matrices, the switching probability is estimated as the mean over 100 realizations.
# These probabilities are calculated with and without additional constraints on the parameter space (see Methods)
AllType_Inte <- Generate_AllType_Inte()

Nrand <- 100 #number of randomizations

## This is the main part of the code, where we calculate the probabilities for all possible topolgies and 1-point transitions with and without constraints
pb <- progress_estimated(Nrand) #progress bar
#sample the switching probability via simulation
result <- 1:Nrand %>% 
  map_dfr(~sample_switch_probability(persistence_type = 'permanence'))

#Here we plot the transition probabilities
result %>% 
  group_by(from, to, type_of_transition) %>%
  summarise_at(c("prob_transition", "prob_persistence_with", "prob_persistence_without"), mean) %>%
  ungroup() %>%
  select(-from, -to) %>% 
  gather(constraint, prob_persistence, -type_of_transition,-prob_transition) %>% 
  filter(type_of_transition %in% c('2->3', '3->2')) %>% # We filter results from antagonistic to mutualistic and vice versa. 1: competition (-,-);  2:antagonistic (+,-); and 3:mutualism (+,+). If the user needs from competition to antagonism, then simply write '1->2'. Note that we cannot write '1->3' from competition to mutualism as it will require a 2-point transition (which can be done with extra code made by the user).  
  ggplot(aes(type_of_transition, prob_transition, color = type_of_transition)) +
  geom_boxplot()+
  geom_jitter()+
  facet_wrap(~constraint)+
  labs(
    x = 'Interaction switch',
    y = 'Transition probability',
    title = 'Transition probability'
  )+
  scale_color_manual(values=c("#39495C", "#FDC03B")) +
  scale_x_discrete(breaks=c("2->3","3->2"),
                   labels=c(expression(bold("A") %->% bold("M")), expression(bold("M") %->% bold("A")))
  )+
  theme_bw()+
  theme(
    strip.background = element_blank(),
    strip.text.x = element_blank(),
    aspect.ratio = 1,
    text = element_text(size=24),
    legend.position = "none",
    panel.spacing = unit(2, "lines")
  )

##Here we plot the persistence probability
result %>% 
  group_by(from, to, type_of_transition) %>%
  summarise_at(c("prob_transition", "prob_persistence_with", "prob_persistence_without"), mean) %>%
  ungroup() %>%
  select(-from, -to) %>% 
  gather(constraint, prob_persistence, -type_of_transition,-prob_transition) %>% 
  filter(type_of_transition %in% c('2->3', '3->2')) %>%
  ggplot(aes(type_of_transition, prob_persistence, color = type_of_transition)) +
  geom_boxplot()+
  geom_jitter()+
  facet_wrap(~constraint)+
  labs(
    x = 'Interaction switch',
    y = 'Persistence probability',
    title = 'Persistence probability'
  )+
  scale_color_manual(values=c("#39495C", "#FDC03B")) +
  scale_x_discrete(breaks=c("2->3","3->2"),
                   labels=c(expression(bold("A") %->% bold("M")), expression(bold("M") %->% bold("A")))
  )+
  theme_bw()+
  theme(
    strip.background = element_blank(),
    strip.text.x = element_blank(),
    aspect.ratio = 1,
    text = element_text(size=24),
    legend.position = "none",
    panel.spacing = unit(2, "lines")
  )

##Here we plot the switching probability
result %>% 
  group_by(from, to, type_of_transition) %>%
  summarise_at(c("prob_transition", "prob_persistence_with", "prob_persistence_without"), mean) %>%
  ungroup() %>%
  select(-from, -to) %>% 
  gather(constraint, prob_persistence, -type_of_transition,-prob_transition) %>% 
  filter(type_of_transition %in% c('2->3', '3->2')) %>%
  ggplot(aes(type_of_transition, prob_persistence * prob_transition, color = type_of_transition)) +
  geom_boxplot()+
  geom_jitter()+
  facet_wrap(~constraint)+
  labs(
    x = 'Interaction switch',
    y = 'Switching probability',
    title = 'Switching probability'
  )+
  scale_color_manual(values=c("#39495C", "#FDC03B")) +
  scale_x_discrete(breaks=c("2->3","3->2"),
                   labels=c(expression(bold("A") %->% bold("M")), expression(bold("M") %->% bold("A")))
  )+
  theme_bw()+
  theme(
    strip.background = element_blank(),
    strip.text.x = element_blank(),
    aspect.ratio = 1,
    text = element_text(size=24),
    legend.position = "none",
    panel.spacing = unit(2, "lines")
  )
