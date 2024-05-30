#ifndef List_MODEL_H
#define List_MODEL_H

#include <Rcpp.h>
#include <numeric>

#include "dgCMatrix.h"
#include "math.h"

using namespace Rcpp;

// C++ interface to access elements of a List/List model

// NOTE: Episode in time-dependent Lists are currently unsupported.
// NOTE: All indices are 0-based.

typedef List List;
typedef List List;

// Access model information
bool is_solved(const List& model);
bool is_converged(const List& model);

// More accessors
NumericVector start_vector(const List& model);
CharacterVector get_states(const List& model);
LogicalVector absorbing_states(const List& model);
CharacterVector get_obs(const List& model);
CharacterVector get_actions(const List& model);
double get_discount(const List& model);

// NA is inf, we return the sum of episode horizons for now
int get_horizon(const List& model);

// Transitions
// Can be a dense matrix or a dgCMatrix
// Available functions are: x_matrix returns a dense matrix, x_prob returns double, and x_row returns a vector
NumericMatrix transition_matrix(const List& model, int action, int episode = -1);
double transition_prob(const List& model, int action, int start_state, 
                              int end_state, int episode = -1);
NumericVector transition_row(const List& model, int action, int start_state, 
                                    int episode = -1);


// Reward
// List has no observation!
// Note: R_index does not apply to episode!!!
// Note: Lists don't use observations (observations) use observation = 0!
NumericMatrix reward_matrix_MDP(const List& model, int action, 
                                       int start_state, int episode = -1);
double reward_val_MDP(const List& model, int action, 
                         int start_state, int end_state,
                         int episode = -1, bool R_index = FALSE);


// returns the List policy as a vector. Index is the state index and the value is the action index.
IntegerVector get_policy_MDP(const List& model);


// get pg and alpha epochs (in case of non converged policies)
// epochs start with 0
int get_pg_index_cpp(const List& model, int epoch);
NumericMatrix get_alpha(const List& model, int epoch = 0);
DataFrame get_pg(const List& model, int epoch = 0);

// terminal value
double terminal_val(const List& model, int state);

// Observations
NumericMatrix observation_matrix(const List& model, int action,
                                        int episode = -1);
double observation_prob(const List& model, int action, int end_state, 
                               int observation, int episode = -1);
NumericVector observation_row(const List& model, int action, int end_state, 
                                     int episode = -1);

// Reward
// Note: R_index does not apply to episode!!!
// TODO add support for episodes
// Available are reward_matrix and reward_val
NumericMatrix reward_matrix(const List& model, int action, int start_state, 
                                   int episode = -1);

// Note: R_index does not apply to episode!!!
// Note: Lists don't use observations (observations) use observation = 0!
double reward_val(const List& model, int action, 
                         int start_state, int end_state, int observation = 0,
                         int episode = -1, bool R_index = FALSE);


// C++ interface
DataFrame reward_cpp(const NumericMatrix& belief, const NumericMatrix& alpha);
NumericVector update_belief_cpp(const List& model, const NumericVector& belief,
                                int action, int observation, int digits);


#endif