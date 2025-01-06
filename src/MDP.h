#ifndef List_MODEL_H
#define List_MODEL_H

#include <Rcpp.h>
#include <numeric>

#include "dgRMatrix.h"
#include "math.h"

using namespace Rcpp;

// C++ interface to access elements of a List/List model

// NOTE: All indices are 0-based.

typedef List List;
typedef List List;

// Access model information
bool is_solved(const List& model);
bool is_converged(const List& model);

// More accessors
NumericVector start_vector(const List& model);
CharacterVector get_states(const List& model);
IntegerVector absorbing_states(const List& model);
CharacterVector get_obs(const List& model);
CharacterVector get_actions(const List& model);
double get_discount(const List& model);

// NA is inf, we return the sum of episode horizons for now
int get_horizon(const List& model);

// Transitions
// Can be a dense matrix or a dgRMatrix
// Available functions are: x_matrix returns a dense matrix, x_prob returns double, and x_row returns a vector
NumericMatrix transition_matrix(const List& model, int action, int episode = -1);
double transition_prob(const List& model, int action, int start_state, 
                              int end_state, int episode = -1);
NumericVector transition_row(const List& model, int action, int start_state, 
                                    int episode = -1);


// Reward
NumericMatrix reward_matrix_MDP(const List& model, int action);

double reward_val_MDP(const List& model, int action, 
                         int start_state, int end_state);


// returns the policy as a vector. Index is the state index and the value is the action index.
IntegerVector get_policy_MDP(const List& model);

#endif