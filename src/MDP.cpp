#include "MDP.h"

#include <Rcpp.h>
using namespace Rcpp;

// C++ interface to access elements of a List/List model

// NOTE: Episode in time-dependent Lists are currently unsupported.
// NOTE: All indices are 0-based.

// Call R functions
Environment pkg = Environment::namespace_env("markovDP");
Function R_start_vector = pkg["start_vector"];
Function R_absorbing_states = pkg["absorbing_states"];

// R function calls are very slow
Function R_transition_matrix = pkg["transition_matrix"];
Function R_reward_matrix = pkg["reward_matrix"];

// Access model information
bool is_solved(const List& model) { 
  return model.containsElementNamed("solution");
}

bool is_converged(const List& model) { 
  return as<LogicalVector>(as<List>(model["solution"])["converged"])[0];
}

// More accessors
NumericVector start_vector(const List& model) {
  return as<NumericVector>(R_start_vector(model));
}  

CharacterVector get_states(const List& model) {
  return as<CharacterVector>(model["states"]);
}  

// TODO: reimplement in C++
IntegerVector absorbing_states(const List& model) {
  //return R_absorbing_states(model);
  return R_absorbing_states(model, _["sparse"] = CharacterVector::create("index"));
}

CharacterVector get_actions(const List& model) {
  return as<CharacterVector>(model["actions"]);
}  

// TODO: Add available_actions()

double get_discount(const List& model) {
  return model["discount"];
} 

// NA is inf, we return the sum of episode horizons for now
int get_horizon(const List& model) {
  NumericVector h = model["horizon"];
  if (!is_finite(h)[0]) 
    return R_NaInt;
  
  return (int) sum(h);
}  

// Transitions
// Can be a dense matrix or a dgCMatrix
// Available functions are: x_matrix returns a dense matrix, x_prob returns double, and x_row returns a vector
NumericMatrix transition_matrix(const List& model, int action, int episode) {
  RObject acts;
  if (episode >= 0)
    acts = as<List>(model["transition_prob"])[episode];
  else
    acts = model["transition_prob"];

  // function
  if (is<Function>(acts)) {
    return as<NumericMatrix>(R_transition_matrix(model, action + 1, 
                                                 _["sparse"] = LogicalVector::create(0)));
  }
 
  // it is a list
  acts = as<List>(acts)[action];
 
  // dense matrix
  if (is<NumericMatrix>(acts)) 
    return as<NumericMatrix>(acts); 
 
  // dgCMatrix
  if (is<S4>(acts))
    return dgCMatrix(as<S4>(acts)).dense();
  
  // uniform/identity
  if (is<CharacterVector>(acts)) {
    int n_states = get_states(model).size();
    if (as<CharacterVector>(acts)[0] == "uniform") {
      NumericVector m(n_states * n_states, 1.0 / n_states); 
      m.attr("dim") = IntegerVector::create(n_states, n_states); 
      return as<NumericMatrix>(m);
    }
    
    if (as<CharacterVector>(acts)[0] == "identity") {
      NumericMatrix m = NumericMatrix::diag(n_states, 1.0);
      return m;
    }
    
    stop("Unknown matrix specifier! Only 'identity' and 'uniform' are allowed.");
  }
  
  stop("transition_matrix: model needs to be normalized with normalize_MDP().");
}

double transition_prob(const List& model, int action, int start_state, 
                       int end_state, int episode) {
  RObject acts;
  if (episode >= 0)
    acts = as<List>(model["transition_prob"])[episode];
  else
    acts = model["transition_prob"];
  
  // function
  if (is<Function>(acts)) {
    return as<double>(R_transition_matrix(model, action + 1, start_state + 1, end_state + 1));
  }
  
  // it is a list
  acts = as<List>(acts)[action];
  
  // dense matrix
  if (is<NumericMatrix>(acts)) 
    return as<NumericMatrix>(acts)(start_state, end_state); 
  
  // dgCMatrix
  if (is<S4>(acts))
    return dgCMatrix(as<S4>(acts)).at(start_state, end_state);
  
  // uniform/identity
  if (is<CharacterVector>(acts)) {
    int n_states = get_states(model).size();
    if (as<CharacterVector>(acts)[0] == "uniform")
      return (1.0 / n_states);
    
    if (as<CharacterVector>(acts)[0] == "identity") {
      if (start_state == end_state)
        return 1.0;
      else
        return 0.0;
    }
    
    stop("Unknown matrix specifier! Only 'identity' and 'uniform' are allowed.");
  }
  
  stop("transition_prob: model needs to be normalized with normalize_MDP().");
}

NumericVector transition_row(const List& model, int action, int start_state, 
                             int episode) {
  RObject acts;
  if (episode >= 0)
    acts = as<List>(model["transition_prob"])[episode];
  else
    acts = model["transition_prob"];
  
  // function
  if (is<Function>(acts)) {      
    return as<NumericVector>(R_transition_matrix(model, action + 1, start_state + 1,
                                                 _["sparse"] = LogicalVector::create(0)));
  }
  
  // it is a list
  acts = as<List>(acts)[action];
  
  // dense matrix
  if (is<NumericMatrix>(acts)) 
    return as<NumericMatrix>(acts).row(start_state); 
  
  // dgCMatrix
  if (is<S4>(acts))
    return dgCMatrix(as<S4>(acts)).row(start_state);
  
  // uniform/identity
  if (is<CharacterVector>(acts)) {
    int n_states = get_states(model).size();
    if (as<CharacterVector>(acts)[0] == "uniform") {
      NumericVector v(n_states, 1.0 / n_states);
      return v;
    }
    
    if (as<CharacterVector>(acts)[0] == "identity") {
      NumericVector v(n_states, 0.0);
      v[start_state] = 1.0;
      return v;
    }
    
    stop("Unknown matrix specifier! Only 'identity' and 'uniform' are allowed.");
  }
  
  stop("transition_row: model needs to be normalized with normalize_MDP().");
}


// Reward
// List has no observation!
NumericMatrix reward_matrix_MDP(const List& model, int action) {
  RObject reward = model["reward"];
  
  if (is<Function>(reward)) {
    return as<NumericMatrix>(R_reward_matrix(model, action + 1,
                                             _["sparse"] = LogicalVector::create(0)));
  }
  
  if (is<DataFrame>(reward)) {
    DataFrame df = as<DataFrame>(reward);
    IntegerVector actions = df[0], start_states = df[1], end_states = df[2];
    NumericVector values = df[3]; 
    
    NumericMatrix rew(get_states(model).size(), get_states(model).size());
    
    for (auto i = 0; i < df.nrows(); ++i) {
      if(!(IntegerVector::is_na(actions[i]) || actions[i] == action))
        continue;
        
      if(!IntegerVector::is_na(start_states[i]) && 
         !IntegerVector::is_na(end_states[i])) {
        rew(start_states[i], end_states[i]) = values[i];
      } else if (!IntegerVector::is_na(start_states[i])) {
        rew(_ , end_states[i]) = NumericVector(rew.rows(), values[i]);
      } else {
        rew(start_states[i], _) = NumericVector(rew.cols(), values[i]);
      }
    }
        
    return rew;
  }
  
  // dgCMatrix
  if (is<S4>(reward))
    return dgCMatrix(as<S4>(reward)).dense();
  
  // it is a matrix
  return as<NumericMatrix>(as<List>(reward)[action]);
}


// Note: R_index does not apply to episode!!!
// Note: Lists don't use observations (observations) use observation = 0!
double reward_val_MDP(const List& model, int action, 
                         int start_state, int end_state,
                         int episode, bool R_index) {
  RObject reward = model["reward"];
  RObject acts;
  
  if (episode >= 0)
    reward = as<List>(reward)[episode];
  
  if (is<Function>(reward)) {
    return as<double>(R_reward_matrix(model, action + 1, start_state + 1,
                                      end_state + 1));
  }
  
  if (is<DataFrame>(reward)) {
    // factors in the data.frame are 1-based!!!
    if (!R_index) {
      action++; start_state++; end_state++;
    }
    
    DataFrame df = as<DataFrame>(reward);
    // find the best matching entry
    IntegerVector actions = df[0], start_states = df[1], end_states = df[2];
    NumericVector values = df[3]; 
    
    for (auto i = df.nrows()-1; i >= 0; --i) {
      if(
          (IntegerVector::is_na(actions[i]) || actions[i] == action) && 
          (IntegerVector::is_na(start_states[i]) || start_states[i] == start_state) &&
          (IntegerVector::is_na(end_states[i]) || end_states[i] == end_state)
        )
        return values[i];
        
    }
    return 0.0;
  }
 
 // it is a list
 acts = as<List>(reward)[action];
 
 // dense matrix
 if (is<NumericMatrix>(acts)) 
   return as<NumericMatrix>(acts)(start_state, end_state); 
 
 // dgCMatrix
 if (is<S4>(acts))
   return dgCMatrix(as<S4>(acts)).at(start_state, end_state);
 
 stop("reward_val: model needs to be normalized with normalize_MDP().");

}  


// returns a vector of actions. Index is the state index and the value is the action index.
IntegerVector get_policy_MDP(const List& model) {
  if (!is_solved(model))
    stop("Unsolved List model. No policy available");
  
  return as<IntegerVector>(as<List>(as<List>(as<List>(model["solution"])["policy"])[0])["action"]) - 1;
}


