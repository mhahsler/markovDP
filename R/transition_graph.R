#' Transition Graph
#'
#' Returns the transition model as an \pkg{igraph} object.
#'
#' The transition model of an MDP is a Markov Chain. This function extracts the transition model as
#' an igraph object.
#'
#' @family MDP
#'
#' @importFrom igraph graph_from_adjacency_matrix %>% E E<- V V<- add_layout_ as_data_frame as_tree graph_from_data_frame induced_subgraph norm_coords
#'
#' @param x object of class [MDP].
#' @param action the name or id of an action or a set of actions. By default the transition model for all actions is returned.
#' @param state_col colors used to represent the states.
#' @param simplify_transitions logical; combine parallel transition arcs into a single arc.
#' @param remove_unavailable_actions logical; don't show arrows for unavailable actions.
#'
#' @returns returns the transition model as an igraph object.
#' @examples
#' data("Maze")
#'
#' g <- transition_graph(Maze)
#' g
#'
#' plot_transition_graph(Maze)
#' plot_transition_graph(Maze,
#'   vertex.size = 20,
#'   edge.label.cex = .1, edge.arrow.size = .5, margin = .5
#' )
#'
#' ## Plot using the igraph library
#' library(igraph)
#' plot(g)
#'
#' # plot with a different layout
#' plot(g,
#'   layout = igraph::layout_with_sugiyama,
#'   vertex.size = 20,
#'   edge.label.cex = .6
#' )
#'
#' ## Use visNetwork (if installed)
#' if (require(visNetwork)) {
#'   g_vn <- toVisNetworkData(g)
#'   nodes <- g_vn$nodes
#'   edges <- g_vn$edges
#'
#'   visNetwork(nodes, edges) %>%
#'     visNodes(physics = FALSE) %>%
#'     visEdges(smooth = list(type = "curvedCW", roundness = .6), arrows = "to")
#' }
#' @export
transition_graph <-
  function(x,
           action = NULL,
           state_col = NULL,
           simplify_transitions = TRUE,
           remove_unavailable_actions = TRUE) {
    state_col <-
      colors_discrete(length(x$states), state_col)

    m <-
      transition_matrix(
        x,
        action = NULL,
        sparse = FALSE
      )

    if (is.null(action)) {
      action <- x$actions
    }

    gs <- sapply(
      action,
      FUN = function(a) {
        g <-
          graph_from_adjacency_matrix(m[[a]], mode = "directed", weighted = TRUE)
        E(g)$label <- a
        df <- as_data_frame(g)

        # remove unavailable actions
        if (remove_unavailable_actions) {
          available <- sapply(seq_len(nrow(df)), FUN = function(i) {
            all(reward_matrix(x, action = df$label[i], start.state = df$from[i], 
                              end.state = df$to[i]) != -Inf)
          })
          df <- df[available, , drop = FALSE]
        }
        df
      }, simplify = FALSE
    )

    g <- graph_from_data_frame(do.call(rbind, gs))
    # make sure the vertices are in the same order as in the description
    g <- igraph::permute(g, match(V(g)$name, x$states))


    E(g)$label <- paste0(E(g)$label, ifelse(E(g)$weight != 1, paste0(" (", round(E(g)$weight, 2), ")"), ""))

    if (simplify_transitions) {
      g <- igraph::simplify(
        g,
        edge.attr.comb = list(
          label = function(x) {
            paste(x, collapse = "/\n")
          },
          "ignore"
        ),
        remove.loops = FALSE
      )
    }

    if (!any(is.na(state_col))) {
      V(g)$color <- state_col
    }

    g
  }

#' @rdname transition_graph
#' @param main a main title for the plot.
#' @param ... further arguments are passed on to `igraph::plot.igraph()`.
#' @export
plot_transition_graph <- function(x,
                                  action = NULL,
                                  state_col = NULL,
                                  simplify_transitions = TRUE,
                                  main = NULL,
                                  ...) {
  g <- transition_graph(
    x,
    action = action,
    state_col = state_col,
    simplify_transitions = simplify_transitions
  )

  if (is.null(main)) {
    main <- paste("Transition Graph:", x$name)
  }

  plot(
    g,
    edge.curved = curve_multiple_directed(g, .8),
    edge.loop.angle = -pi / 4, main = main,
    ...
  )
}

#' @rdname transition_graph
#' @param graph The input graph.
#' @param start The curvature at the two extreme edges.
#' @export
curve_multiple_directed <- function(graph, start = 0.3) {
  el <- igraph::as_edgelist(graph, names = FALSE)
  o <- apply(el, 1, order)[1, ]
  el <-
    apply(
      el,
      1,
      FUN = function(x) {
        paste(sort(x), collapse = ":")
      }
    )
  cu <- stats::ave(
    rep(NA, length(el)),
    el,
    FUN = function(x) {
      if (length(x) == 1) {
        return(0)
      } else {
        return(seq(-start, start, length = length(x)))
      }
    }
  )

  cu[o == 2] <- cu[o == 2] * -1
  cu
}
