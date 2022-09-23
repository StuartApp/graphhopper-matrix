package com.graphhopper.routing.matrix.algorithm;

import com.graphhopper.routing.AlgorithmOptions;
import com.graphhopper.routing.querygraph.QueryGraph;
import com.graphhopper.routing.querygraph.QueryRoutingCHGraph;
import com.graphhopper.routing.util.TraversalMode;
import com.graphhopper.storage.RoutingCHGraph;
import com.graphhopper.util.Helper;
import com.graphhopper.util.Parameters;


public class MatrixRoutingAlgorithmFactory {
    private final QueryRoutingCHGraph routingCHGraph;
    private final RoutingCHGraph graphNoVirtualNodes;

    public MatrixRoutingAlgorithmFactory(RoutingCHGraph routingCHGraph, QueryGraph queryGraph) {
        this.routingCHGraph = new QueryRoutingCHGraph(routingCHGraph, queryGraph);
        this.graphNoVirtualNodes = routingCHGraph;
    }


    public MatrixAlgorithm createAlgo(AlgorithmOptions opts) {
        // CHECK hints instructions and calc_points are false

        String defaultAlgo = Parameters.Algorithms.DIJKSTRA_MANY_TO_MANY;
        String algo = opts.getAlgorithm();
        if (Helper.isEmpty(algo))
            algo = defaultAlgo;
        if (Parameters.Algorithms.DIJKSTRA_MANY_TO_MANY.equals(algo)) {
            if(opts.getTraversalMode() == TraversalMode.NODE_BASED){
               return new ManyToManyNode(routingCHGraph,graphNoVirtualNodes);
            }else{
                return new ManyToManyEdge(routingCHGraph,graphNoVirtualNodes);
            }
        } else {
            throw new IllegalArgumentException("Algorithm " + algo + " not supported for Matrix calculation.");
        }
    }

}
