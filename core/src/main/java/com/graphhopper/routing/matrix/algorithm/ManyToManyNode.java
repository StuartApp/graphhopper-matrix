package com.graphhopper.routing.matrix.algorithm;

import com.graphhopper.routing.matrix.MatrixEntry;
import com.graphhopper.routing.querygraph.QueryRoutingCHGraph;
import com.graphhopper.routing.util.TraversalMode;
import com.graphhopper.storage.*;

public class ManyToManyNode extends AbstractManyToMany {

    private TraversalMode traversalMode = TraversalMode.NODE_BASED;


    public ManyToManyNode(QueryRoutingCHGraph graph,RoutingCHGraph graphNoVirtualNodes){

        super(graph,graphNoVirtualNodes);

        if (graph.hasTurnCosts())
            throw new IllegalStateException("Weightings supporting turn costs cannot be used with node-based traversal mode");
    }

    @Override
    protected int getOrigEdgeId(RoutingCHEdgeIteratorState edge, boolean reverse) {
        return edge.getEdge();
    }


   @Override
    protected double calcWeight(RoutingCHEdgeIteratorState iter, MatrixEntry currEdge,
                                boolean reverse, boolean accumulate){

        if(accumulate){
            return iter.getWeight(reverse) + currEdge.getWeightOfVisitedPath();
        }else{
            return iter.getWeight(reverse);
        }

    }

    @Override
    protected long calcTime(RoutingCHEdgeIteratorState iter, MatrixEntry currEdge,
                            boolean reverse, boolean accumulate){

        if(accumulate){
            return iter.getTime(reverse) + currEdge.time;
        }else{
            return iter.getTime(reverse);
        }
    }

    @Override
    protected double calcDistance(RoutingCHEdgeIteratorState iter, MatrixEntry currEdge, boolean accumulate){


        if(accumulate){
            return iter.getDistance() + currEdge.distance;
        }else{
            return iter.getDistance();
        }
    }

    @Override
    public String getName(){
        return getClass().getSimpleName();
    }
}