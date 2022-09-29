package com.graphhopper.routing.matrix.algorithm;

import com.graphhopper.routing.matrix.MatrixEntry;
import com.graphhopper.routing.matrix.SBIEntry;
import com.graphhopper.routing.querygraph.QueryRoutingCHGraph;
import com.graphhopper.routing.util.TraversalMode;
import com.graphhopper.storage.*;
import com.graphhopper.util.EdgeIteratorState;

public class ManyToManyNode extends AbstractManyToMany {

    private TraversalMode traversalMode = TraversalMode.NODE_BASED;


    public ManyToManyNode(QueryRoutingCHGraph graph,RoutingCHGraph graphNoVirtualNodes){

        super(graph,graphNoVirtualNodes);

        if (graph.hasTurnCosts())
            throw new IllegalStateException("Weightings supporting turn costs cannot be used with node-based traversal mode");
    }

    @Override
    protected int getTraversalId(RoutingCHEdgeIteratorState state, int origEdgeId,Boolean reverse){
        return traversalMode.createTraversalId(state.getBaseNode(),state.getAdjNode(),state.getEdge(),reverse);
    }

    @Override
    protected int getTraversalId(EdgeIteratorState state, Boolean reverse){
        return traversalMode.createTraversalId(state,reverse);
    }


    @Override
    protected double calcWeight(RoutingCHEdgeIteratorState iter, int incomingEdge, boolean reverse){
        return iter.getWeight(reverse);
    }

    protected double calcWeightNoVirtual(RoutingCHEdgeIteratorState iter, int incomingEdge, boolean reverse){
        return iter.getWeight(reverse);
    }


    @Override
    protected int getIncomingEdge(SBIEntry entry) {
        return super.getIncomingEdge(entry);
    }


    @Override
    protected long calcTime(RoutingCHEdgeIteratorState iter, int incomingEdge, boolean reverse){
        return iter.getTime(reverse);
    }

    @Override
    protected long calcTimeNoVirtual(RoutingCHEdgeIteratorState iter, int incomingEdge, boolean reverse){
        return iter.getTime(reverse);
    }

    @Override
    protected double calcDistance(RoutingCHEdgeIteratorState iter){
        return iter.getDistance();
    }

    @Override
    public String getName(){
        return getClass().getSimpleName();
    }
}