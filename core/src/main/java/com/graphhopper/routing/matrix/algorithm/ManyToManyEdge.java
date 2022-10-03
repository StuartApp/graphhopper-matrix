package com.graphhopper.routing.matrix.algorithm;


import com.graphhopper.routing.matrix.MatrixEntry;
import com.graphhopper.routing.matrix.SBIEntry;
import com.graphhopper.routing.querygraph.QueryRoutingCHGraph;
import com.graphhopper.routing.util.TraversalMode;
import com.graphhopper.storage.*;
import com.graphhopper.util.EdgeIteratorState;


public class ManyToManyEdge extends AbstractManyToMany {

    private TraversalMode traversalMode = TraversalMode.EDGE_BASED;

    public ManyToManyEdge(QueryRoutingCHGraph graph, RoutingCHGraph graphNoVirtualNodes){
        super(graph,graphNoVirtualNodes);

        if (!graph.isEdgeBased()) {
            throw new IllegalArgumentException("Edge-based CH algorithms only work with edge-based CH graphs");
        }
    }

    @Override
    protected int getTraversalId(RoutingCHEdgeIteratorState edge, Boolean reverse){

        return traversalMode.createTraversalId(edge.getBaseNode(),edge.getAdjNode(),edge.getEdge(),reverse);
    }

    @Override
    protected int getTraversalId(EdgeIteratorState edge, Boolean reverse){
        return traversalMode.createTraversalId(edge,reverse);
    }


    private double calcWeight(RoutingCHEdgeIteratorState edgeState, Boolean reverse, int prevOrNextEdgeId){

        double edgeWeight = edgeState.getWeight(reverse);
        final int origEdgeId = reverse ? edgeState.getOrigEdgeLast() : edgeState.getOrigEdgeFirst();
        double turnCosts = reverse
                ? graph.getTurnWeight(origEdgeId, edgeState.getBaseNode(), prevOrNextEdgeId)
                : graph.getTurnWeight(prevOrNextEdgeId, edgeState.getBaseNode(), origEdgeId);

        return edgeWeight + turnCosts;

    }

    private double calcWeightNoVirtual(RoutingCHEdgeIteratorState edgeState, Boolean reverse, int prevOrNextEdgeId){

        double edgeWeight = edgeState.getWeight(reverse);
        final int origEdgeId = reverse ? edgeState.getOrigEdgeLast() : edgeState.getOrigEdgeFirst();
        double turnCosts = reverse
                ? graphNoVirtualNodes.getTurnWeight(origEdgeId, edgeState.getBaseNode(), prevOrNextEdgeId)
                : graphNoVirtualNodes.getTurnWeight(prevOrNextEdgeId, edgeState.getBaseNode(), origEdgeId);

        return edgeWeight + turnCosts;
    }

    @Override
    protected double calcWeight(RoutingCHEdgeIteratorState iter, int incomingEdge, boolean reverse){
        return calcWeight(iter,reverse,incomingEdge);
    }

    protected double calcWeightNoVirtual(RoutingCHEdgeIteratorState iter, int incomingEdge, boolean reverse){
        return calcWeightNoVirtual(iter,reverse,incomingEdge);
    }

    private long calcTime(RoutingCHEdgeIteratorState edgeState, Boolean reverse, int prevOrNextEdgeId){

        long time = edgeState.getTime(reverse);
        int origEdgeId;
        if(reverse){
            origEdgeId = edgeState.getOrigEdgeLast();
        }else{
            origEdgeId = edgeState.getOrigEdgeFirst();
        }
        long turnCost;
        if(reverse){
            turnCost = weighting.calcTurnMillis(origEdgeId,edgeState.getBaseNode(),prevOrNextEdgeId);
        }else{
            turnCost = weighting.calcTurnMillis(prevOrNextEdgeId,edgeState.getBaseNode(),origEdgeId);
        }

        return time + turnCost;
    }

    private long calcTimeNoVirtual(RoutingCHEdgeIteratorState edgeState, Boolean reverse, int prevOrNextEdgeId){

        long time = edgeState.getTime(reverse);
        int origEdgeId;
        if(reverse){
            origEdgeId = edgeState.getOrigEdgeLast();
        }else{
            origEdgeId = edgeState.getOrigEdgeFirst();
        }
        long turnCost;
        if(reverse){
            turnCost = weightingNoVirtualNodes.calcTurnMillis(origEdgeId,edgeState.getBaseNode(),prevOrNextEdgeId);
        }else{
            turnCost = weightingNoVirtualNodes.calcTurnMillis(prevOrNextEdgeId,edgeState.getBaseNode(),origEdgeId);
        }

        return time + turnCost;
    }


    @Override
    protected long calcTime(RoutingCHEdgeIteratorState iter, int incomingEdge, boolean reverse){
        return calcTime(iter,reverse,incomingEdge);
    }

    @Override
    protected long calcTimeNoVirtual(RoutingCHEdgeIteratorState iter,int incomingEdge, boolean reverse){
        return calcTimeNoVirtual(iter,reverse,incomingEdge);
    }

    @Override
    protected int getOrigEdgeId(RoutingCHEdgeIteratorState edge, boolean reverse) {
        return reverse ? edge.getOrigEdgeFirst() : edge.getOrigEdgeLast();
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