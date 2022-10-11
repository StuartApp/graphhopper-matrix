package com.graphhopper.routing.matrix.algorithm;


import com.graphhopper.routing.matrix.MatrixEntry;
import com.graphhopper.routing.querygraph.QueryRoutingCHGraph;
import com.graphhopper.routing.util.TraversalMode;
import com.graphhopper.routing.weighting.Weighting;
import com.graphhopper.storage.*;
import com.graphhopper.util.GHUtility;


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


        //int origEdgeId = getOrigEdgeId(edge,reverse);
        //int baseNode = getOtherNode(origEdgeId,edge.getAdjNode());

        //return GHUtility.createEdgeKey(baseNode, edge.getAdjNode(), origEdgeId, reverse);

        return GHUtility.createEdgeKey(edge.getBaseNode(), edge.getAdjNode(), edge.getEdge(), reverse);
    }

    @Override
    protected int getOrigEdgeId(RoutingCHEdgeIteratorState edge, boolean reverse) {
        return reverse ? edge.getOrigEdgeFirst() : edge.getOrigEdgeLast();
    }


    private double calcWeight(RoutingCHEdgeIteratorState edgeState, Boolean reverse, int prevOrNextEdgeId){

        double edgeWeight = edgeState.getWeight(reverse);
        int origEdgeId = reverse ? edgeState.getOrigEdgeLast() : edgeState.getOrigEdgeFirst();

        double turnCosts = reverse
                    ? graph.getTurnWeight(origEdgeId, edgeState.getBaseNode(), prevOrNextEdgeId)
                    : graph.getTurnWeight(prevOrNextEdgeId, edgeState.getBaseNode(), origEdgeId);

        return edgeWeight + turnCosts;
    }

    @Override
    protected double calcWeight(RoutingCHEdgeIteratorState iter, MatrixEntry currEdge,
                                boolean reverse, boolean accumulate){

        double w = calcWeight(iter,reverse,currEdge.originalEdge);
        if(accumulate){
            return  w + currEdge.getWeightOfVisitedPath();
        }else{
            return w;
        }
    }

    private long calcTime(RoutingCHEdgeIteratorState edgeState, Boolean reverse,
                          int prevOrNextEdgeId){

        long time = edgeState.getTime(reverse);
        int origEdgeId = getOrigEdgeId(edgeState,reverse);

        long turnCost;

        if(reverse){
            turnCost = weighting.calcTurnMillis(origEdgeId,edgeState.getBaseNode(),prevOrNextEdgeId);
        }else{
            turnCost = weighting.calcTurnMillis(prevOrNextEdgeId,edgeState.getBaseNode(),origEdgeId);
        }

        return time + turnCost;
    }

    @Override
    protected long calcTime(RoutingCHEdgeIteratorState iter, MatrixEntry currEdge,
                            boolean reverse, boolean accumulate){

        long time = calcTime(iter,reverse,currEdge.originalEdge);
        if(accumulate){
             return  time + currEdge.time;
         }else{
             return time;
         }

    }

    protected int getIncomingEdge(MatrixEntry entry) {
        return entry.edge;
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