package com.graphhopper.routing.matrix.algorithm;

import com.carrotsearch.hppc.*;
import com.graphhopper.coll.GHIntObjectHashMap;
import com.graphhopper.routing.matrix.*;
import com.graphhopper.routing.querygraph.QueryRoutingCHGraph;
import com.graphhopper.routing.weighting.Weighting;
import com.graphhopper.storage.*;
import com.graphhopper.storage.index.Snap;

import java.util.*;

/*
 * Implements Simultaneous Bucket Initialization Many-to-Many algorithm (Theo Wieland)
 * https://i11www.iti.kit.edu/_media/teaching/theses/ba_wieland22.pdf
 */

public abstract class AbstractManyToMany implements MatrixAlgorithm {


    protected RoutingCHGraph graph;
    protected RoutingCHGraph graphNoVirtualNodes;

    protected Weighting weighting;
    protected Weighting noVirtualWeighting;

    protected RoutingCHEdgeExplorer inEdgeExplorer;
    protected RoutingCHEdgeExplorer inEdgeExplorerNoVirtual;
    protected RoutingCHEdgeExplorer outEdgeExplorer;
    protected RoutingCHEdgeExplorer outEdgeExplorerNoVirtual;

    protected CHEdgeFilter sbiLevelEdgeFilter;

    protected boolean alreadyRun = false;

    protected int size;
    protected int maxVisitedNodes = Integer.MAX_VALUE;
    protected int visitedNodes = 0;
    protected int maxNodes;
    protected PriorityQueue<RankedNode> heap;
    protected SBIAlgorithm sbiAlgo;
    IntSet traversed;


    public AbstractManyToMany(QueryRoutingCHGraph graph, RoutingCHGraph graphNoVirtualNodes) {

        this.graph = graph;
        this.weighting = graph.getWrappedWeighting();
        this.inEdgeExplorer = graph.createInEdgeExplorer();
        this.outEdgeExplorer = graph.createOutEdgeExplorer();
        this.maxNodes = graph.getBaseGraph().getBaseGraph().getNodes();

        this.graphNoVirtualNodes = graphNoVirtualNodes;
        this.noVirtualWeighting = graphNoVirtualNodes.getWeighting();
        this.inEdgeExplorerNoVirtual = graphNoVirtualNodes.createInEdgeExplorer();
        this.outEdgeExplorerNoVirtual = graphNoVirtualNodes.createOutEdgeExplorer();

        this.traversed = new IntHashSet();


        this.sbiLevelEdgeFilter = new CHEdgeFilter() {

            @Override
            public boolean accept(RoutingCHEdgeIteratorState edgeState) {

                int base = edgeState.getBaseNode();
                int adj = edgeState.getAdjNode();

                // always accept virtual edges, see #288
                if (base >= maxNodes || adj >= maxNodes) return true;

                // minor performance improvement: shortcuts in wrong direction are disconnected, so no need to exclude them
                if (edgeState.isShortcut()) return true;

                return graph.getLevel(base) <= graph.getLevel(adj);
            }
        };

        this.size = Math.min(Math.max(200, graph.getNodes() / 10), 150_000);
        this.heap = new PriorityQueue<>(size);
    }

    @Override
    public DistanceMatrix calcMatrix(List<Snap> sources, List<Snap> targets) {

        checkAlreadyRun();

        DistanceMatrix matrix = new DistanceMatrix(sources.size(), targets.size());
        this.sbiAlgo = new SBIAlgorithm(graph,heap,matrix);
        IntObjectMap<IntArrayList> targetIdxsNodes = new GHIntObjectHashMap<>(targets.size());

        //Backward
        int idxTarget = 0;
        while (idxTarget < targets.size()) {

            int targetClosestNode = targets.get(idxTarget).getClosestNode();
            sbiAlgo.addTarget(targetClosestNode,idxTarget);

            //Find Buckets
            findInitialNodesBackward(targets.get(idxTarget), idxTarget);
            idxTarget++;
        }

        backward();

        //Reset collections for forward
        this.heap.clear();
        this.traversed.clear();

        //Forward
        int idxSource = 0;
        IntIntMap processedSources = new IntIntHashMap();
        List<IntIntPair> cloneResults = new ArrayList<>();
        while (idxSource < sources.size()) {
            //If we have n sources pointing to the same closestNode, we don't want to calculate both, we will copy the result of
            //one to the other
            int closestNode = sources.get(idxSource).getClosestNode();

            sbiAlgo.addSource(closestNode,idxSource);

            if (processedSources.containsKey(closestNode)) {
                cloneResults.add(new IntIntPair(processedSources.get(closestNode), idxSource));
            } else {
                findInitialNodesForward(sources.get(idxSource), idxSource);
                processedSources.put(closestNode, idxSource);
            }

            idxSource++;
        }

        forward();

        cloneResults.stream().forEach(pair -> matrix.copyResult(pair.source, pair.target));

        return matrix;
    }

    protected void checkAlreadyRun() {
        if (alreadyRun) throw new IllegalStateException("Create a new instance per call");
        alreadyRun = true;
    }


    protected abstract int getOrigEdgeId(RoutingCHEdgeIteratorState edge, boolean reverse);


    protected int getOrigEdgeIdV2(RoutingCHEdgeIteratorState edge, boolean reverse){
        return (reverse) ? edge.getOrigEdgeLast() : edge.getOrigEdgeFirst();
    }

    protected abstract double calcWeight(RoutingCHEdgeIteratorState iter, MatrixEntry currEdge,
                                         boolean reverse, boolean accumulate);

    protected abstract long calcTime(RoutingCHEdgeIteratorState iter, MatrixEntry currEdge,
                                     boolean reverse, boolean accumulate);

    protected abstract double calcDistance(RoutingCHEdgeIteratorState iter, MatrixEntry currEdge,
                                           boolean accumulate);
    @Override
    public int getVisitedNodes() {
        return visitedNodes;
    }

    @Override
    public void setMaxVisitedNodes(int numberOfNodes) {
        this.maxVisitedNodes = numberOfNodes;
    }

    private boolean isVirtual(int node) {
        return node >= maxNodes;
    }

    private int getLevel(int node) {
        if (isVirtual(node)) {
            return 0;
        } else {
            return graph.getLevel(node);
        }
    }

    private void findInitialNodesBackward(Snap snap, int idx) {

        int closestNode = snap.getClosestNode();
        int closestLevel = getLevel(closestNode);
        IntSet processed = new IntHashSet();

        boolean closestIsVirtual = isVirtual(closestNode);

        if(!closestIsVirtual && !this.traversed.contains(closestNode)){

            boolean noAccessibleNodes = true;

            MatrixEntry current = new MatrixEntry(closestNode, 0, 0, 0, closestLevel);
            RoutingCHEdgeIterator downIterator = inEdgeExplorer.setBaseNode(closestNode);

            while (downIterator.next()) {

                if(!sbiLevelEdgeFilter.accept(downIterator)){
                    continue;
                }

                int adjNode = downIterator.getAdjNode();
                double weight = calcWeight(downIterator, current, true, true);
                boolean isVirtualAdj = isVirtual(adjNode);
                int levelAdj = getLevel(adjNode);
                int origEdgeId = getOrigEdgeIdV2(downIterator, true);

                if (weight == Double.POSITIVE_INFINITY || isVirtualAdj) {
                    continue;
                }

                double distance = calcDistance(downIterator, current, true);
                long time = calcTime(downIterator, current, true, true);

                sbiAlgo.addInitialOutVertex(adjNode,closestNode,idx,weight,time,distance,true,
                        origEdgeId,downIterator.getOrigEdgeFirst(), downIterator.getOrigEdgeLast());
                heap.add(new RankedNode(adjNode,levelAdj,false));
                traversed.add(adjNode);
                noAccessibleNodes = false;
            }

            heap.add(new RankedNode(closestNode,closestLevel,noAccessibleNodes));
            this.traversed.add(closestNode);

        }else{
            MatrixEntry initial = new MatrixEntry(closestNode, 0, 0, 0, getLevel(closestNode));
            Deque<MatrixEntry> queue = new ArrayDeque<>();
            queue.add(initial);

            while (!queue.isEmpty()) {

                MatrixEntry current = queue.poll();
                int baseNode = current.adjNode;

                RoutingCHEdgeIterator downIterator = inEdgeExplorer.setBaseNode(baseNode);
                while (downIterator.next()) {

                    int adjNode = downIterator.getAdjNode();

                    int origEdgeId = getOrigEdgeIdV2(downIterator, true);
                    double weight = calcWeight(downIterator, current, true, true);
                    boolean isVirtualAdj = isVirtual(adjNode);

                    if (weight == Double.POSITIVE_INFINITY) {
                        continue;
                    }


                    double distance = calcDistance(downIterator, current, true);
                    long time = calcTime(downIterator, current, true, true);

                    MatrixEntry entry = new MatrixEntry(downIterator.getEdge(), origEdgeId, adjNode, current.adjNode, getLevel(adjNode), weight, time, distance);

                    if (isVirtualAdj && adjNode != closestNode && !processed.contains(adjNode)) {
                        processed.add(adjNode);
                        queue.add(entry);
                    } else if (adjNode != closestNode && !processed.contains(adjNode)) {
                        if(!this.traversed.contains(adjNode)){
                            heap.add(new RankedNode(adjNode,getLevel(adjNode),false));
                            this.traversed.add(adjNode);
                        }
                    }

                    sbiAlgo.addInitialOutVertex(adjNode,closestNode,idx,weight,time,distance,true, downIterator.getOrigEdge(),
                            downIterator.getOrigEdgeFirst(),downIterator.getOrigEdgeLast());
                }
            }
        }


    }

    private void findInitialNodesForward(Snap snap, int idx) {

        int closestNode = snap.getClosestNode();
        IntSet processed = new IntHashSet();
        boolean closestIsVirtual = isVirtual(closestNode);
        int closestLevel = getLevel(closestNode);

        if(!closestIsVirtual && !this.traversed.contains(closestNode)){

            MatrixEntry current = new MatrixEntry(closestNode, 0, 0, 0, closestLevel);
            RoutingCHEdgeIterator downIterator = outEdgeExplorer.setBaseNode(closestNode);
            boolean noAccessibleNodes = true;
            while (downIterator.next()) {

                if(!sbiLevelEdgeFilter.accept(downIterator)){
                    continue;
                }

                int adjNode = downIterator.getAdjNode();
                double weight = calcWeight(downIterator, current, true, true);
                boolean isVirtualAdj = isVirtual(adjNode);
                int levelAdj = getLevel(adjNode);
                int origEdgeId = getOrigEdgeIdV2(downIterator, true);

                if (weight == Double.POSITIVE_INFINITY || isVirtualAdj) {
                    continue;
                }


                double distance = calcDistance(downIterator, current, true);
                long time = calcTime(downIterator, current, true, true);
                sbiAlgo.addInitialOutVertex(adjNode,closestNode,idx,weight,time,distance,false,
                        origEdgeId, downIterator.getOrigEdgeFirst(), downIterator.getOrigEdgeLast());


                heap.add(new RankedNode(adjNode,levelAdj,false));
                traversed.add(adjNode);
                noAccessibleNodes = false;
            }

            heap.add(new RankedNode(closestNode,getLevel(closestNode),noAccessibleNodes));
            this.traversed.add(closestNode);
        }else{
            MatrixEntry initial = new MatrixEntry(closestNode, 0, 0, 0, getLevel(closestNode));

            Deque<MatrixEntry> queue = new ArrayDeque<>();
            queue.add(initial);

            while (!queue.isEmpty()) {

                MatrixEntry current = queue.poll();
                int baseNode = current.adjNode;

                RoutingCHEdgeIterator outIterator = outEdgeExplorer.setBaseNode(baseNode);
                while (outIterator.next()) {

                    int adjNode = outIterator.getAdjNode();
                    int origEdgeId = getOrigEdgeIdV2(outIterator, false);

                    boolean isVirtualAdj = isVirtual(adjNode);

                    double weight = calcWeight(outIterator, current, false, true);

                    if (Double.isInfinite(weight)) {
                        continue;
                    }

                    double distance = calcDistance(outIterator, current, true);
                    long time = calcTime(outIterator, current, false, true);

                    MatrixEntry entry = new MatrixEntry(outIterator.getEdge(), origEdgeId, adjNode, current.adjNode, getLevel(adjNode), weight, time, distance);

                    if (isVirtualAdj && adjNode != closestNode && !processed.contains(adjNode)) {
                        processed.add(adjNode);
                        queue.add(entry);
                    } else if (adjNode != closestNode && !processed.contains(adjNode)) {
                        heap.add(new RankedNode(adjNode,getLevel(adjNode),false));
                        this.traversed.add(adjNode);
                        sbiAlgo.addInitialOutVertex(adjNode,closestNode,idx,weight,time,distance,false,
                                outIterator.getOrigEdge(), outIterator.getOrigEdgeFirst(), outIterator.getOrigEdgeLast());
                    }

                }
            }

        }

    }

    private void backward(){

        while (!heap.isEmpty()) {
            RankedNode current = heap.poll();
            int currentNode = current.node;
            boolean isTerminal = current.noAccessibleNodes;

            NodeOutVertices currentNodeOuts =
                    obtainOutVerticesForCurrentNode(currentNode,inEdgeExplorerNoVirtual,true,isTerminal);
            simultaneousBucketInitialization(currentNode,true,currentNodeOuts);

        }
    }

    private void forward() {

        while (!heap.isEmpty()) {
            RankedNode current = heap.poll();
            int currentNode = current.node;
            boolean isTerminal = current.noAccessibleNodes;

            NodeOutVertices currentNodeOuts2 =
                    obtainOutVerticesForCurrentNode(currentNode,outEdgeExplorerNoVirtual,false,isTerminal);
            simultaneousBucketInitialization(currentNode,false,currentNodeOuts2);

        }
    }

    private void simultaneousBucketInitialization(int currentNode,boolean reverse,NodeOutVertices outs) {

        sbiAlgo.addOuts(currentNode,outs,reverse);
        if(!reverse){
            sbiAlgo.findRoutes(currentNode,outs);
        }

    }

    private NodeOutVertices obtainOutVerticesForCurrentNode(int currentNode,
                                                            RoutingCHEdgeExplorer explorer,
                                                            boolean reverse, boolean noAccessibleNodes){
        NodeOutVertices outs = new NodeOutVertices();

        RoutingCHEdgeIterator iter = explorer.setBaseNode(currentNode);
        while (iter.next()) {

            if (!this.sbiLevelEdgeFilter.accept(iter) && !noAccessibleNodes) {
                continue;
            }

            double weight = iter.getWeight(reverse);

            if (Double.isInfinite(weight)) {
                continue;
            }

            int origEdgeId = getOrigEdgeId(iter,reverse);

            int adjNode = iter.getAdjNode();
            long time = iter.getTime(reverse);
            double distance = iter.getDistance();
            Vertex v = new Vertex(currentNode,iter.getEdge(),adjNode,origEdgeId,weight,time,distance, iter.getOrigEdgeFirst(),iter.getOrigEdgeLast());
            outs.add(v);

            //Add to the heap
            if(!traversed.contains(adjNode)){
                traversed.add(adjNode);
                RankedNode rankedNode = new RankedNode(adjNode,getLevel(adjNode),false);
                heap.add(rankedNode);
            }
        }

        return outs;

    }
}