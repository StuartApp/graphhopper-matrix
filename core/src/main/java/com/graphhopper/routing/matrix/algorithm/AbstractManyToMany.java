package com.graphhopper.routing.matrix.algorithm;

import com.carrotsearch.hppc.*;
import com.carrotsearch.hppc.procedures.IntObjectProcedure;
import com.carrotsearch.hppc.procedures.IntProcedure;
import com.graphhopper.coll.GHIntObjectHashMap;
import com.graphhopper.routing.matrix.*;
import com.graphhopper.routing.querygraph.QueryRoutingCHGraph;
import com.graphhopper.routing.weighting.Weighting;
import com.graphhopper.storage.*;
import com.graphhopper.storage.index.Snap;
import com.graphhopper.util.PairingUtils;

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

    protected IntObjectMap<IntObjectMap<BucketEntry>> backwardBuckets;
    protected IntObjectMap<IntObjectMap<BucketEntry>> forwardBuckets;

    IntObjectMap<MatrixEntry> bestWeightMap;

    protected boolean alreadyRun = false;

    protected int size;
    protected int maxVisitedNodes = Integer.MAX_VALUE;
    protected int visitedNodes = 0;

    protected int maxNodes;

    protected PriorityQueue<RankedNode> heap;

    protected LongDoubleMap tentativeWeights;

    IntObjectMap<ObjectArrayList<Vertex>> verticesNode;

    Vertex minVertex;


    IntSet traversed;

    boolean DEBUG = true;

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

        this.bestWeightMap = new IntObjectHashMap<>();

        this.sbiLevelEdgeFilter = new CHEdgeFilter() {

            @Override
            public boolean accept(RoutingCHEdgeIteratorState edgeState) {

                int base = edgeState.getBaseNode();
                int adj = edgeState.getAdjNode();

                // always accept virtual edges, see #288
                if (base >= maxNodes || adj >= maxNodes) return true;

                // minor performance improvement: shortcuts in wrong direction are disconnected, so no need to exclude them
                if (edgeState.isShortcut()) return true;

                //if (base == adj) return false;

                return graph.getLevel(base) <= graph.getLevel(adj);
            }
        };

        this.size = Math.min(Math.max(200, graph.getNodes() / 10), 150_000);
        this.backwardBuckets = new GHIntObjectHashMap<>(size);
        this.forwardBuckets = new GHIntObjectHashMap<>(size);

        this.heap = new PriorityQueue<>(size);
        this.tentativeWeights = new LongDoubleHashMap(100);
        this.verticesNode = new IntObjectHashMap<>();

    }

    @Override
    public DistanceMatrix calcMatrix(List<Snap> sources, List<Snap> targets) {

        checkAlreadyRun();

        DistanceMatrix matrix = new DistanceMatrix(sources.size(), targets.size());
        IntObjectMap<IntArrayList> targetIdxsNodes = new GHIntObjectHashMap<>(targets.size());

        //Backward
        int idxTarget = 0;
        while (idxTarget < targets.size()) {

            int targetClosestNode = targets.get(idxTarget).getClosestNode();

            //Find Buckets
            findInitialNodesBackward(targets.get(idxTarget), idxTarget);

            //Avoid iterate over the same node two times
            if (!targetIdxsNodes.containsKey(targetClosestNode)) {
                IntArrayList a = new IntArrayList();
                a.add(idxTarget);
                targetIdxsNodes.put(targetClosestNode, a);
            } else {
                targetIdxsNodes.get(targetClosestNode).add(idxTarget);
            }

            idxTarget++;
        }

        if(DEBUG){
            System.out.println("############################## BACKWARD");
        }

        backward(targets);

        //Reset collections for forward
        this.heap.clear();
        this.verticesNode.clear();
        this.traversed.clear();

        //Forward
        int idxSource = 0;
        IntIntMap processedSources = new IntIntHashMap();
        List<IntIntPair> cloneResults = new ArrayList<>();
        while (idxSource < sources.size()) {
            //If we have n sources pointing to the same closestNode, we don't want to calculate both, we will copy the result of
            //one to the other
            int closestNode = sources.get(idxSource).getClosestNode();
            if (processedSources.containsKey(closestNode)) {
                cloneResults.add(new IntIntPair(processedSources.get(closestNode), idxSource));
            } else {
                findInitialNodesForward(sources.get(idxSource), idxSource, targetIdxsNodes, matrix);
                processedSources.put(closestNode, idxSource);
            }

            idxSource++;
        }

        if(DEBUG){
            System.out.println("############################## FORWARD");
        }

        forward(processedSources, matrix, targetIdxsNodes);

        cloneResults.stream().forEach(pair -> matrix.copyResult(pair.source, pair.target));

        return matrix;
    }

    protected void checkAlreadyRun() {
        if (alreadyRun) throw new IllegalStateException("Create a new instance per call");
        alreadyRun = true;
    }

    protected int getOtherNode(int edge, int node) {
        return graph.getBaseGraph().getOtherNode(edge, node);
    }

    protected abstract int getTraversalId(RoutingCHEdgeIteratorState edge, Boolean reverse);

    protected abstract int getOrigEdgeId(RoutingCHEdgeIteratorState edge, boolean reverse);

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

        System.out.println("Backward Init ********************** : " + closestNode + " virtual:" + closestIsVirtual);

        //During the initialization, we add pair (t, 0) to B(t) for each t ∈ T. By adding
        //such a pair, we indicate that t can be reached from t, with a shortest path of length zero.
        IntObjectMap<BucketEntry> bucketTargets = backwardBuckets.get(closestNode);
        if(bucketTargets == null){
            bucketTargets = new IntObjectHashMap<>();
            backwardBuckets.put(closestNode,bucketTargets);
        }

        bucketTargets.put(closestNode, new BucketEntry(0, 0, 0, idx));

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

                if (weight == Double.POSITIVE_INFINITY || isVirtualAdj) {
                    continue;
                }
                System.out.println("Node: " + adjNode + " w:" + weight + " virtual:" + isVirtualAdj + " level:" + levelAdj);


                double distance = calcDistance(downIterator, current, true);
                long time = calcTime(downIterator, current, true, true);
                addVertices(downIterator.getEdge(),closestNode, adjNode, weight, time,distance, -1, -1 ,-1);

                IntObjectMap<BucketEntry> buckets = backwardBuckets.get(adjNode);
                if(buckets == null){
                    buckets = new IntObjectHashMap<>();
                    backwardBuckets.put(adjNode,buckets);
                }
                buckets.put(closestNode, new BucketEntry(weight, time, distance, idx));
                System.out.println("Saving Bucket:" + adjNode + " - " + closestNode);

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

                    int origEdgeId = getOrigEdgeId(downIterator, true);
                    double weight = calcWeight(downIterator, current, true, true);
                    boolean isVirtualAdj = isVirtual(adjNode);

                    System.out.println("Node: " + adjNode + " w:" + weight + " virtual:" + isVirtualAdj);

                    if (weight == Double.POSITIVE_INFINITY) {
                        continue;
                    }


                    double distance = calcDistance(downIterator, current, true);
                    long time = calcTime(downIterator, current, true, true);

                    MatrixEntry entry = new MatrixEntry(downIterator.getEdge(), origEdgeId, adjNode, current.adjNode, getLevel(adjNode), weight, time, distance);

                    if (isVirtualAdj && adjNode != closestNode && !processed.contains(adjNode)) {
                        processed.add(adjNode);
                        queue.add(entry);
                        addVertices(downIterator.getEdge(),closestNode, adjNode, weight, time,distance, -1, -1 ,-1);
                    } else if (adjNode != closestNode && !processed.contains(adjNode)) {
                        if(!this.traversed.contains(adjNode)){
                            heap.add(new RankedNode(adjNode,getLevel(adjNode),false));
                            this.traversed.add(adjNode);
                        }

                        addVertices(downIterator.getEdge(),closestNode, adjNode, weight, time,distance, -1, -1, -1);
                    }

                    IntObjectMap<BucketEntry> buckets = backwardBuckets.get(adjNode);
                    if(buckets == null){
                        buckets = new IntObjectHashMap<>();
                        backwardBuckets.put(adjNode,buckets);
                    }
                    buckets.put(closestNode, new BucketEntry(weight, time, distance, idx));
                    System.out.println("Saving Bucket:" + adjNode + " - " + closestNode);
                }
            }
        }


    }

    private void findInitialNodesForward(Snap snap, int idx,
                                         IntObjectMap<IntArrayList> targets, DistanceMatrix dm) {

        int closestNode = snap.getClosestNode();
        IntSet processed = new IntHashSet();
        boolean closestIsVirtual = isVirtual(closestNode);
        int closestLevel = getLevel(closestNode);

        System.out.println("Forward Init ********************** : " + closestNode + " virtual:" + closestIsVirtual);

        //During the initialization, we add pair (t, 0) to B(t) for each t ∈ T. By adding
        //such a pair, we indicate that t can be reached from t, with a shortest path of length zero.

        IntObjectMap<BucketEntry> bucketTargets = forwardBuckets.get(closestNode);
        if(bucketTargets == null){
            bucketTargets = new IntObjectHashMap<>();
            forwardBuckets.put(closestNode,bucketTargets);
        }
        bucketTargets.put(closestNode, new BucketEntry(0, 0, 0, idx));

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

                if (weight == Double.POSITIVE_INFINITY || isVirtualAdj) {
                    continue;
                }
                System.out.println("Node: " + adjNode + " w:" + weight + " virtual:" + isVirtualAdj + " level:" + levelAdj);


                double distance = calcDistance(downIterator, current, true);
                long time = calcTime(downIterator, current, true, true);
                addVertices(downIterator.getEdge(),closestNode, adjNode, weight, time,distance, -1, -1 ,-1);

                IntObjectMap<BucketEntry> buckets = forwardBuckets.get(adjNode);
                if(buckets == null){
                    buckets = new IntObjectHashMap<>();
                    forwardBuckets.put(adjNode,buckets);
                }
                buckets.put(closestNode, new BucketEntry(weight, time, distance, idx));
                System.out.println("Saving Bucket:" + adjNode + " - " + closestNode);

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
                    int origEdgeId = getOrigEdgeId(outIterator, false);

                    boolean isVirtualAdj = isVirtual(adjNode);

                    double weight = calcWeight(outIterator, current, false, true);

                    System.out.println(" Iter:" + adjNode + " w:" + weight);

                    if (Double.isInfinite(weight)) {
                        continue;
                    }

                    if (targets.containsKey(adjNode)) {
                        long uniqueId = PairingUtils.pair(closestNode, adjNode);

                        final double savedWeight = tentativeWeights.get(uniqueId);

                        if ((savedWeight == 0.0 || (weight < savedWeight)) && closestNode != adjNode) {

                            double distance = calcDistance(outIterator, current, true);
                            long time = calcTime(outIterator, current, false, true);

                            tentativeWeights.put(uniqueId, weight);

                            targets.get(adjNode).forEach(new IntProcedure() {
                                @Override
                                public void apply(int i) {
                                    dm.setCell(idx, i, distance, time);
                                }
                            });
                        }
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
                        addVertices(outIterator.getEdge(),closestNode, adjNode, weight, time,distance, origEdgeId, outIterator.getOrigEdgeFirst(), outIterator.getOrigEdgeLast());
                    }

                    IntObjectMap<BucketEntry> buckets = forwardBuckets.get(adjNode);
                    if(buckets == null){
                        buckets = new IntObjectHashMap<>();
                        forwardBuckets.put(adjNode,buckets);
                    }
                    buckets.put(closestNode, new BucketEntry(weight, time, distance, idx));
                }
            }

        }

    }


    private NodeInVertices obtainInsVerticesForCurrentNode(int currentNode, IntObjectMap<IntObjectMap<BucketEntry>> buckets){

        NodeInVertices nodeInVertices = new NodeInVertices();
        ObjectArrayList<Vertex> inVertices = verticesNode.get(currentNode);

        if(inVertices != null){
            Object[] objectBuffer = inVertices.buffer;
            int inVerticesSize = inVertices.size();

            for (int i = 0; i < inVerticesSize; i++) {

                Vertex in = (Vertex) objectBuffer[i];
                //System.out.println("IN: " + in);
                nodeInVertices.addOrigEdgeIds(in.origEdgeId);

                IntObjectMap<BucketEntry> b = buckets.get(in.baseNode);
                if(b != null){
                    b.forEach(new IntObjectProcedure<BucketEntry>() {
                        @Override
                        public void apply(int i, BucketEntry bucketEntry) {
                            //System.out.println(in.baseNode + " <-> " + bucketEntry);
                            nodeInVertices.add(in,bucketEntry,i);
                        }
                    });
                }
            }
        }

        return nodeInVertices;

    }

    private NodeOutVertices obtainOutVerticesForCurrentNode(int currentNode,
                                                            NodeInVertices ins,
                                                            RoutingCHEdgeExplorer explorer,
                                                            boolean reverse,boolean noAccessibleNodes){
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

            //If an in vertex origEdgeId is equal of this vertex origEdgeId we discard it
            if(ins.containsEdge(origEdgeId)){
                continue;
            }

            int adjNode = iter.getAdjNode();
            long time = iter.getTime(reverse);
            double distance = iter.getDistance();
            Vertex v = new Vertex(currentNode,iter.getEdge(),adjNode,origEdgeId,weight,time,distance, iter.getOrigEdgeFirst(),iter.getOrigEdgeLast());
            outs.add(v);
            //System.out.println("Out: " + v);
        }

        return outs;

    }



    private void copyDistancesToTerminalsToCurrentNodeBuckets(int currentNode, NodeInVertices in,
                                                              boolean reverse,
                                                              IntObjectMap<IntArrayList> targets,
                                                              DistanceMatrix dm,
                                                              IntObjectMap<IntObjectMap<BucketEntry>> buckets){

        IntObjectMap<BucketEntry> currentNodeBuckets = buckets.get(currentNode);
        if(currentNodeBuckets == null){
            currentNodeBuckets = new IntObjectHashMap();
            buckets.put(currentNode, currentNodeBuckets);
        }else{
            System.out.println("Current: " + currentNodeBuckets);
        }

        for (int i = 0; i < in.minSize(); i++) {

            VertexTerminal vt = in.minValues()[i];

            double w = vt.weight;
            double d = vt.distance;
            long  t = vt.time;

            BucketEntry entry = new BucketEntry(w, t, d, vt.terminalIdx);
            currentNodeBuckets.put(vt.terminal,entry);

            if(DEBUG){
                System.out.println("Saving Bucket:" + currentNode + " - " +  vt.terminal +
                        " ->  w:" + entry.weight + " t: " + entry.time + " d: " + entry.distance);
            }

            if (!reverse) {
                saveBestPath(vt.terminal, vt.terminalIdx, currentNode, entry, targets, dm);
            }
        }

    }

    private void saveAccessibleOutVertices(int currentNode,NodeInVertices ins,NodeOutVertices outs, boolean reverse){

        IntObjectMap<OutVertexCost> mins = new IntObjectHashMap();

        if(outs.hasOutValues()){

            if(ins.isEmpty()){
                //Save all the outputs
                TurnCost tc = new TurnCost(0,0);
                for (int i = 0; i < outs.outsSize(); i++) {
                    Vertex out = outs.outValues()[i];
                    OutVertexCost min = mins.get(out.origEdgeId);
                    if(min == null || min.cost() > out.weight){
                        OutVertexCost outCost = new OutVertexCost(out,tc);
                        mins.put(out.origEdgeId,outCost);
                    }
                }
            }else{

                System.out.println("Ins size: " + currentNode + " " + ins.size());
                for (int i = 0; i < ins.size(); i++) {

                    VertexTerminal in = ins.values()[i];
                    int inOrigEdgeId = reverse ? in.origEdgeFirst : in.origEdgeLast;

                    System.out.println("In " + in);

                    for (int ii = 0; ii < outs.outsSize(); ii++) {

                        Vertex out = outs.outValues()[ii];
                        int outOrigEdgeId = reverse ? out.origEdgeLast : out.origEdgeFirst;

                        double cost = reverse
                                ? graph.getTurnWeight(outOrigEdgeId, currentNode, inOrigEdgeId)
                                : graph.getTurnWeight(inOrigEdgeId, currentNode, outOrigEdgeId);

                        if(isVertexAccessible(cost)){
                            TurnCost tc = new TurnCost(in.weight,cost);
                            OutVertexCost outCost = new OutVertexCost(out,tc);
                            OutVertexCost min = mins.get(out.origEdgeId);
                            if(min == null || min.cost() > outCost.cost()){
                                mins.put(out.origEdgeId,outCost);
                            }
                        }

                            //Check Selfs
                            if(outs.hasSelfLoops()){
                                for (int iii = 0; iii < outs.selfSize(); iii++) {
                                    Vertex self = outs.selfValues()[iii];
                                    int selfOrigEdgeId = reverse ? self.origEdgeLast : self.origEdgeFirst;

                                    double costInToSelf = reverse
                                            ? graph.getTurnWeight(selfOrigEdgeId, currentNode, inOrigEdgeId)
                                            : graph.getTurnWeight(inOrigEdgeId, currentNode, selfOrigEdgeId);

                                    if(isVertexAccessible(costInToSelf)){
                                        double costSelfToOut = reverse
                                                ? graph.getTurnWeight(outOrigEdgeId, currentNode,selfOrigEdgeId )
                                                : graph.getTurnWeight(selfOrigEdgeId, currentNode, outOrigEdgeId);

                                        if(isVertexAccessible((costSelfToOut))){
                                            Vertex accessibleOut = out.withSelf(self);
                                            TurnCost tc = new TurnCost(accessibleOut.weight,costInToSelf + costSelfToOut);
                                            OutVertexCost outCost = new OutVertexCost(out,tc);

                                            OutVertexCost min = mins.get(out.origEdgeId);
                                            if(min == null || min.cost() > outCost.cost()){
                                                mins.put(out.origEdgeId,outCost);
                                            }
                                        }
                                    }
                                }
                            }
                    }
                }
            }

            //Save mins outs
            mins.forEach(new IntObjectProcedure<OutVertexCost>() {
                @Override
                public void apply(int i, OutVertexCost outCost) {
                    Vertex out = outCost.getOut();
                    System.out.println("Out: " + out);

                    ObjectArrayList<Vertex> adjNodeIns = verticesNode.get(out.adjNode);
                    if(adjNodeIns == null){
                        adjNodeIns = new ObjectArrayList<>();
                        verticesNode.put(out.adjNode,adjNodeIns);
                    }
                    adjNodeIns.add(out);

                    //Add to the heap
                    if(!traversed.contains(out.adjNode)){
                        traversed.add(out.adjNode);
                        RankedNode rankedNode = new RankedNode(out.adjNode,getLevel(out.adjNode),false);
                        heap.add(rankedNode);
                    }
                }
            });
        }
    }

    private void simultaneousBucketInitialization(int currentNode,
                                                  NodeInVertices ins, NodeOutVertices outs,
                                                  boolean reverse,
                                                  IntObjectMap<IntObjectMap<BucketEntry>> buckets,
                                                  IntObjectMap<IntArrayList> targets,
                                                  DistanceMatrix dm) {

        if(ins.isNotEmpty()) {
            //Copy distances to buckets in the current node
            copyDistancesToTerminalsToCurrentNodeBuckets(currentNode, ins, reverse, targets, dm, buckets);
        }

        //Save Accessible Out Vertices
        saveAccessibleOutVertices(currentNode, ins, outs, reverse);
    }

    private boolean isVertexAccessible(double cost) {
        return Double.isFinite(cost);
    }

    private long turnCostInMillis(double turnCost){
        return (long) turnCost * 1000;
    }


    private void addVertices(int edge,int baseNode, int adjNode, double weight, long time, double distance, int origEdgeId, int origEdgeFirst, int origEdgeLast) {

        ObjectArrayList<Vertex> vertexs = this.verticesNode.get(adjNode);

        if(vertexs == null){
            vertexs = new ObjectArrayList<>();
            this.verticesNode.put(adjNode,vertexs);
        }

        Vertex v = new Vertex(baseNode,edge,adjNode,origEdgeId,weight,time,distance,origEdgeFirst,origEdgeLast);

        if(!v.isSelfLoop()){
            System.out.println(adjNode + " Add Vertex:" + v);
            vertexs.add(v);
        }
    }

    private void backward(List<Snap> targets){

        if(DEBUG){
            System.out.println("######### BACKWARD ################");
        }

        while (!heap.isEmpty()) {
            RankedNode current = heap.poll();
            int currentNode = current.node;
            boolean isTerminal = current.noAccessibleNodes;

            if(DEBUG){
                System.out.println("## " + currentNode + " order:" + current.level);
            }

            NodeInVertices currentNodeIns = obtainInsVerticesForCurrentNode(currentNode,backwardBuckets);
            NodeOutVertices currentNodeOuts =
                    obtainOutVerticesForCurrentNode(currentNode,currentNodeIns,inEdgeExplorerNoVirtual,true,isTerminal);

           simultaneousBucketInitialization(currentNode,currentNodeIns,currentNodeOuts,true,backwardBuckets,null,null);
        }
    }

    private void forward(IntIntMap processedSources, DistanceMatrix dm, IntObjectMap<IntArrayList> targets) {

        if(DEBUG){
            System.out.println("######### FORWARD ################");
        }

        while (!heap.isEmpty()) {
            RankedNode current = heap.poll();
            int currentNode = current.node;
            boolean isTerminal = current.noAccessibleNodes;

            if(DEBUG){
                System.out.println("## " + currentNode);
            }



            NodeInVertices currentNodeIns = obtainInsVerticesForCurrentNode(currentNode,forwardBuckets);
            NodeOutVertices currentNodeOuts =
                    obtainOutVerticesForCurrentNode(currentNode,currentNodeIns,outEdgeExplorerNoVirtual,false, isTerminal);

            if(currentNodeIns.isEmpty()){
                if(processedSources.containsKey(currentNode)){
                    saveBestPath(currentNode,processedSources.get(currentNode),currentNode,new BucketEntry(0,0,0,0),targets,dm);
                }
            }

            simultaneousBucketInitialization(currentNode,currentNodeIns,currentNodeOuts,false,forwardBuckets,targets,dm);
        }
    }

    private void saveBestPath(int sourceNode, int idxSource, int currNode, BucketEntry forwardEntry,
                              IntObjectMap<IntArrayList> targets, DistanceMatrix dm) {

        final IntObjectMap<BucketEntry> backwardEntries = backwardBuckets.get(currNode);


        if (backwardEntries != null) {

            backwardEntries.forEach(new IntObjectProcedure<BucketEntry>() {
                @Override
                public void apply(int target, BucketEntry entry) {

                    if (sourceNode != target) {

                        long uniqueId = PairingUtils.pair(sourceNode, target);

                        final double savedWeight = tentativeWeights.get(uniqueId);
                        final double currentWeight = forwardEntry.weight + entry.weight;

                        if (savedWeight == 0.0 || (currentWeight < savedWeight)) {

                            if(DEBUG){
                                System.out.println(idxSource + "- " + entry.idx + " Forward: " + forwardEntry.weight + " Backward: " + entry.weight);
                                System.out.println(currNode + " $$$$$ Current:" + currentWeight + " Saved: " + savedWeight);
                            }

                            final long time = forwardEntry.time + entry.time;
                            final double distance = forwardEntry.distance + entry.distance;
                            tentativeWeights.put(uniqueId, currentWeight);

                            if(DEBUG){
                                System.out.println(" $$$$$ Saved Distance:" + distance + " Time: " + time);
                            }

                            IntArrayList targetsIdxs = targets.get(target);
                            int[] buffer = targetsIdxs.buffer;
                            int size = targetsIdxs.size();

                            for (int i = 0; i < size; i++) {
                                dm.setCell(idxSource, buffer[i], distance, time);
                            }
                        }
                    }
                }
            });
        }
    }
}