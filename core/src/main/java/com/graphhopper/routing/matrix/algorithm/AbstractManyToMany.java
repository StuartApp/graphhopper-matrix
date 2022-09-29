package com.graphhopper.routing.matrix.algorithm;

import com.carrotsearch.hppc.*;
import com.carrotsearch.hppc.procedures.IntIntProcedure;
import com.carrotsearch.hppc.procedures.IntObjectProcedure;
import com.carrotsearch.hppc.procedures.IntProcedure;
import com.carrotsearch.hppc.procedures.ObjectProcedure;
import com.graphhopper.coll.GHIntObjectHashMap;
import com.graphhopper.routing.matrix.*;
import com.graphhopper.routing.querygraph.QueryRoutingCHGraph;
import com.graphhopper.routing.weighting.Weighting;
import com.graphhopper.storage.*;
import com.graphhopper.storage.index.Snap;
import com.graphhopper.util.EdgeIteratorState;
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
    protected Weighting weightingNoVirtualNodes;

    protected RoutingCHEdgeExplorer inEdgeExplorer;
    protected RoutingCHEdgeExplorer inEdgeExplorerNoVirtual;
    protected RoutingCHEdgeExplorer outEdgeExplorer;
    protected RoutingCHEdgeExplorer outEdgeExplorerNoVirtual;

    protected CHEdgeFilter sbiLevelEdgeFilter;

    protected IntObjectMap<IntObjectMap<BucketEntry>> backwardBuckets;
    protected IntObjectMap<IntObjectMap<BucketEntry>> forwardBuckets;

    protected boolean alreadyRun = false;

    protected int size;
    protected int maxVisitedNodes = Integer.MAX_VALUE;
    protected int visitedNodes = 0;

    protected int maxNodes;

    protected PriorityQueue<RankedNode> heap;

    protected LongDoubleMap tentativeWeights;

    protected IntSet visited = new IntHashSet();

    IntObjectMap<ObjectArrayList<SBIEntry>> inVertices;
    IntObjectMap<ObjectArrayList<PruningVertex>> prunningVertices;
    IntObjectMap<ObjectArrayList<SBIEntry>> outVertices;

    IntDoubleMap bestTraversalWeight;

    IntSet nodesAdded;
    IntIntHashMap terminals;

    public AbstractManyToMany(QueryRoutingCHGraph graph, RoutingCHGraph graphNoVirtualNodes) {

        this.graph = graph;
        this.weighting = graph.getWrappedWeighting();
        this.weightingNoVirtualNodes = graphNoVirtualNodes.getWeighting();
        this.inEdgeExplorer = graph.createInEdgeExplorer();
        this.outEdgeExplorer = graph.createOutEdgeExplorer();
        this.maxNodes = graph.getBaseGraph().getBaseGraph().getNodes();

        this.graphNoVirtualNodes = graphNoVirtualNodes;
        this.inEdgeExplorerNoVirtual = graphNoVirtualNodes.createInEdgeExplorer();
        this.outEdgeExplorerNoVirtual = graphNoVirtualNodes.createOutEdgeExplorer();

        this.nodesAdded = new IntHashSet();
        this.terminals = new IntIntHashMap();
        this.bestTraversalWeight = new IntDoubleHashMap();

        this.sbiLevelEdgeFilter = new CHEdgeFilter() {

            @Override
            public boolean accept(RoutingCHEdgeIteratorState edgeState) {

                int base = edgeState.getBaseNode();
                int adj = edgeState.getAdjNode();

                if (base == adj) return false;

                // always accept virtual edges, see #288
                if (base >= maxNodes || adj >= maxNodes) return true;

                // minor performance improvement: shortcuts in wrong direction are disconnected, so no need to exclude them
                if (edgeState.isShortcut()) return true;

                return graph.getLevel(base) <= graph.getLevel(adj);
            }
        };

        this.size = Math.min(Math.max(200, graph.getNodes() / 10), 150_000);
        this.backwardBuckets = new GHIntObjectHashMap<>(size);
        this.forwardBuckets = new GHIntObjectHashMap<>(size);

        this.heap = new PriorityQueue<>(size);
        this.tentativeWeights = new LongDoubleHashMap(100);

        this.inVertices = new IntObjectHashMap<>(this.graph.getNodes());
        this.prunningVertices = new IntObjectHashMap<>(this.graph.getNodes());
        this.outVertices = new IntObjectHashMap<>(this.graph.getNodes());

    }

    @Override
    public DistanceMatrix calcMatrix(List<Snap> sources, List<Snap> targets) {

        //TODO - To debug purposes
        /*
        #######################################
Forward
748099 -> 251.02300000000002 : 251.02300000000002
2096490 -> 415.50300000000004 : 164.48000000000002
1730181 -> 555.494 : 139.99099999999999
2094484 -> 561.2948955243202 : 5.80089552432014
Backward
748099 -> 1292.48489552432 : 731.1899999999999
2904594 -> 6495.7658955243205 : 5203.281000000001
537 -> 8233.81689552432 : 1738.0509999999995
2966295 -> 9245.26689552432 : 1011.4500000000007
139912 -> 9642.285895524321 : 397.01900000000023
 -> 9654.79689552432 : 12.510999999998603

         */
/*
        int base = 2966295;
        RoutingCHEdgeIterator iter = inEdgeExplorerNoVirtual.setBaseNode(base);
        while (iter.next()) {

            if(sbiLevelEdgeFilter.accept(iter)){
                int adjNode = iter.getAdjNode();
                double weight = calcWeight(iter, true, base);
                double distance = calcDistance(iter);
                long time = calcTime(iter, true, base);
                int originId = getOrigEdgeId(iter,true);
                int traversalId = getTraversalId(iter, originId,true);

                if(adjNode == 537){
                    System.out.println("#### "+ base + " -> " + adjNode + " distance: " + distance + "(" + weight + ")" + " --- " + iter.isShortcut() + " -->" + traversalId + "(" + iter.getEdge() + " - " + iter.getOrigEdgeFirst() + "/" + iter.getOrigEdgeLast());
                }


            }

        }
        //System.exit(1);


 */


                //

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

        //System.out.println("############################## BACKWARD");

        simultaneousBucketInitializationBackward(targets);

        //Reset collections for forward
        this.heap.clear();
        this.nodesAdded.clear();
        this.terminals.clear();
        this.prunningVertices.clear();

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

        //System.out.println("############################## FORWARD");

        simultaneousBucketInitializationForward(sources, matrix, targetIdxsNodes);

        cloneResults.stream().forEach(pair -> matrix.copyResult(pair.source, pair.target));

        return matrix;
    }

    protected void checkAlreadyRun() {
        if (alreadyRun) throw new IllegalStateException("Create a new instance per call");
        alreadyRun = true;
    }

    protected abstract int getTraversalId(RoutingCHEdgeIteratorState edge, int origEdgeId, Boolean reverse);

    protected abstract int getTraversalId(EdgeIteratorState state,Boolean reverse);

    protected int getOrigEdgeId(RoutingCHEdgeIteratorState edge, boolean reverse) {
        return reverse ? edge.getOrigEdgeFirst() : edge.getOrigEdgeLast();
    }

    protected int getOrigEdgeId(RoutingCHEdgeIterator edge, boolean reverse) {
        return reverse ? edge.getOrigEdgeFirst() : edge.getOrigEdgeLast();
    }

    protected abstract double calcWeight(RoutingCHEdgeIteratorState iter, int incomingEdge, boolean reverse);
    protected abstract double calcWeightNoVirtual(RoutingCHEdgeIteratorState iter, int incomingEdge, boolean reverse);


    protected abstract long calcTime(RoutingCHEdgeIteratorState iter,int incomingEdge, boolean reverse);
    protected abstract long calcTimeNoVirtual(RoutingCHEdgeIteratorState iter,int incomingEdge, boolean reverse);

    protected int getIncomingEdge(SBIEntry entry) {
        return entry.incEdge;
    }

    protected abstract double calcDistance(RoutingCHEdgeIteratorState iter);

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

    private void findInitialNodesBackward(Snap snap, int idx) {

        System.out.println("Backward Init **********************");

        int closestNode = snap.getClosestNode();
        IntSet processed = new IntHashSet();
        int snapTraversalId = getTraversalId(snap.getClosestEdge(), true);

        System.out.println("Closest Node Backward:" + closestNode);

        //During the initialization, we add pair (t, 0) to B(t) for each t ∈ T. By adding
        //such a pair, we indicate that t can be reached from t, with a shortest path of length zero.
        IntObjectMap<BucketEntry> bucketTargets = new GHIntObjectHashMap<>();
        bucketTargets.put(closestNode, new BucketEntry(0, 0, 0, idx));
        backwardBuckets.put(closestNode, bucketTargets);

        //Closest Node is Virtual
        if (isVirtual(closestNode)) {

            Deque<VirtualNodeEntry> queue = new ArrayDeque<>();
            queue.add(new VirtualNodeEntry(snapTraversalId,closestNode,-1, 0, 0, 0));

            while (!queue.isEmpty()) {

                VirtualNodeEntry current = queue.poll();
                int baseNode = current.node;

                RoutingCHEdgeIterator downIterator = inEdgeExplorer.setBaseNode(baseNode);
                while (downIterator.next()) {

                    int adjNode = downIterator.getAdjNode();
                    int originId = getOrigEdgeId(downIterator,true);
                    int traversalId = getTraversalId(downIterator, originId,true);

                    System.out.println( baseNode + " <--" + "(" + downIterator.getEdge() + ")-- " + adjNode + " Traversal: " + traversalId);

                    boolean isVirtualAdj = isVirtual(adjNode);

                    if (isVirtualAdj && adjNode != closestNode && !processed.contains(adjNode)) {
                        processed.add(adjNode);
                        double weight = calcWeight(downIterator, current.edge, true) + current.weight;
                        double distance = calcDistance(downIterator) + current.distance;
                        long time = calcTime(downIterator, current.edge,true) + current.time;
                        queue.add(new VirtualNodeEntry(traversalId,adjNode, downIterator.getEdge(), weight, time, distance));
                    } else if (adjNode != closestNode && !processed.contains(adjNode)) {
                        heap.add(new RankedNode(traversalId,downIterator.getEdge(),adjNode, graph.getLevel(adjNode),true));
                        nodesAdded.add(traversalId);

                        //Down Vertices
                        ObjectArrayList<SBIEntry> dVertices = inVertices.get(adjNode);
                        if (dVertices == null) {
                            dVertices = new ObjectArrayList<>();
                            inVertices.put(adjNode, dVertices);
                        }

                        double weight = calcWeight(downIterator, current.edge, true) + current.weight;

                        if (weight < Double.POSITIVE_INFINITY) {
                            double distance = calcDistance(downIterator) + current.distance;
                            long time = calcTime(downIterator, current.edge, true) + current.time;
                            SBIEntry entry = new SBIEntry(traversalId,downIterator.getEdge(),closestNode, weight, time, distance);
                            System.out.println(entry);
                            dVertices.add(entry);

                            IntObjectMap<BucketEntry> nodeTargets = backwardBuckets.get(adjNode);
                            if(nodeTargets == null){
                                nodeTargets = new IntObjectHashMap<>();
                                backwardBuckets.put(adjNode,nodeTargets);
                            }
                            BucketEntry b = new BucketEntry(0,0,0,idx);
                            System.out.println(" Saving Bucket : " + adjNode + " ->" + closestNode + " : " + b);
                            nodeTargets.put(closestNode,b);
                        }
                    }
                }
            }

        } else {
            //Closest Node is not virtual

            int inEdges = 0;
            RoutingCHEdgeIterator uniqueIn = null;
            RoutingCHEdgeIterator downIterator = inEdgeExplorer.setBaseNode(closestNode);
            while (downIterator.next()) {

                boolean accept = this.sbiLevelEdgeFilter.accept(downIterator);
                if (accept) {
                    break;
                } else {
                    uniqueIn = downIterator;
                }
            }

            if (inEdges == 0 && uniqueIn != null) {
                int traversalId = getTraversalId(uniqueIn,snap.getClosestEdge().getEdge(),true);

                heap.add(new RankedNode(traversalId,uniqueIn.getEdge(),uniqueIn.getAdjNode(), graph.getLevel(uniqueIn.getAdjNode()),true));
                nodesAdded.add(uniqueIn.getAdjNode());

                double weight = calcWeight(uniqueIn, snap.getClosestEdge().getEdge() ,true);

                //Down Vertices
                ObjectArrayList<SBIEntry> dVertices = inVertices.get(traversalId);
                if (dVertices == null) {
                    dVertices = new ObjectArrayList<>();
                    inVertices.put(traversalId, dVertices);
                }

                if (weight < Double.POSITIVE_INFINITY) {
                    double distance = calcDistance(uniqueIn);
                    long time = calcTime(uniqueIn, getOrigEdgeId(uniqueIn,true), true);
                    dVertices.add(new SBIEntry(traversalId, downIterator.getAdjNode(), uniqueIn.getEdge(), weight, time, distance));
                }
            } else {
                heap.add(new RankedNode(snapTraversalId,snap.getClosestEdge().getEdge(),closestNode,graph.getLevel(closestNode),true));
                nodesAdded.add(closestNode);
            }
        }
        System.out.println("Heap");
        heap.forEach(System.out::println);

    }

    private void findInitialNodesForward(Snap snap, int idx,
                                         IntObjectMap<IntArrayList> targets, DistanceMatrix dm) {

        int closestNode = snap.getClosestNode();
        IntSet processed = new IntHashSet();
        int snapTraversalId = getTraversalId(snap.getClosestEdge(), false);

        //During the initialization, we add pair (t, 0) to B(t) for each t ∈ T. By adding
        //such a pair, we indicate that t can be reached from t, with a shortest path of length zero.
        IntObjectMap<BucketEntry> bucketTargets = new GHIntObjectHashMap<>();
        bucketTargets.put(closestNode, new BucketEntry(0, 0, 0, idx));
        forwardBuckets.put(closestNode, bucketTargets);

        Deque<VirtualNodeEntry> queue = new ArrayDeque<>();
        queue.add(new VirtualNodeEntry(snapTraversalId,closestNode,snap.getClosestEdge().getEdge(), 0, 0, 0));

        //Process Queue

        while (!queue.isEmpty()) {

            VirtualNodeEntry current = queue.poll();
            int baseNode = current.node;

            RoutingCHEdgeIterator outIterator = outEdgeExplorer.setBaseNode(baseNode);
            while (outIterator.next()) {

                int adjNode = outIterator.getAdjNode();
                int edge = outIterator.getEdge();

                boolean isVirtualAdj = isVirtual(adjNode);

                int originId = getOrigEdgeId(outIterator,true);

                double w = outIterator.getWeight(false); //TODO - Refactor this

                double weight = w + current.weight;

                if (weight < Double.POSITIVE_INFINITY) {


                    if (targets.containsKey(adjNode)) {
                        long uniqueId = PairingUtils.pair(closestNode, adjNode);

                        final double savedWeight = tentativeWeights.get(uniqueId);

                        if ((savedWeight == 0.0 || (weight < savedWeight)) && closestNode != adjNode) {

                            double distance = calcDistance(outIterator) + current.distance;
                            long time = calcTime(outIterator, current.edge, false) + current.time;

                            tentativeWeights.put(uniqueId, weight);

                            targets.get(adjNode).forEach(new IntProcedure() {
                                @Override
                                public void apply(int i) {
                                    dm.setCell(idx, i, distance, time);
                                }
                            });
                        }
                    }

                    if (isVirtualAdj && adjNode != closestNode && !processed.contains(adjNode)) {

                        double distance = calcDistance(outIterator) + current.distance;
                        long time = calcTime(outIterator, snap.getClosestEdge().getEdge() ,false) + current.time;
                        processed.add(adjNode);
                        queue.add(new VirtualNodeEntry(snapTraversalId,adjNode,current.edge, weight, time, distance));
                    } else if (adjNode != closestNode && !processed.contains(adjNode)) {

                        int traversalId = getTraversalId(outIterator, originId,false);

                        this.heap.add(new RankedNode(traversalId,edge,adjNode, graph.getLevel(adjNode),true));
                        this.nodesAdded.add(adjNode);

                        //Out Vertices
                        ObjectArrayList<SBIEntry> nodeOutVertices = outVertices.get(adjNode);
                        if (nodeOutVertices == null) {
                            nodeOutVertices = new ObjectArrayList<>();
                            outVertices.put(adjNode, nodeOutVertices);
                        }

                        double distance = calcDistance(outIterator) + current.distance;
                        long t = outIterator.getTime(false); //TODO - Refactor this
                        long time = t + current.time;

                        nodeOutVertices.add(new SBIEntry(traversalId, outIterator.getEdge(),closestNode,weight, time, distance));

                        IntObjectMap<BucketEntry> nodeTargets = forwardBuckets.get(adjNode);
                        if(nodeTargets == null){
                            nodeTargets = new IntObjectHashMap<>();
                            forwardBuckets.put(adjNode,nodeTargets);
                        }
                        BucketEntry b = new BucketEntry(0,0,0,idx);
                        System.out.println(" Saving Bucket : " + adjNode + " ->" + closestNode + " : " + b);
                        nodeTargets.put(closestNode,b);
                    }
                }
            }
        }
    }

    private void initializeVertices(RankedNode baseNode, RoutingCHEdgeExplorer explorer,
                                    IntObjectMap<ObjectArrayList<SBIEntry>> vertices, boolean reverse) {

        RoutingCHEdgeIterator downIterator = explorer.setBaseNode(baseNode.adjNode);

        while (downIterator.next()) {

            int adjNode = downIterator.getAdjNode();
            int edge = downIterator.getEdge();
            int traversalId = getTraversalId(downIterator, baseNode.edge,reverse);


                int adjRank = graph.getLevel(adjNode);
                boolean accept = this.sbiLevelEdgeFilter.accept(downIterator);

                if ((accept && !nodesAdded.contains(traversalId))) {
                    nodesAdded.add(traversalId);
                    heap.add(new RankedNode(traversalId,edge,adjNode, adjRank,false));
                }

                if (accept) {

                    double weight = calcWeightNoVirtual(downIterator, baseNode.edge, reverse);;

                    System.out.println("    Id:" + traversalId + " --" + edge + "--> " + adjNode + " weight: " + weight );

                    if (weight < Double.POSITIVE_INFINITY) {

                        ObjectArrayList<SBIEntry> dVertices = vertices.get(adjNode);
                        if (dVertices == null) {
                            dVertices = new ObjectArrayList<>();
                            vertices.put(adjNode, dVertices);
                        }

                        double storedWeigh = bestTraversalWeight.get(traversalId);
                        if(storedWeigh == 0.0 || storedWeigh > weight){
                            double distance = calcDistance(downIterator);

                            long time;

                            if(baseNode.virtual){
                                time = downIterator.getTime(reverse);
                            }else{
                                time = calcTimeNoVirtual(downIterator, baseNode.edge,reverse);
                            }

                            dVertices.add(new SBIEntry(baseNode.traversalId, edge, baseNode.adjNode, weight, time, distance));
                            bestTraversalWeight.put(traversalId,weight);
                        }
                    }
                }
        }

    }

    private void initializePruningVertices(RankedNode baseNode, RoutingCHEdgeExplorer explorer, boolean reverse) {

        RoutingCHEdgeIterator upIterator = explorer.setBaseNode(baseNode.adjNode);
        while (upIterator.next()) {
            boolean accept = this.sbiLevelEdgeFilter.accept(upIterator);

            if (accept) {

                int originId = getOrigEdgeId(upIterator,reverse);
                int traversalId = getTraversalId(upIterator, originId,reverse);

                ObjectArrayList<PruningVertex> uVertices = prunningVertices.get(traversalId);
                if (uVertices == null) {
                    uVertices = new ObjectArrayList<>();
                    prunningVertices.put(traversalId, uVertices);
                }

                double weight = calcWeight(upIterator,baseNode.edge ,reverse);
                if (weight < Double.POSITIVE_INFINITY) {
                    uVertices.add(new PruningVertex(traversalId, weight));
                }
            }
        }
    }

    private void applyRetrospectivePruningAlgorithm(RankedNode baseNode, double[] weights,
                                                    IntObjectMap<IntObjectMap<BucketEntry>> buckets) {

        ObjectArrayList<PruningVertex> upList = prunningVertices.get(baseNode.edge); //TODO - Review that

        if (upList != null) {
            upList.forEach(new ObjectProcedure<PruningVertex>() {
                @Override
                public void apply(PruningVertex vertex) {
                    int w = vertex.vertex;
                    double weight = vertex.weight;

                    IntObjectMap<BucketEntry> wBuckets = buckets.get(w);
                    if (wBuckets != null) {
                        wBuckets.forEach(new IntObjectProcedure<BucketEntry>() {
                            @Override
                            public void apply(int i, BucketEntry wBucket) {
                                if (wBucket.weight > weight + weights[wBucket.idx]) {
                                    wBuckets.remove(i);
                                }
                            }
                        });
                    }
                }
            });
        }
    }

    private void discoverBucketEntriesToCopy(RankedNode baseNode, IntObjectMap<ObjectArrayList<SBIEntry>> vertices,
                                             double[] distances, long[] times, double[] weights,
                                             IntObjectMap<IntObjectMap<BucketEntry>> buckets,
                                             IntIntHashMap terminals) {

        System.out.println("     Discover Buckets: ID: " + baseNode.adjNode);
        ObjectArrayList<SBIEntry> downList = vertices.get(baseNode.adjNode);
        if (downList != null) {
            downList.forEach(new ObjectProcedure<SBIEntry>() {
                @Override
                public void apply(SBIEntry vertex) {
                    int w = vertex.adj;
                    double weight = vertex.weight;

                    System.out.println("     W: " + w);

                    IntObjectMap<BucketEntry> bucketTargets = buckets.get(w);
                    if (bucketTargets != null) {

                        bucketTargets.forEach(new IntObjectProcedure<BucketEntry>() {
                            @Override
                            public void apply(int t, BucketEntry entry) {

                                int targetIdx = entry.idx;
                                double pathWeight = entry.weight;

                                terminals.putIfAbsent(t, targetIdx);

                                double currentWeight = pathWeight + weight;
                                if (currentWeight < weights[targetIdx] || weights[targetIdx] == 0) {
                                      System.out.println("Bucket: " + baseNode + " -> " + vertex.adj + "(" + w
                                              + ") -> " + entry.distance + " - "  + vertex.distance);
                                    weights[targetIdx] = currentWeight;
                                    distances[targetIdx] = entry.distance + vertex.distance;
                                    times[targetIdx] = entry.time + vertex.time;
                                }
                            }
                        });
                    }
                }
            });
        }
    }

    private void copyBucketsEntriesToCurrentVertex(RankedNode current,
                                                   IntObjectMap<IntObjectMap<BucketEntry>> buckets,
                                                   double[] distances, long[] times, double[] weights, boolean saveBestPath,
                                                   DistanceMatrix dm, IntObjectMap<IntArrayList> targets,
                                                   IntIntHashMap terminals) {
        terminals.forEach(new IntIntProcedure() {
            @Override
            public void apply(int source, int sourceIdx) {
                IntObjectMap<BucketEntry> b = buckets.get(current.adjNode);

                if (b == null) {
                    b = new IntObjectHashMap<>();
                    buckets.put(current.adjNode, b);
                }

                BucketEntry entry = new BucketEntry(weights[sourceIdx], times[sourceIdx], distances[sourceIdx], sourceIdx);
                b.put(source, entry);
                System.out.println("Copy Bucket: " + current.edge + " --> " + current.adjNode  + "(" + current.traversalId + ")");
                System.out.println("Copy Bucket:" + current.adjNode + " - " +  source + " -> " + entry);
                if (saveBestPath) {
                    saveBestPath(source, sourceIdx, current.adjNode, entry, targets, dm);
                }

                weights[sourceIdx] = Double.POSITIVE_INFINITY;
                distances[sourceIdx] = Double.POSITIVE_INFINITY;
                times[sourceIdx] = Long.MAX_VALUE;
            }
        });
    }

    private void simultaneousBucketInitializationBackward(List<Snap> targets) {

        IntIntHashMap terminals = new IntIntHashMap();

        double[] distances = new double[targets.size()];
        long[] times = new long[targets.size()];
        double[] weights = new double[targets.size()];

        System.out.println("##################################");
        System.out.println("BACKWARD");

        while (!heap.isEmpty()) {

            RankedNode rn = heap.poll();

            System.out.println(" >>>>>> Processing Backward: " + rn.adjNode + " id: " + rn.traversalId);
            System.out.println(" >>>>>>");

            initializeVertices(rn, inEdgeExplorerNoVirtual, inVertices, true);

            //initializePruningVertices(rn, outEdgeExplorerNoVirtual, false);

            discoverBucketEntriesToCopy(rn, inVertices, distances, times, weights, backwardBuckets, terminals);

            //applyRetrospectivePruningAlgorithm(rn, weights, backwardBuckets);

            copyBucketsEntriesToCurrentVertex(rn, backwardBuckets, distances, times, weights, false, null, null, terminals);

            terminals.clear();

        }


        System.out.println("??####################### " + backwardBuckets.get(8187));
    }

    private void simultaneousBucketInitializationForward(List<Snap> sources,
                                                         DistanceMatrix dm, IntObjectMap<IntArrayList> targets) {

        System.out.println("##################################");
        System.out.println("FORWARD");
        double[] distances = new double[sources.size()];
        long[] times = new long[sources.size()];
        double[] weights = new double[sources.size()];

        while (!heap.isEmpty()) {

            RankedNode rn = heap.poll();
            System.out.println(" Processing Forward: " + rn.adjNode + " id: " + rn.traversalId);

            initializeVertices(rn, outEdgeExplorerNoVirtual, outVertices, false);

            //initializePruningVertices(rn, inEdgeExplorerNoVirtual, true);


            discoverBucketEntriesToCopy(rn, outVertices, distances, times, weights, forwardBuckets, terminals);

            //applyRetrospectivePruningAlgorithm(rn, weights, forwardBuckets);

            copyBucketsEntriesToCurrentVertex(rn, forwardBuckets, distances, times, weights,
                    true, dm, targets, terminals);

            terminals.clear();
        }
    }

    private void saveBestPath(int sourceNode, int idxSource, int currNode, BucketEntry forwardEntry,
                              IntObjectMap<IntArrayList> targets, DistanceMatrix dm) {

        final IntObjectMap<BucketEntry> backwardEntries = backwardBuckets.get(currNode);
        System.out.println(" Saving: " + currNode + " - " + backwardEntries);

        if (backwardEntries != null) {

            backwardEntries.forEach(new IntObjectProcedure<BucketEntry>() {
                @Override
                public void apply(int target, BucketEntry entry) {

                    if (sourceNode != target) {

                        long uniqueId = PairingUtils.pair(sourceNode, target);

                        final double savedWeight = tentativeWeights.get(uniqueId);
                        final double currentWeight = forwardEntry.weight + entry.weight;

                        if (savedWeight == 0.0 || (currentWeight < savedWeight)) {

                            System.out.println(" Current:" + currentWeight + " Saved: " + savedWeight);

                            final long time = forwardEntry.time + entry.time;
                            final double distance = forwardEntry.distance + entry.distance;
                            tentativeWeights.put(uniqueId, currentWeight);

                            IntArrayList targetsIdxs = targets.get(target);
                            int[] buffer = targetsIdxs.buffer;
                            int size = targetsIdxs.size();

                            for (int i = 0; i < size; i++) {
                                System.out.println("Saved Distance: " + distance);
                                dm.setCell(idxSource, buffer[i], distance, time);
                            }
                        }
                    }
                }
            });
        }
    }
}