package com.graphhopper.routing.matrix.algorithm;

import com.carrotsearch.hppc.*;
import com.carrotsearch.hppc.procedures.IntIntProcedure;
import com.carrotsearch.hppc.procedures.IntObjectProcedure;
import com.carrotsearch.hppc.procedures.IntProcedure;
import com.carrotsearch.hppc.procedures.ObjectProcedure;
import com.graphhopper.coll.GHIntObjectHashMap;
import com.graphhopper.routing.matrix.*;
import com.graphhopper.routing.querygraph.QueryRoutingCHGraph;
import com.graphhopper.routing.subnetwork.SubnetworkStorage;
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

    IntObjectMap<ObjectArrayList<Vertex>> inVertices;
    IntObjectMap<ObjectArrayList<PruningVertex>> prunningVertices;
    IntObjectMap<ObjectArrayList<Vertex>> outVertices;

    IntSet nodesAdded;
    IntIntHashMap terminals;

    public AbstractManyToMany(QueryRoutingCHGraph graph, RoutingCHGraph graphNoVirtualNodes) {

        this.graph = graph;
        this.weighting = graph.getWrappedWeighting();
        this.inEdgeExplorer = graph.createInEdgeExplorer();
        this.outEdgeExplorer = graph.createOutEdgeExplorer();
        this.maxNodes = graph.getBaseGraph().getBaseGraph().getNodes();

        this.graphNoVirtualNodes = graphNoVirtualNodes;
        this.inEdgeExplorerNoVirtual = graphNoVirtualNodes.createInEdgeExplorer();
        this.outEdgeExplorerNoVirtual = graphNoVirtualNodes.createOutEdgeExplorer();

        this.nodesAdded = new IntHashSet();
        this.terminals = new IntIntHashMap();

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

    protected int getOrigEdgeId(RoutingCHEdgeIteratorState edge, boolean reverse) {
        return reverse ? edge.getOrigEdgeFirst() : edge.getOrigEdgeLast();
    }

    protected abstract double calcWeight(RoutingCHEdgeIteratorState iter, MatrixEntry currEdge, boolean reverse);

    protected abstract long calcTime(RoutingCHEdgeIteratorState iter, MatrixEntry currEdge, boolean reverse);

    protected int getIncomingEdge(MatrixEntry entry) {
        return entry.edge;
    }

    protected abstract double calcDistance(RoutingCHEdgeIteratorState iter, MatrixEntry currEdge);


    private double calcWeight(RoutingCHEdgeIteratorState edgeState, Boolean reverse, int prevOrNextEdgeId) {

        double edgeWeight = edgeState.getWeight(reverse);
        final int origEdgeId = reverse ? edgeState.getOrigEdgeLast() : edgeState.getOrigEdgeFirst();
        double turnCosts = reverse
                ? graph.getTurnWeight(origEdgeId, edgeState.getBaseNode(), prevOrNextEdgeId)
                : graph.getTurnWeight(prevOrNextEdgeId, edgeState.getBaseNode(), origEdgeId);
        return edgeWeight + turnCosts;

    }


    private long calcTime(RoutingCHEdgeIteratorState edgeState, Boolean reverse, int prevOrNextEdgeId) {

        long time = edgeState.getTime(reverse);
        long turnCost;
        int origEdgeId;
        if (reverse) {
            origEdgeId = edgeState.getOrigEdgeLast();
            turnCost = weighting.calcTurnMillis(origEdgeId, edgeState.getBaseNode(), prevOrNextEdgeId);
        } else {
            origEdgeId = edgeState.getOrigEdgeFirst();
            turnCost = weighting.calcTurnMillis(prevOrNextEdgeId, edgeState.getBaseNode(), origEdgeId);
        }

        return time + turnCost;
    }


    protected double calcDistance(RoutingCHEdgeIteratorState iter) {
        return iter.getDistance();
    }

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
        //System.out.println("Backward Init **********************");

        int closestNode = snap.getClosestNode();
        IntSet processed = new IntHashSet();

        //During the initialization, we add pair (t, 0) to B(t) for each t ∈ T. By adding
        //such a pair, we indicate that t can be reached from t, with a shortest path of length zero.
        IntObjectMap<BucketEntry> bucketTargets = new GHIntObjectHashMap<>();
        bucketTargets.put(closestNode, new BucketEntry(0, 0, 0, idx));
        backwardBuckets.put(closestNode, bucketTargets);

        //Closest Node is Virtual
        if (isVirtual(closestNode)) {

            Deque<VirtualNodeEntry> queue = new ArrayDeque<>();
            queue.add(new VirtualNodeEntry(closestNode, 0, 0, 0));

            while (!queue.isEmpty()) {

                VirtualNodeEntry current = queue.poll();
                int baseNode = current.node;


                RoutingCHEdgeIterator downIterator = inEdgeExplorer.setBaseNode(baseNode);
                while (downIterator.next()) {

                    int adjNode = downIterator.getAdjNode();
                    //System.out.println("Base:" + baseNode + " --> " + adjNode);

                    boolean isVirtualAdj = isVirtual(adjNode);


                    if (isVirtualAdj && adjNode != closestNode && !processed.contains(adjNode)) {
                        processed.add(adjNode);
                        double weight = calcWeight(downIterator, true, baseNode) + current.weight;
                        double distance = calcDistance(downIterator) + current.distance;
                        long time = calcTime(downIterator, true, baseNode) + current.time;
                        queue.add(new VirtualNodeEntry(adjNode, weight, time, distance));
                    } else if (adjNode != closestNode && !processed.contains(adjNode)) {
                        heap.add(new RankedNode(adjNode, graph.getLevel(adjNode)));
                        nodesAdded.add(adjNode);

                        //Down Vertices
                        ObjectArrayList<Vertex> dVertices = inVertices.get(adjNode);
                        if (dVertices == null) {
                            dVertices = new ObjectArrayList<>();
                            inVertices.put(adjNode, dVertices);
                        }

                        double weight = calcWeight(downIterator, true, baseNode) + current.weight;

                        if (weight < Double.POSITIVE_INFINITY) {
                            double distance = calcDistance(downIterator) + current.distance;
                            long time = calcTime(downIterator, true, baseNode) + current.time;
                            dVertices.add(new Vertex(closestNode, adjNode, weight, time, distance));
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
                heap.add(new RankedNode(uniqueIn.getAdjNode(), graph.getLevel(uniqueIn.getAdjNode())));
                nodesAdded.add(uniqueIn.getAdjNode());

                double weight = calcWeight(uniqueIn, true, closestNode);

                //Down Vertices
                ObjectArrayList<Vertex> dVertices = inVertices.get(uniqueIn.getAdjNode());
                if (dVertices == null) {
                    dVertices = new ObjectArrayList<>();
                    inVertices.put(uniqueIn.getAdjNode(), dVertices);
                }

                if (weight < Double.POSITIVE_INFINITY) {
                    double distance = calcDistance(uniqueIn);
                    long time = calcTime(uniqueIn, true, closestNode);
                    dVertices.add(new Vertex(closestNode,uniqueIn.getAdjNode(), weight, time, distance));
                }
            } else {
                heap.add(new RankedNode(closestNode, graph.getLevel(closestNode)));
                nodesAdded.add(closestNode);
            }
        }
    }

    private void findInitialNodesForward(Snap snap, int idx,
                                         IntObjectMap<IntArrayList> targets, DistanceMatrix dm) {

        int closestNode = snap.getClosestNode();
        IntSet processed = new IntHashSet();

        //During the initialization, we add pair (t, 0) to B(t) for each t ∈ T. By adding
        //such a pair, we indicate that t can be reached from t, with a shortest path of length zero.
        IntObjectMap<BucketEntry> bucketTargets = new GHIntObjectHashMap<>();
        bucketTargets.put(closestNode, new BucketEntry(0, 0, 0, idx));
        forwardBuckets.put(closestNode, bucketTargets);

        Deque<VirtualNodeEntry> queue = new ArrayDeque<>();
        queue.add(new VirtualNodeEntry(closestNode, 0, 0, 0));

        //Process Queue

        while (!queue.isEmpty()) {

            VirtualNodeEntry current = queue.poll();
            int baseNode = current.node;

            RoutingCHEdgeIterator outIterator = outEdgeExplorer.setBaseNode(baseNode);
            while (outIterator.next()) {

                int adjNode = outIterator.getAdjNode();

                boolean isVirtualAdj = isVirtual(adjNode);

                double weight = calcWeight(outIterator, false, baseNode) + current.weight;

                if (weight < Double.POSITIVE_INFINITY) {


                    if (targets.containsKey(adjNode)) {
                        long uniqueId = PairingUtils.pair(closestNode, adjNode);

                        final double savedWeight = tentativeWeights.get(uniqueId);

                        if ((savedWeight == 0.0 || (weight < savedWeight)) && closestNode != adjNode) {

                            double distance = calcDistance(outIterator) + current.distance;
                            long time = calcTime(outIterator, false, baseNode) + current.time;

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
                        long time = calcTime(outIterator, false, baseNode) + current.time;
                        processed.add(adjNode);
                        queue.add(new VirtualNodeEntry(adjNode, weight, time, distance));
                    } else if (adjNode != closestNode && !processed.contains(adjNode)) {
                        this.heap.add(new RankedNode(adjNode, graph.getLevel(adjNode)));
                        this.nodesAdded.add(adjNode);

                        //Out Vertices
                        ObjectArrayList<Vertex> nodeOutVertices = outVertices.get(adjNode);
                        if (nodeOutVertices == null) {
                            nodeOutVertices = new ObjectArrayList<>();
                            outVertices.put(adjNode, nodeOutVertices);
                        }

                        double distance = calcDistance(outIterator) + current.distance;
                        long time = calcTime(outIterator, false, baseNode) + current.time;

                        nodeOutVertices.add(new Vertex(closestNode, adjNode, weight, time, distance));
                    }
                }
            }
        }
    }

    private void initializeVertices(int baseNode, RoutingCHEdgeExplorer explorer,
                                    IntObjectMap<ObjectArrayList<Vertex>> vertices, boolean reverse) {

        RoutingCHEdgeIterator downIterator = explorer.setBaseNode(baseNode);
        IntObjectMap<Vertex> vertexs = new IntObjectHashMap<>();

        while (downIterator.next()) {

            int adjNode = downIterator.getAdjNode();

                int adjRank = graph.getLevel(adjNode);
                boolean accept = this.sbiLevelEdgeFilter.accept(downIterator);

                if ((accept && !nodesAdded.contains(adjNode))) {
                    nodesAdded.add(adjNode);
                    heap.add(new RankedNode(adjNode, adjRank));
                }

                if (accept) {

                    double weight = calcWeight(downIterator, reverse, baseNode);

                    if (weight < Double.POSITIVE_INFINITY) {

                        int originId = getOrigEdgeId(downIterator,reverse);
                        int traversalId = getTraversalId(downIterator, originId,reverse);

                        Vertex storedVertex = vertexs.get(traversalId);
                        if(storedVertex == null || storedVertex.weight > weight){
                            double distance = calcDistance(downIterator);
                            long time = calcTime(downIterator, reverse, baseNode);
                            vertexs.put(traversalId,new Vertex(baseNode, adjNode, weight, time, distance));
                        }
                    }
                }
        }

        //Because we saw that from a node can exists multiples edges to adj node, we take into consideration
        //only the edges with less weight
        vertexs.forEach(new IntObjectProcedure<Vertex>() {
            @Override
            public void apply(int traversalId, Vertex vertex) {
                ObjectArrayList<Vertex> dVertices = vertices.get(vertex.adj);
                if (dVertices == null) {
                    dVertices = new ObjectArrayList<>();
                    vertices.put(vertex.adj, dVertices);
                }
                dVertices.add(vertex);
            }
        });

    }

    private void initializePruningVertices(int baseNode, RoutingCHEdgeExplorer explorer, boolean reverse) {

        RoutingCHEdgeIterator upIterator = explorer.setBaseNode(baseNode);
        while (upIterator.next()) {
            int adjNode = upIterator.getAdjNode();
            boolean accept = this.sbiLevelEdgeFilter.accept(upIterator);

            if (accept) {
                ObjectArrayList<PruningVertex> uVertices = prunningVertices.get(adjNode);
                if (uVertices == null) {
                    uVertices = new ObjectArrayList<>();
                    prunningVertices.put(adjNode, uVertices);
                }

                double weight = calcWeight(upIterator, reverse, baseNode);
                if (weight < Double.POSITIVE_INFINITY) {
                    uVertices.add(new PruningVertex(baseNode, weight));
                }
            }
        }
    }

    private void applyRetrospectivePruningAlgorithm(int baseNode, double[] weights,
                                                    IntObjectMap<IntObjectMap<BucketEntry>> buckets) {

        ObjectArrayList<PruningVertex> upList = prunningVertices.get(baseNode);

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

    private void discoverBucketEntriesToCopy(int baseNode, IntObjectMap<ObjectArrayList<Vertex>> vertices,
                                             double[] distances, long[] times, double[] weights,
                                             IntObjectMap<IntObjectMap<BucketEntry>> buckets,
                                             IntIntHashMap terminals) {


        ObjectArrayList<Vertex> downList = vertices.get(baseNode);
        if (downList != null) {
            downList.forEach(new ObjectProcedure<Vertex>() {
                @Override
                public void apply(Vertex vertex) {
                    int w = vertex.base;
                    double weight = vertex.weight;

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
                                    //if(find.contains(w)) {
                                      System.out.println("Bucket: " + baseNode + " -> " + w + " -> " + entry.distance + " - "  + vertex.distance);
                                    //}
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

    private void copyBucketsEntriesToCurrentVertex(int baseNode,
                                                   IntObjectMap<IntObjectMap<BucketEntry>> buckets,
                                                   double[] distances, long[] times, double[] weights, boolean saveBestPath,
                                                   DistanceMatrix dm, IntObjectMap<IntArrayList> targets,
                                                   IntIntHashMap terminals) {
        terminals.forEach(new IntIntProcedure() {
            @Override
            public void apply(int source, int sourceIdx) {
                IntObjectMap<BucketEntry> b = buckets.get(baseNode);

                if (b == null) {
                    b = new IntObjectHashMap<>();
                    buckets.put(baseNode, b);
                }

                BucketEntry entry = new BucketEntry(weights[sourceIdx], times[sourceIdx], distances[sourceIdx], sourceIdx);
                b.put(source, entry);
                System.out.println("Saving Bucket:" + baseNode + " - " +  source + " -> " + entry);
                if (saveBestPath) {
                    saveBestPath(source, sourceIdx, baseNode, entry, targets, dm);
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
            int baseNode = rn.node;
            System.out.println(" Processing Backward " + baseNode);

            initializeVertices(baseNode, inEdgeExplorerNoVirtual, inVertices, true);

            //initializePruningVertices(baseNode, outEdgeExplorerNoVirtual, false);

            discoverBucketEntriesToCopy(baseNode, inVertices, distances, times, weights, backwardBuckets, terminals);

            //applyRetrospectivePruningAlgorithm(baseNode, weights, backwardBuckets);

            copyBucketsEntriesToCurrentVertex(baseNode, backwardBuckets, distances, times, weights, false, null, null, terminals);

            if(baseNode == 8187){
                System.out.println("####################### " + backwardBuckets.get(8187));
            }

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
            int baseNode = rn.node;
            System.out.println("Processing Forward: " + baseNode);

            initializeVertices(baseNode, outEdgeExplorerNoVirtual, outVertices, false);

            initializePruningVertices(baseNode, inEdgeExplorerNoVirtual, true);


            discoverBucketEntriesToCopy(baseNode, outVertices, distances, times, weights, forwardBuckets, terminals);

            applyRetrospectivePruningAlgorithm(baseNode, weights, forwardBuckets);

            copyBucketsEntriesToCurrentVertex(baseNode, forwardBuckets, distances, times, weights,
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