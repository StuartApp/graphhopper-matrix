package com.graphhopper.routing.matrix;

import com.carrotsearch.hppc.*;
import com.carrotsearch.hppc.procedures.IntProcedure;
import com.graphhopper.storage.RoutingCHGraph;
import com.graphhopper.util.PairingUtils;
import java.util.PriorityQueue;

public class SBIAlgorithm {

    RoutingCHGraph graph;
    PriorityQueue<RankedNode> heap;

    IntObjectMap<IntArrayList> targetIndexesNodes;
    IntObjectMap<IntArrayList> sourcesIndexesNodes;

    IntObjectMap<NodeTerminals> nodeTerminals;

    IntSet traversed;

    DistanceMatrix dm;

    LongDoubleMap shortestRoutes;


    //IntDoubleMap forwardComparator = new IntDoubleHashMap();
    //IntIntMap forwardComparatorOrder = new IntIntHashMap();
    //IntDoubleMap backwardComparator = new IntDoubleHashMap();
    //IntIntMap backwardComparatorOrder = new IntIntHashMap();


    public SBIAlgorithm(RoutingCHGraph graph, PriorityQueue<RankedNode> heap, DistanceMatrix dm ) {

        this.graph = graph;
        this.heap = heap;
        this.targetIndexesNodes = new IntObjectHashMap<>();
        this.sourcesIndexesNodes = new IntObjectHashMap<>();
        this.nodeTerminals = new IntObjectHashMap<>();
        this.dm = dm;
        this.shortestRoutes = new LongDoubleHashMap();
        this.traversed = new IntHashSet();

        /*
        this.forwardComparator.put(2966296,8010.0);
        this.forwardComparator.put(2966296,69626.0);
        this.forwardComparator.put(150209,18021.0);
        this.forwardComparator.put(490,10150.0);

        this.forwardComparatorOrder.put(2966296,2966296);
        this.forwardComparatorOrder.put(150209,2966296);
        this.forwardComparatorOrder.put(490,150209);


        this.backwardComparator.put(2966296,229917.0);
        this.backwardComparator.put(1612521,162879.0);
        this.backwardComparator.put(65376,63344.0);
        this.backwardComparator.put(284578,28668.0);
        this.backwardComparator.put(284580,2042.0);

        this.backwardComparatorOrder.put(1612521,2966296);
        this.backwardComparatorOrder.put(65376,1612521);
        this.backwardComparatorOrder.put(284578,65376);
        this.backwardComparatorOrder.put(284580,284578);

         */

    }

    public void addSource(int sourceNode, int sourceIdx){

        IntArrayList sourceIdxs = sourcesIndexesNodes.get(sourceNode);
        if(sourceIdxs == null){
            sourceIdxs = new IntArrayList();
            sourcesIndexesNodes.put(sourceNode,sourceIdxs);
            initializeTerminals(sourceNode,sourceIdx,false);
            findRoutes(sourceNode,null);
        }

        sourceIdxs.add(sourceIdx);
    }

    public void addTarget(int targetNode, int targetIdx){

        IntArrayList targetIdxs = targetIndexesNodes.get(targetNode);
        if(targetIdxs == null){
            targetIdxs = new IntArrayList();
            targetIndexesNodes.put(targetNode,targetIdxs);
            initializeTerminals(targetNode,targetIdx,true);
        }

        targetIdxs.add(targetIdx);
    }

    private NodeTerminals obtainNodeTerminals(int node){
        NodeTerminals terminals = nodeTerminals.get(node);
        if(terminals == null){
            terminals = new NodeTerminals();
            nodeTerminals.put(node,terminals);
        }

        return terminals;
    }

    private void initializeTerminals(int terminal, int terminalIdx, boolean reverse){

        NodeTerminals terminals = obtainNodeTerminals(terminal);
        Terminal initial = new Terminal(0,0,0,terminal,terminalIdx);
        terminals.addInitialTerminal(initial,reverse,-1,-1,-1);
    }

    //NO REMOVE - Used for debugging routes errors purposes
    private void checkResult(VertexWithTerminals in,
                             int node, int adjNode, double weight, long time, double distance, boolean reverse){
        /*
        if(reverse){
            if(backwardComparatorOrder.containsKey(node)){
                if(backwardComparatorOrder.get(node) == adjNode){
                    double expected = backwardComparator.get(adjNode);
                    if(expected == time){
                        System.out.println("BACKWARD: " + node + " -> " + adjNode + " : " + time + " = " + expected);
                        in.getTerminalsList().forEach(new IntObjectProcedure<Terminal>() {
                            @Override
                            public void apply(int i, Terminal terminal) {
                                System.out.println(" Target : " + terminal.nodeIdx + " - " + terminal.time);
                            }
                        });
                    }
                }
            }
        }else{
            if(forwardComparatorOrder.containsKey(node)){
                if(forwardComparatorOrder.get(node) == adjNode){
                    double expected = forwardComparator.get(adjNode);
                    if(expected == time){
                        System.out.println("FORWARD: " + in.adjNode + " -> " + adjNode + " : " + time + " = " + expected);
                        in.getTerminalsList().forEach(new IntObjectProcedure<Terminal>() {
                            @Override
                            public void apply(int i, Terminal terminal) {
                                System.out.println(" Source : " + terminal.nodeIdx + " - " + terminal.time);
                            }
                        });
                    }
                }
            }

        }

         */
    }

    public void addInitialOutVertex(int node, int terminal, int terminalIdx, double weight, long time,
                                    double distance, boolean reverse, int origEdgeId, int origEdgeFirst, int origEdgeLast) {

        NodeTerminals terminals = obtainNodeTerminals(node);
        Terminal initial = new Terminal(weight,time,distance,terminal,terminalIdx);
        terminals.addInitialTerminal(initial,reverse,origEdgeId,origEdgeFirst,origEdgeLast);
        addToHeap(node);
        if(reverse) {
            findRoutes(node,null);
        }
    }

    private void addOutVertex(VertexWithTerminals in, Vertex out, boolean reverse){

        //NOT REMOVE!!! checkResult(in, out.baseNode, out.adjNode,out.weight,out.time,out.distance,reverse);

        NodeTerminals terminals = obtainNodeTerminals(out.adjNode);
        terminals.addTerminal(in,out,reverse);
    }

    private VertexWithTerminals[] EMPTY_INS = {};
    private VertexWithTerminals[] getIns( int node, boolean reverse){

        NodeTerminals nt = nodeTerminals.get(node);

        if(nt == null){
            return EMPTY_INS;
        }else{
            return nodeTerminals.get(node).getTerminals(reverse);
        }
    }

    private double calculateTurnCost(int inVertexId, int outVertexId, int node, boolean reverse){

            return reverse
                    ? graph.getTurnWeight(outVertexId, node, inVertexId)
                    : graph.getTurnWeight(inVertexId, node, outVertexId);
    }

    private boolean isAccessible(double turnCost){
        return Double.isFinite(turnCost);
    }

    public void addOuts(int node, NodeOutVertices outs, boolean reverse){

        VertexWithTerminals[] ins = getIns(node,reverse);
        int insSize = ins.length;

        if(outs.hasOutValues()){

            for (int i = 0; i < insSize; i++) {
                VertexWithTerminals in = ins[i];
                int inOrigEdgeId = reverse ? in.origEdgeFirst : in.origEdgeLast;

                for (int ii = 0; ii < outs.outsSize(); ii++) {

                    Vertex out = outs.outValues()[ii];

                    //If an in vertex origEdgeId is equal of this vertex origEdgeId we discard it
                    if(in.origEdgeId == out.origEdgeId){
                        continue;
                    }

                    int outOrigEdgeId = reverse ? out.origEdgeLast : out.origEdgeFirst;

                    double cost = calculateTurnCost(inOrigEdgeId,outOrigEdgeId,node,reverse);


                    //Out is accessible
                    if(isAccessible(cost)){

                        Vertex outWithCost = out.withTurnCost(cost);
                        addOutVertex(in,outWithCost,reverse);
                        addToHeap(out.adjNode);
                    }

                    //Check Selfs
                    if(outs.hasSelfLoops()){

                        for (int iii = 0; iii < outs.selfSize(); iii++) {
                            Vertex self = outs.selfValues()[iii];

                            if(in.origEdgeId == self.origEdgeId){
                                continue;
                            }

                            int selfOrigEdgeId = reverse ? self.origEdgeLast : self.origEdgeFirst;

                            double costInToSelf = calculateTurnCost(inOrigEdgeId,selfOrigEdgeId,node,reverse);

                            //In to self is accessible
                            if(isAccessible(costInToSelf)){

                                double costSelfToOut = calculateTurnCost(selfOrigEdgeId,outOrigEdgeId,node,reverse);
                                //Self to out is accessible
                                if(isAccessible((costSelfToOut))){
                                    Vertex outWithSelfCost = out.withSelf(self).withSelfTurnCost(costInToSelf,costSelfToOut);
                                    addOutVertex(in,outWithSelfCost,reverse);
                                    addToHeap(out.adjNode);
                                }
                            }
                        }

                    }
                }
            }
        }
    }

    private void addToHeap(int node){
        if(traversed.contains(node)){
            RankedNode rankedNode = new RankedNode(node,graph.getLevel(node),false);
            heap.add(rankedNode);
            addTraversed(node);
        }
    }

    public void addTraversed(int node){
        traversed.add(node);
    }

    public boolean containsTraversed(int node){
        return traversed.contains(node);
    }

    public void clearTraversed(){
        traversed.clear();
    }

    private void saveToDistanceMatrix(int sourceIdx, int targetNode, long time, double distance){
        IntArrayList targets = targetIndexesNodes.get(targetNode);
            targets.forEach(new IntProcedure() {
                @Override
                public void apply(int targetIdx) {
                    dm.setCell(sourceIdx,targetIdx, distance, time);
                }
            });
    }

    private void saveShortRoutes(VertexWithTerminals forward, VertexWithTerminals backward, int node){

        int fEdgeId = forward.origEdgeLast;
        int bEdgeId = backward.origEdgeFirst;

        double turnCost = calculateTurnCost(fEdgeId,bEdgeId,node,false);

        if(isAccessible(turnCost)){
            Terminal[] fterminals = forward.getTerminals();
            Terminal[] bterminals = backward.getTerminals();

            int fterminalsSize = fterminals.length;
            int bterminalsSize = bterminals.length;

            for(int i = 0; i < fterminalsSize; i++ ){
                Terminal ft = fterminals[i];

                for(int ii = 0; ii < bterminalsSize; ii++ ){
                    Terminal bt = bterminals[ii];

                    int sourceNode = ft.node;
                    int sourceIdx = ft.nodeIdx;
                    int targetNode = bt.node;

                    long key = PairingUtils.pair(sourceNode, targetNode);
                    double possible = ft.weight + bt.weight + turnCost;
                    double current = shortestRoutes.get(key);

                    if(sourceNode == targetNode){
                        saveToDistanceMatrix(sourceIdx,targetNode,0,0);
                    }else if (current == 0.0 || current > possible){
                        shortestRoutes.put(key,possible);
                        long time = ft.time + bt.time;
                        double distance = ft.distance + bt.distance;
                        saveToDistanceMatrix(sourceIdx,targetNode,time,distance);
                    }
                }
            }
        }

    }

    private void saveShortRoutesWithSelfLoop(VertexWithTerminals forward, VertexWithTerminals backward, int node,
                                             NodeOutVertices outs){

        int fEdgeId = forward.origEdgeLast;
        int bEdgeId = backward.origEdgeFirst;

        Vertex[] selfs = outs.selfValues();
        int selfsSize = outs.selfSize();

        for(int i = 0; i < selfsSize; i++){

            Vertex self = selfs[i];
            double turnCostFromForwardToSelf =  calculateTurnCost(fEdgeId,self.origEdgeFirst,node,false);
            if(isAccessible(turnCostFromForwardToSelf)){
                double turnCostFromSelfToBackward =  calculateTurnCost(self.origEdgeLast,bEdgeId,node,false);
                if(isAccessible(turnCostFromSelfToBackward)){

                    Terminal[] fterminals = forward.getTerminals();
                    Terminal[] bterminals = backward.getTerminals();

                    int fterminalsSize = fterminals.length;
                    int bterminalsSize = bterminals.length;

                    for(int ii = 0; ii < fterminalsSize; ii++ ){
                        Terminal ft = fterminals[ii];

                        for(int iii = 0; iii < bterminalsSize; iii++ ){
                            Terminal bt = bterminals[iii];

                            int sourceNode = ft.node;
                            int sourceIdx = ft.nodeIdx;
                            int targetNode = bt.node;

                            long key = PairingUtils.pair(sourceNode, targetNode);
                            double possible = ft.weight + bt.weight + + self.weight + turnCostFromForwardToSelf + turnCostFromSelfToBackward;
                            double current = shortestRoutes.get(key);

                            if (current == 0.0 || current > possible){
                                shortestRoutes.put(key,possible);
                                long time = ft.time + bt.time + self.time;
                                double distance = ft.distance + bt.distance + self.distance;
                                saveToDistanceMatrix(sourceIdx,targetNode,time,distance);
                            }
                        }
                    }
                }
            }


        }
    }

    public void findRoutes(int node,NodeOutVertices outs){

        NodeTerminals nodeTerminals = obtainNodeTerminals(node);

        if(nodeTerminals.hasPossibleShortRoutes()){

            VertexWithTerminals[] backwardTerminals = nodeTerminals.getBackwardTerminals();
            VertexWithTerminals[] forwardTerminals = nodeTerminals.getForwardTerminals();

            int forwardSize = forwardTerminals.length;
            int backwardSize = backwardTerminals.length;

            for(int i = 0; i < forwardSize; i++){
                VertexWithTerminals f = forwardTerminals[i];

                for(int ii = 0; ii < backwardSize; ii++){
                    VertexWithTerminals b = backwardTerminals[ii];
                    saveShortRoutes(f,b,node);
                    if(outs != null && outs.hasSelfLoops()){
                        saveShortRoutesWithSelfLoop(f,b,node,outs);
                    }
                }
            }
        }
    }

}
