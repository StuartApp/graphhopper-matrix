package com.graphhopper.routing.matrix;

import com.carrotsearch.hppc.IntHashSet;
import com.carrotsearch.hppc.IntSet;
import com.carrotsearch.hppc.LongObjectHashMap;
import com.graphhopper.util.PairingUtils;

public class NodeInVertices {

    private IntSet origEdgeIds;

    private LongObjectHashMap<VertexTerminal> inVertices;

    private LongObjectHashMap<VertexTerminal> minsVertices;

    private VertexTerminal[] ins;
    private VertexTerminal[] mins;

    public NodeInVertices(){
        this.origEdgeIds = new IntHashSet();
        this.inVertices = new LongObjectHashMap<>();
        this.minsVertices = new LongObjectHashMap<>();
    }

    public void add(Vertex in, BucketEntry entry, int i){

        if(ins == null){
            Long idx = PairingUtils.pair(i,in.origEdgeId);
            VertexTerminal vt = new VertexTerminal(in,i,entry);
            VertexTerminal current = inVertices.get(idx);
            if(current == null || current.weight > vt.weight){
                inVertices.put(idx,vt);
                addMin(vt);
            }
        }else{
            throw new IllegalStateException("You cannot add more vertex once ins are generated");
        }

    }

    private void addMin(VertexTerminal in){

        int terminalNode = in.terminal;
        VertexTerminal min = minsVertices.get(terminalNode);
        if(min == null || in.weight < min.weight){
            minsVertices.put(terminalNode,in);
        }
    }

    public boolean containsEdge(int edgeId){
        return origEdgeIds.contains(edgeId);
    }

    public void addOrigEdgeIds(int origId){
        this.origEdgeIds.add(origId);
    }

    public boolean isNotEmpty(){
        return !this.inVertices.isEmpty();
    }

    public boolean isEmpty(){
        return this.inVertices.isEmpty();
    }

    public VertexTerminal[] values(){
        if(isNotEmpty()){

            if(ins == null){
                this.ins = this.inVertices.values().toArray(VertexTerminal.class);
            }
            return ins;
        }

        throw new IllegalStateException("VertexTerminal without values");
    }

    public int size(){
            return this.inVertices.size();
    }

    public VertexTerminal[] minValues(){
        if(!this.minsVertices.isEmpty()){
            if(this.mins == null){
                this.mins = this.minsVertices.values().toArray(VertexTerminal.class);
            }
            return this.mins;
        }

        throw new IllegalStateException("VertexTerminal without min values");
    }

    public int minSize(){
            return this.minsVertices.size();
    }

}