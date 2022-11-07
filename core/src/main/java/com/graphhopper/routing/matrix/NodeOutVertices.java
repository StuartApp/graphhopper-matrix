package com.graphhopper.routing.matrix;

import com.carrotsearch.hppc.LongObjectHashMap;


public class NodeOutVertices {


    private LongObjectHashMap<Vertex> outs;
    private LongObjectHashMap<Vertex> selfs;

    private Vertex[] outsBuffer;
    private Vertex[] selfBuffer;

    public NodeOutVertices(){
        this.outs = new LongObjectHashMap<>();
        this.selfs = new LongObjectHashMap<>();
    }

    public void add(Vertex v){

        if(v.isSelfLoop()){
            Vertex self = selfs.get(v.origEdgeId);
            if(self == null || self.weight > v.weight){
                selfs.put(v.origEdgeId,v);
            }
        }else{
            Vertex out = outs.get(v.origEdgeId);
            if(out == null || out.weight > v.weight){
                outs.put(v.origEdgeId,v);
            }
        }
    }

    public boolean hasNoSelfLoops(){
        return this.selfs.isEmpty();
    }

    public boolean hasSelfLoops(){
        return !this.selfs.isEmpty();
    }

    public boolean hasOutValues(){
        return !this.outs.isEmpty();
    }

    public boolean hasNotOutValues(){
        return this.outs.isEmpty();
    }

    public Vertex[] outValues(){
        if(hasOutValues()){
            if(outsBuffer == null){
                this.outsBuffer = this.outs.values().toArray(Vertex.class);
            }
            return this.outsBuffer;
        }

        throw new IllegalStateException("Out Vertices without values");
    }

    public int outsSize(){
        return this.outs.size();
    }

    public Vertex[] selfValues(){
        if(hasSelfLoops()){
            if(selfBuffer == null){
                this.selfBuffer = this.selfs.values().toArray(Vertex.class);
            }
            return this.selfBuffer;
        }

        throw new IllegalStateException("Self Vertices without values");
    }

    public int selfSize(){
        return this.selfs.size();
    }

}
