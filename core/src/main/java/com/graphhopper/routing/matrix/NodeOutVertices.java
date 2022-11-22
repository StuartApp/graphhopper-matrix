package com.graphhopper.routing.matrix;

import com.carrotsearch.hppc.ObjectArrayList;


public class NodeOutVertices {

    private ObjectArrayList<Vertex> outs;
    private ObjectArrayList<Vertex> selfs;

    private Vertex[] outsBuffer;
    private Vertex[] selfBuffer;

    public NodeOutVertices(){
        this.outs = new ObjectArrayList<>();
        this.selfs = new ObjectArrayList<>();
    }

    public void add(Vertex v){

        if(v.isSelfLoop()){
            selfs.add(v);
        }else{
            outs.add(v);
        }
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
                this.outsBuffer = this.outs.toArray(Vertex.class);
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
                this.selfBuffer = this.selfs.toArray(Vertex.class);
            }
            return this.selfBuffer;
        }

        throw new IllegalStateException("Self Vertices without values");
    }

    public int selfSize(){
        return this.selfs.size();
    }

}
