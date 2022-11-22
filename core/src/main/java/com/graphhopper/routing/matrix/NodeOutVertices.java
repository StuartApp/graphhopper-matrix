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
            Vertex self = selfs.get(v.edge);
            if(self == null || self.weight > v.weight){
                selfs.put(v.edge,v);
            }
        }else{


            Vertex out = outs.get(v.edge);

            if(v.baseNode == 2966296 && v.adjNode == 2318131 && out != null){
                System.out.println(v.edge + " - current:" + out.weight + " vertex: " + v.weight);
            }else if(v.baseNode == 2966296 && v.adjNode == 2318131 && out == null){
                    System.out.println(v.edge + " - current:" + 0 + " vertex: " + v.weight);

            }

            if(out == null || out.weight > v.weight){

                outs.put(v.edge,v);
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
