package com.graphhopper.routing.matrix;

public class VertexTerminal {
    public int baseNode;
    public int adjNode;
    public int edge;
    public int origEdgeId;
    public double weight;
    public long time;
    public double distance;
    public int terminal;
    public int terminalIdx;
    public int origEdgeFirst;
    public int origEdgeLast;

    public double bucketWeight;
    public long bucketTime;
    public double bucketDistance;
    public double vertexWeight;
    public long vertexTime;
    public double vertexDistance;

    public VertexTerminal(Vertex vert, int terminal, BucketEntry entry) {
        this.baseNode = vert.baseNode;
        this.edge = vert.edge;
        this.adjNode = vert.adjNode;
        this.origEdgeId = vert.origEdgeId;
        this.weight = vert.weight + entry.weight;
        this.time = vert.time + entry.time;
        this.distance = vert.distance + entry.distance;
        this.terminal = terminal;
        this.terminalIdx = entry.idx;
        this.origEdgeFirst = vert.origEdgeFirst;
        this.origEdgeLast = vert.origEdgeLast;

        this.bucketWeight = entry.weight;
        this.bucketTime = entry.time;
        this.bucketDistance = entry.distance;
        this.vertexWeight = vert.weight;
        this.vertexTime = vert.time;
        this.vertexDistance = vert.distance;
    }

    public boolean isSelfLoop(){
        return this.baseNode == this.adjNode;
    }


    @Override
    public String toString() {
        return "VertexTerminal{" +
                " baseNode=" + baseNode +
                ", adjNode=" + adjNode +
                ", edge=" + edge +
                ", origEdgeId=" + origEdgeId +
                ", weight=" + weight +
                ", time=" + time +
                ", distance=" + distance +
                ", terminal=" + terminal +
                ", entry w=" + bucketWeight +
                ", vert w=" + vertexWeight +
                '}';
    }
}
