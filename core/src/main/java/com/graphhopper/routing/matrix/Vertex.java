package com.graphhopper.routing.matrix;

import java.util.Objects;

public class Vertex {
    public int baseNode;
    public int adjNode;
    public int edge;
    public int origEdgeId;
    public double weight;
    public long time;
    public double distance;
    public int origEdgeFirst;
    public int origEdgeLast;


    public Vertex(int base,int edge, int adj,int origEdgeId,double weight, long time, double distance, int origEdgeFirst, int origEdgeLast) {
        this.baseNode = base;
        this.edge = edge;
        this.adjNode = adj;
        this.origEdgeId = origEdgeId;
        this.weight = weight;
        this.time = time;
        this.distance = distance;
        this.origEdgeFirst = origEdgeFirst;
        this.origEdgeLast = origEdgeLast;
    }

    public boolean isSelfLoop(){
        return this.baseNode == this.adjNode;
    }

    public Vertex withSelf(Vertex self){
        if(!self.isSelfLoop()){
            throw new IllegalStateException("Trying to add not self vertex");
        }

        return new Vertex(this.baseNode,this.edge,this.adjNode,this.origEdgeId,
                this.weight + self.weight, this.time + self.time, this.distance + self.distance, this.origEdgeFirst, this.origEdgeLast);
    }

    @Override
    public String toString() {
        return "Vertex{" +
                "base=" + baseNode +
                " edge=" + edge +
                " adj =" + adjNode +
                " origEdgeId=" + origEdgeId +
                ", weight=" + weight +
                ", time=" + time +
                ", distance=" + distance +
                ", origFirst=" + origEdgeFirst +
                ", origLast=" + origEdgeLast +
                '}';
    }

    @Override
    public boolean equals(Object o) {
        if (this == o) return true;
        if (o == null || getClass() != o.getClass()) return false;
        Vertex vertex = (Vertex) o;
        return baseNode == vertex.baseNode && edge == vertex.edge && Double.compare(vertex.weight, weight) == 0;
    }

    @Override
    public int hashCode() {
        return Objects.hash(baseNode, edge, weight);
    }
}
