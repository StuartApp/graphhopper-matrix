package com.graphhopper.routing.matrix;

import java.util.Objects;

public class DownEdge {
    public int baseNode;
    public double weight;
    public long time;
    public double distance;


    public DownEdge(int base, double weight, long time, double distance) {
        this.baseNode = base;
        this.weight = weight;
        this.time = time;
        this.distance = distance;
    }

    @Override
    public String toString() {
        return "VertexV2{" +
                "baseNode=" + baseNode +
                ", weight=" + weight +
                ", time=" + time +
                ", distance=" + distance +
                '}';
    }

    @Override
    public boolean equals(Object o) {
        if (this == o) return true;
        if (o == null || getClass() != o.getClass()) return false;
        DownEdge vertexV2 = (DownEdge) o;
        return baseNode == vertexV2.baseNode && Double.compare(vertexV2.weight, weight) == 0;
    }

    @Override
    public int hashCode() {
        return Objects.hash(baseNode, weight);
    }
}