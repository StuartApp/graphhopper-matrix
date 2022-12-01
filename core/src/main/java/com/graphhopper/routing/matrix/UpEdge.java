package com.graphhopper.routing.matrix;

import java.util.Objects;

public class UpEdge {
    public int baseNode;
    public double weight;


    public UpEdge(int base, double weight) {
        this.baseNode = base;
        this.weight = weight;
    }

    @Override
    public String toString() {
        return "UpEdge{" +
                "baseNode=" + baseNode +
                ", weight=" + weight +
                '}';
    }

    @Override
    public boolean equals(Object o) {
        if (this == o) return true;
        if (o == null || getClass() != o.getClass()) return false;
        UpEdge vertexV2 = (UpEdge) o;
        return baseNode == vertexV2.baseNode && Double.compare(vertexV2.weight, weight) == 0;
    }

    @Override
    public int hashCode() {
        return Objects.hash(baseNode, weight);
    }
}