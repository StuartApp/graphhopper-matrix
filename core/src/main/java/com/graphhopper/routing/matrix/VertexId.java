package com.graphhopper.routing.matrix;

import java.util.Objects;

public class VertexId {
    public int edge;
    public int origEdgeId;


    public VertexId(int edge, int origEdgeId) {
        this.edge = edge;
        this.origEdgeId = origEdgeId;
    }
}
